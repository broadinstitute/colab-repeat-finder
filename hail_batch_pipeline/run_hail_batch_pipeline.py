import collections
import json
import logging
import os
from pprint import pformat
import hailtop.fs as hfs

from step_pipeline import pipeline, Backend, Localize, Delocalize

DOCKER_IMAGE = "weisburd/colab-repeat-finder@sha256:ce854aae2fb906de9b8d6e7dc877087f0d6b804de52339c80642a70a81b4be8c"
STR_ANALYSIS_DOCKER_IMAGE = "weisburd/str-analysis@sha256:40eaad32db567bfde2267ef29c5e931abe69cf6eaefd35ff2316a75d08174151"

logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)

DEFAULT_BASE_OUTPUT_DIR = "gs://bw-proj/gnomad-bw/colab-repeat-finder"

#REFERENCE_GENOME_FASTA = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
HG38_REFERENCE_GENOME_FASTA = "gs://str-truth-set/hg38/ref/hg38.fa.gz"
CHM13_REFERENCE_GENOME_FASTA = "gs://str-truth-set/chm13/ref/chm13v2.0.fa.gz"

def run(cmd):
    print(cmd)
    os.system(cmd)


def parse_args(batch_pipeline):
    p = batch_pipeline.get_config_arg_parser()

    p.add_argument("--reference-genome", choices=["hg38", "chm13"], default="hg38")
    p.add_argument("--min-span", type=int, default=9)
    p.add_argument("--min-repeats", type=int, default=3)
    p.add_argument("--min-motif-size", type=int, default=1)
    p.add_argument("--max-motif-size", type=int, default=1_000)  #default=18000//3)
    p.add_argument("--use-nonpreemptibles", action="store_true")
    p.add_argument("--batch-size", help="Interval size in base pairs to process per job", type=int, default=500_000)
    p.add_argument("--output-prefix", help="Output filename prefix")
    p.add_argument("--output-dir", help="Google Cloud Storage output directory")
    p.add_argument("--start-batch-i", type=int, help="Start processing from this batch", default=0)
    p.add_argument("-n", type=int, help="Run only this many batches")

    args = batch_pipeline.parse_known_args()

    return args


def main():

    bp = pipeline("colab repeat finder", backend=Backend.HAIL_BATCH_SERVICE, config_file_path="~/.step_pipeline")

    args = parse_args(bp)
    p = bp.get_config_arg_parser()

    if args.reference_genome == "hg38":
        reference_genome_fasta = HG38_REFERENCE_GENOME_FASTA
    elif args.reference_genome == "chm13":
        reference_genome_fasta = CHM13_REFERENCE_GENOME_FASTA
    else:
        p.error(f"Invalid reference genome: {args.reference_genome}")

    if not args.output_prefix:
        args.output_prefix = f"{args.reference_genome}_repeats"
    if not args.output_dir:
        args.output_dir = os.path.join(DEFAULT_BASE_OUTPUT_DIR, args.reference_genome)
        
    # compute batch sizes
    batches = []
    with hfs.open(f"{reference_genome_fasta}.fai") as fai_index_file:
        for line in fai_index_file:
            fields = line.strip().split("\t")
            chrom = fields[0]
            chrom_size = int(fields[1])
            if len(chrom) > 5:
                # skip the extra contigs
                continue

            for start_0based in range(0, chrom_size, args.batch_size):
                batches.append((chrom, start_0based, min(chrom_size, start_0based + args.batch_size)))

    common_prefix = (
        f"{args.output_prefix}.motifs_{args.min_motif_size}_to_{args.max_motif_size}bp"
        f".repeats_{args.min_repeats}x_and_spans_{args.min_span}bp"
    )
    output_dir = os.path.join(args.output_dir, common_prefix)
    output_paths_glob = os.path.join(output_dir, f"**/*.bed*")
    bp.precache_file_paths(output_paths_glob)
    steps = []
    files_to_download_when_done = []
    for i, (chrom, start_0based, end) in enumerate(batches):
        if i < args.start_batch_i:
            continue
        if args.n and i >= args.start_batch_i + args.n:
            break

        output_prefix = f"{common_prefix}.{chrom}_{start_0based:09d}_{end:09d}"

        s1 = bp.new_step(
            f"{output_prefix} [motifs: {args.min_motif_size}-{args.max_motif_size}]",
            arg_suffix=f"crf",
            image=DOCKER_IMAGE,
            step_number=1,
            cpu=1,
            storage="5Gi",
            memory="standard",
            preemptible=not args.use_nonpreemptibles,
            localize_by=Localize.COPY, #Localize.GSUTIL_COPY, #HAIL_BATCH_CLOUDFUSE,
            delocalize_by=Delocalize.COPY,
            output_dir=output_dir)

        steps.append(s1)

        # process input files
        local_fasta, _, _ = s1.inputs(reference_genome_fasta, f"{reference_genome_fasta}.fai", f"{reference_genome_fasta}.gzi")

        s1.command("set -ex")
        s1.command(f"""time python3 /perfect_repeat_finder.py \
            --min-motif-size {args.min_motif_size} \
            --max-motif-size {args.max_motif_size} \
            --interval {chrom}:{start_0based}-{end} \
            --min-span {args.min_span} \
            --min-repeats {args.min_repeats} \
            --output-prefix {output_prefix} \
            {local_fasta}
        """)

        s1.command(f"bgzip {output_prefix}.bed")
        s1.command(f"tabix {output_prefix}.bed.gz")

        s1.output(f"{output_prefix}.bed.gz")
        s1.output(f"{output_prefix}.bed.gz.tbi")

    concatenated_output_file = f"{common_prefix}.concatenated"
    print(f"Generating combined output file: {concatenated_output_file}")
    s2 = bp.new_step(
        name="combine-step",
        step_number=2,
        image=STR_ANALYSIS_DOCKER_IMAGE,
        storage="20Gi",
        cpu=1,
        localize_by=Localize.COPY,
        delocalize_by=Delocalize.COPY,
        output_dir=output_dir,
    )

    for s1 in steps:
        s2.depends_on(s1)
        
    s2.command("set -ex")
    s2.command(f"gcloud storage cp {os.path.join(output_dir, f'**/{common_prefix}.chr*.bed.gz')} .")
    s2.command(f"ls {common_prefix}.chr*.bed.gz | wc -l")

    s2.command(f"for i in {common_prefix}.chr*.bed.gz; do zcat $i >> {concatenated_output_file}.unsorted.bed; done")
    # sort the concatenated bed file
    s2.command(f"sort -k1,1 -k2,2n {concatenated_output_file}.unsorted.bed | uniq | bgzip > {concatenated_output_file}.bed.gz")
    s2.command(f"tabix {concatenated_output_file}.bed.gz")
    s2.command(f"gunzip -c {concatenated_output_file}.bed.gz | wc -l")
    s2.output(f"{concatenated_output_file}.bed.gz")
    s2.output(f"{concatenated_output_file}.bed.gz.tbi")

    s3 = bp.new_step(
        name="merge-step",
        step_number=3,
        image=STR_ANALYSIS_DOCKER_IMAGE,
        storage="20Gi",
        cpu=1,
        memory="highmem",
        output_dir=output_dir,
    )

    local_concatenated_bed, _ = s3.use_previous_step_outputs_as_inputs(s2)
    s3.command(f"""python3 -u -m str_analysis.merge_loci \
        --verbose \
        --output-format BED \
        --output-prefix {common_prefix} \
        {local_concatenated_bed}
    """)

    s3.command("ls -l")
    s3.command(f"gunzip -c {common_prefix}.bed.gz | wc -l")

    s3.output(f"{common_prefix}.bed.gz")
    s3.output(f"{common_prefix}.bed.gz.tbi")

    files_to_download_when_done.extend([
        (os.path.join(output_dir, f"{common_prefix}.bed.gz"), "results"),
        (os.path.join(output_dir, f"{common_prefix}.bed.gz.tbi"), "results")
    ])
    bp.run()


    # download results
    for remote_path, destination_dir in files_to_download_when_done:
        if not os.path.isdir(destination_dir):
            print(f"Creating local directory: {destination_dir}")
            os.mkdir(destination_dir)
        print(f"Downloading {remote_path} to {destination_dir}/")
        run(f"gsutil -m cp {remote_path} {destination_dir}")


if __name__ == "__main__":
    main()

