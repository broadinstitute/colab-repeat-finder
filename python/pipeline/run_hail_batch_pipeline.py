import collections
import json
import logging
import os
from pprint import pformat
import hailtop.fs as hfs

from step_pipeline import pipeline, Backend, Localize, Delocalize

DOCKER_IMAGE = "weisburd/colab-repeat-finder@sha256:fd956adffeab6f49fbbdbc1930d7e6952a81e1ad450da3252ba1597e6fb6706a"

logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)

REFERENCE_GENOME_FASTA = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"

OUTPUT_BASE_DIR = f"gs://bw-proj/gnomad-bw/colab-repeat-finder/hg38"


def run(cmd):
    print(cmd)
    os.system(cmd)


def parse_args(batch_pipeline):
    p = batch_pipeline.get_config_arg_parser()

    p.add_argument("--min-span", type=int, default=9)
    p.add_argument("--min-repeats", type=int, default=3)
    p.add_argument("--min-motif-size", type=int, default=1)
    p.add_argument("--max-motif-size", type=int, default=1_000) #default=18000//3)
    p.add_argument("--batch-size", help="Interval size in base pairs to process per job", type=int, default=500_000)
    p.add_argument("--output-prefix", default="hg38_repeats")
    p.add_argument("-n", type=int, help="Test-run this many batches")

    args = batch_pipeline.parse_known_args()

    return args


def main():

    bp = pipeline("colab repeat finder", backend=Backend.HAIL_BATCH_SERVICE, config_file_path="~/.step_pipeline")

    args = parse_args(bp)
    p = bp.get_config_arg_parser()

    # compute batch sizes
    batches = []
    with hfs.open(f"{REFERENCE_GENOME_FASTA}.fai") as fai_index_file:
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
    output_dir = os.path.join(OUTPUT_BASE_DIR, common_prefix)
    bp.precache_file_paths(os.path.join(output_dir, f"**/*.bed*"))
    steps = []
    files_to_download_when_done = []
    for i, (chrom, start_0based, end) in enumerate(batches):
        if args.n and i >= args.n:
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
            localize_by=Localize.COPY, #HAIL_BATCH_CLOUDFUSE,
            delocalize_by=Delocalize.COPY,
            output_dir=output_dir)

        steps.append(s1)

        # process input files
        local_fasta, _ = s1.inputs(REFERENCE_GENOME_FASTA, f"{REFERENCE_GENOME_FASTA}.fai")

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


    combined_output_file = f"{common_prefix}.bed"
    print(f"Generating combined output file: {combined_output_file}")
    s2 = bp.new_step(
        name="combine_job",
        step_number=2,
        image=DOCKER_IMAGE,
        storage="20Gi",
        cpu=1,
        localize_by=Localize.COPY,
        delocalize_by=Delocalize.COPY,
        output_dir=output_dir,
    )

    for step in steps:
        local_bed_path, _ = s2.use_previous_step_outputs_as_inputs(step)
        s2.command(f"zcat {local_bed_path} >> {combined_output_file}")

    s2.command(f"bgzip {combined_output_file}")
    s2.command(f"tabix {combined_output_file}.gz")

    s2.output(f"{combined_output_file}.gz")
    s2.output(f"{combined_output_file}.gz.tbi")
    
    files_to_download_when_done.extend([
        (os.path.join(output_dir, f"{combined_output_file}.gz"), "results"),
        (os.path.join(output_dir, f"{combined_output_file}.gz.tbi"), "results")
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

