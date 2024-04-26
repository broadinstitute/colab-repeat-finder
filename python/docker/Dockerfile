FROM python:3.9-slim

RUN apt-get update && apt-get install -y \
	ca-certificates \
	wget \
	curl \
    tabix \
	pkg-config


COPY utils /utils
COPY repeat_finder.py /repeat_finder.py
COPY repeat_finder_pure_repeats.py /repeat_finder_pure_repeats.py

COPY requirements.txt /requirements.txt
RUN python3 -m pip install -r /requirements.txt

ENTRYPOINT ["/bin/bash"]