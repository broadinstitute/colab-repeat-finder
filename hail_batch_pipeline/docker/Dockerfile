FROM python:3.9-slim

RUN apt-get update && apt-get install -y \
	ca-certificates \
	wget \
	curl \
    tabix \
	pkg-config


COPY utils /utils
COPY perfect_repeat_finder.py /perfect_repeat_finder.py

COPY requirements.txt /requirements.txt
RUN python3 -m pip install -r /requirements.txt

ENTRYPOINT ["/bin/bash"]