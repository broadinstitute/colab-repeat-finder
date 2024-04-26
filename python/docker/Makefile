TAG1 = weisburd/colab-repeat-finder:v1
TAG2 = weisburd/colab-repeat-finder:latest

all: build push

build:
	cd .. && docker build -t $(TAG1) -t $(TAG2) -f docker/Dockerfile .

push:
	docker push $(TAG1)
	docker push $(TAG2)
