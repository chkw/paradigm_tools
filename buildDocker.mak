
DOCKER_HUB_ID=stuartlab

# build a docker image
paradigm-tools:
	docker build -f Dockerfile_paradigm_tools --tag $(DOCKER_HUB_ID)/paradigm-tools .
