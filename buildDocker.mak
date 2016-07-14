
DOCKER_HUB_ID=stuartlab

# build a docker image
paradigm-tools-python:
	docker build -f Dockerfile_python_tools --tag $(DOCKER_HUB_ID)/paradigm-tools-python .
