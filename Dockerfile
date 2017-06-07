FROM fedora:latest
MAINTAINER Steffen Vogel <stvogel@eonerc.rwth-aachen.de>

RUN dnf -y install \
	gcc-c++ \
	make \
	eigen3-devel \
	doxygen \
	python-pandas \
    python3-numpy

ADD https://villas.fein-aachen.org/packages/villas.repo /etc/yum.repos.d/

RUN dnf -y install villas-node-devel
