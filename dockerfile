# Using Linux ubuntu 20.04
FROM ubuntu:20.04

WORKDIR /usr/src/app

RUN apt update && apt install -y \
	gcc \
	g++ \
	make \
	cmake \
	libjsoncpp-dev

