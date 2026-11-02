# Using Linux ubuntu 20.04
FROM ubuntu:22.04

WORKDIR /usr/src/app

RUN apt update && apt install -y \
	gcc \
	g++ \
	make \
	cmake \
	libjsoncpp-dev

RUN cmake -S . -B build && cmake --build build
