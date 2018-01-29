FROM ubuntu:latest
RUN apt-get update
RUN apt-get install -y build-essential zlib1g-dev

RUN cd /bin
COPY generic-metal-2011-03-25.tar.gz ./
RUN tar -xzf ./generic-metal-2011-03-25.tar.gz && rm generic-metal-2011-03-25.tar.gz && cd generic-metal && make all && mv executables/metal /bin && cd /