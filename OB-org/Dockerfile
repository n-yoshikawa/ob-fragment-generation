FROM ubuntu:16.04
RUN apt-get update && apt-get install -y \
  build-essential \
  cmake \
  libcairo2-dev \
  libeigen2-dev \
  wget \
  zlib1g-dev \
  git
RUN wget https://github.com/openbabel/openbabel/archive/openbabel-2-4-1.tar.gz && \
  tar xf openbabel-2-4-1.tar.gz && \
  cd openbabel-openbabel-2-4-1/ && \
  mkdir build && cd build/ && \
  cmake .. && make && make install
RUN cd /home/ && git clone https://github.com/n-yoshikawa/openbabel-paper.git
