FROM ubuntu:16.04
RUN apt-get update && apt-get install -y \
  build-essential \
  cmake \
  libcairo2-dev \
  libeigen2-dev \
  wget \
  zlib1g-dev \
  git
RUN git clone https://github.com/n-yoshikawa/openbabel.git && \
  cd openbabel/ && \
  git checkout fragment-based-builder-cod && \
  mkdir build && cd build/ && \
  cmake .. && make && make install
RUN cd /home/ && git clone https://github.com/n-yoshikawa/openbabel-paper.git
