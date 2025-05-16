FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    libtool \
    autoconf \
    automake \
    software-properties-common \
    && rm -rf /var/lib/apt/lists/*

# Install GCC
RUN add-apt-repository universe && apt update && apt install -y gcc-8 g++-8 \
    && update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-8 100 \
    && update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-8 100 \
    && update-alternatives --config gcc \
    && update-alternatives --config g++

# Install OpenMPI
WORKDIR /tmp
RUN wget https://download.open-mpi.org/release/open-mpi/v4.0/openmpi-4.0.1.tar.gz && \
    tar -xzf openmpi-4.0.1.tar.gz && \
    cd openmpi-4.0.1 && \
    ./configure && make all && make install && ldconfig

# Install FFTW
RUN apt-get install -y libfftw3-dev

RUN apt-get purge -y \
    wget \
    libtool \
    autoconf \
    automake \
    software-properties-common \
    && apt-get clean

WORKDIR /project
