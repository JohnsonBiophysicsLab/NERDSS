# Use an Intel-based Ubuntu image
FROM --platform=linux/amd64 ubuntu:20.04

# Avoid prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Update and install necessary tools
RUN apt-get update && apt-get upgrade -y && apt-get install -y \
    build-essential \
    cmake \
    git \
    gdb \
    google-perftools \
    libgoogle-perftools-dev \
    libgsl-dev \
    libopenmpi-dev \
    openmpi-bin \
    libomp-dev \
    && rm -rf /var/lib/apt/lists/*

# Compile and install Google Test
WORKDIR /usr/src/gtest
RUN git clone https://github.com/google/googletest.git . \
    && mkdir build \
    && cd build \
    && cmake .. \
    && make \
    && make install \
    && cd .. \
    && rm -rf build .git

# Create a new user
RUN useradd -ms /bin/bash myuser

# Create and set ownership of the working directory
RUN mkdir /nerdss_mpi && chown myuser:myuser /nerdss_mpi

# Switch to the new user
USER myuser

# Set up your working directory
WORKDIR /nerdss_mpi
