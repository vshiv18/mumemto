# Use a base image with the desired dependencies
FROM ubuntu:latest

# Install system dependencies
RUN apt-get update && apt-get install -y \
    g++ \
    gcc \
    cmake \
    git \
    wget \
    zlib1g-dev \
    libboost-all-dev

# Install Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda && \
    rm Miniconda3-latest-Linux-x86_64.sh

# Add Miniconda to the PATH
ENV PATH="/opt/conda/bin:${PATH}"

# Create a new conda environment and activate it
RUN conda create -n myenv python=3.8 && \
    echo "conda activate myenv" >> ~/.bashrc
SHELL ["/bin/bash", "--login", "-c"]

# Install Python packages
RUN conda install -n myenv -y \
    tqdm \
    matplotlib

# Clone the repository
RUN git clone https://github.com/vshiv18/mumemto && \
    cd mumemto && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make install

# Add the build directory to the PATH
ENV PATH="/mumemto/build:${PATH}"