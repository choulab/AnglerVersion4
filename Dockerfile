FROM python:3.9-slim

RUN apt-get update && \
    apt-get install -y build-essential wget git && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /tmp

#install blast binaries
RUN wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.13.0+-x64-linux.tar.gz && \
    tar -xvf ncbi-blast-2.13.0+-x64-linux.tar.gz && \
    cp ncbi-blast-2.13.0+/bin/blastn /usr/bin/ && \
    rm -rf ncbi-blast-2.13.0+

#install muscle v3
RUN wget https://drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86linux64.tar.gz && \
    tar -xvf muscle3.8.31_i86linux64.tar.gz && \
    chmod +x muscle3.8.31_i86linux64 && \
    mv muscle3.8.31_i86linux64 /usr/local/bin/muscle

# install nupack from github 
RUN git clone https://github.com/beliveau-lab/NUPACK.git && \
    python3 -m pip install NUPACK/src/package/nupack-4.0.0.23-cp39-cp39-manylinux2014_x86_64.whl && \
    rm -rf NUPACK

COPY requirements.txt ./requirements.txt

# install remaining dependencies
RUN pip install -r requirements.txt

WORKDIR /code

COPY . .

