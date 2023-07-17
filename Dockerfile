FROM python:3.9-slim

RUN apt-get update && \
    apt-get install -y build-essential wget git && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /angler4

#install blast binaries
RUN wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.13.0/ncbi-blast-2.13.0+-x64-linux.tar.gz && \
    tar -xvf ncbi-blast-2.13.0+-x64-linux.tar.gz && \
    rm ncbi-blast-2.13.0+-x64-linux.tar.gz && \
    cp ncbi-blast-2.13.0+/bin/blastn /usr/bin/ && \
    rm -rf ncbi-blast-2.13.0+

#install muscle v3
RUN wget https://drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86linux64.tar.gz && \
    tar -xvf muscle3.8.31_i86linux64.tar.gz && \
    rm muscle3.8.31_i86linux64.tar.gz && \
    chmod +x muscle3.8.31_i86linux64 && \
    mv muscle3.8.31_i86linux64 /usr/local/bin/muscle

# install poetry
ENV POETRY_HOME /opt/poetry
RUN python3 -m venv $POETRY_HOME && \
    $POETRY_HOME/bin/pip install poetry==1.2.0 && \
    mkdir /poetry-env

# local user needs to create the venv
ENV VIRTUAL_ENV /poetry-env
RUN python3 -m venv $VIRTUAL_ENV

ENV PATH="$VIRTUAL_ENV/bin:$POETRY_HOME/bin:$PATH"

# pull in NUPACK source and install python dependencies
RUN git clone https://github.com/choulab/NUPACK.git

COPY ./app /angler4

RUN poetry install
