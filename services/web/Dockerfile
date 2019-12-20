#FROM ubuntu:18.04 as blastbuild
#
#ARG version
#ARG num_procs=8
#LABEL Description="NCBI BLAST" Vendor="NCBI/NLM/NIH" Version=${version} Maintainer=camacho@ncbi.nlm.nih.gov
#
#USER root
#WORKDIR /root/
#
#RUN apt-get -y -m update && apt-get install -y build-essential curl libidn11 libnet-perl perl-doc liblmdb-dev
#
#RUN curl -s ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/${version}/ncbi-blast-${version}+-src.tar.gz | \
# tar xzf - && \
# cd ncbi-blast-${version}+-src/c++ && \
# ./configure --with-mt --with-strip --with-optimization --with-dll --with-experimental=Int8GI --with-flat-makefile --with-openmp --without-vdb --without-gnutls --without-gcrypt --prefix=/blast && \
# cd ReleaseMT/build && \
# make -j ${num_procs} -f Makefile.flat blastdb_aliastool.exe blastdbcheck.exe blastdbcmd.exe blast_formatter.exe blastn.exe blastp.exe blastx.exe convert2blastmask.exe deltablast.exe dustmasker.exe makeblastdb.exe makembindex.exe makeprofiledb.exe psiblast.exe rpsblast.exe rpstblastn.exe segmasker.exe tblastn.exe tblastx.exe windowmasker.exe
#

# pull official base image
FROM ubuntu:18.04

ARG BLAST_VERSION

# set work directory
WORKDIR /usr/src/app

# set environment variables
ENV PYTHONDONTWRITEBYTECODE 1
ENV PYTHONUNBUFFERED 1
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8


RUN apt-get update \
    && apt-get install -y python3-pip python3-dev \
    && cd /usr/local/bin \
    && ln -s /usr/bin/python3 python \
    && pip3 install --upgrade pip
## add bash
#RUN apk update \
#    && apk add --no-cache bash
#
## install psycopg2 dependencies
#RUN apk update \
#    && apk add postgresql-dev gcc python3-dev musl-dev
#
## install pysam dependencies
#RUN apk update \
#    && apk add zlib-dev bzip2-dev xz-dev curl-dev make

# install pybigwig dependencies
RUN apt-get update \
    && apt-get install -y libcurl4-openssl-dev zlib1g-dev

# install blastdbcmd dependencies
RUN apt-get update \
    && apt-get install -y liblmdb-dev

# install samtools dependencies
RUN apt-get update \
    && apt-get install -y wget libncurses5 libncurses5-dev libbz2-dev liblzma-dev && \
    wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2 && \
    tar xjf samtools-1.10.tar.bz2 && \
    cd samtools-1.10 && \
    ./configure && \
    make && \
    make install

# install requirements
COPY ./requirements.txt /usr/src/app/requirements.txt
RUN pip install -r requirements.txt

# use blastcmd

RUN wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/${BLAST_VERSION}/ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz  && \
    tar xzf ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz && \
    mv ncbi-blast-${BLAST_VERSION}+ /blast
#RUN mkdir -p /blast/bin /blast/lib
#COPY --from=blastbuild /root/ncbi-blast-${version}+-src/c++/ReleaseMT/bin /blast/bin
#COPY --from=blastbuild /root/ncbi-blast-${version}+-src/c++/ReleaseMT/lib /blast/lib

# Set path for samtools and blast
ENV PATH="samtools-1.10.1/bin:/blast/bin:${PATH}"

# Create directories to mount blast & pavian data
RUN mkdir -p /blast_db/nt /blast_db/refseq /pavian/out /pavian/in

# copy project
COPY . /usr/src/app/

# run entrypoint.sh
ENTRYPOINT ["/usr/src/app/entrypoint.sh"]