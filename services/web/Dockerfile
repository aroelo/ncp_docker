FROM ubuntu:18.04

ARG BLAST_VERSION

# create directory for the app user
RUN mkdir -p /home/app

# create the app user & group
RUN groupadd -r app && useradd -r -g app app

# create the appropriate directories
ENV HOME=/home/app
ENV APP_HOME=/home/app/web
RUN mkdir $APP_HOME
WORKDIR $APP_HOME

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

# netcat needed for entrypoint file
RUN apt-get update \
    && apt-get install -y netcat

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

# install bedtools to create bigwig
RUN apt-get update \
    && apt-get install -y bedtools \
    && wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig \
    && chmod +x bedGraphToBigWig

# install requirements
COPY ./requirements.txt .
RUN pip install -r requirements.txt

# use blastcmd
RUN wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/${BLAST_VERSION}/ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz  && \
    tar xzf ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz && \
    mv ncbi-blast-${BLAST_VERSION}+ /blast

# Set path for samtools, blast & bedGraphToBigWig
ENV PATH="/home/app/web/samtools-1.10:/blast/bin:/home/app/web:${PATH}"

# Create directories to mount blast & pavian data
RUN mkdir -p /blast_db/nt /blast_db/refseq /pavian/out /pavian/in

# copy project
#COPY . $APP_HOME

# chown all the files to the app user
RUN chown -R app:app $APP_HOME

# change to the app user
USER app

# run entrypoint.sh
ENTRYPOINT ["/home/app/web/entrypoint.sh"]
