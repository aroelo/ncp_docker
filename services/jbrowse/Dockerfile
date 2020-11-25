FROM nginx
ARG DOCKER_UID
ARG DOCKER_GID
ARG JBROWSE_PORT
RUN echo DOCKER_UID:$DOCKER_UID
RUN echo JBROWSE_PORT:$JBROWSE_PORT
RUN echo DOCKER_GID:$DOCKER_GID
RUN groupadd -r $GID && useradd -r -g $GID $UID
ENV DEBIAN_FRONTEND noninteractive

RUN mkdir -p /usr/share/man/man1 /usr/share/man/man7

RUN apt-get -qq update --fix-missing
RUN apt-get --no-install-recommends -y install build-essential zlib1g-dev libxml2-dev libexpat-dev postgresql-client libpq-dev libpng-dev wget unzip perl-doc ca-certificates

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-4.6.14-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /conda/ && \
    rm ~/miniconda.sh

ENV PATH="/conda/bin:${PATH}"

RUN conda install -y --override-channels --channel iuc --channel conda-forge --channel bioconda --channel defaults jbrowse=1.16.6

RUN rm -rf /usr/share/nginx/html && ln -s /conda/opt/jbrowse/ /usr/share/nginx/html && \
    ln -s /jbrowse/data /conda/opt/jbrowse/data && \
    sed -i '/include += {dataRoot}\/tracks.conf/a include += {dataRoot}\/datasets.conf' /conda/opt/jbrowse/jbrowse.conf && \
    sed -i '/include += {dataRoot}\/tracks.conf/a include += {dataRoot}\/../datasets.conf' /conda/opt/jbrowse/jbrowse.conf

ARG JBROWSE_PORT
# Replace default port 80 with port 5004
RUN sed -i 's/listen       80;/listen       '$JBROWSE_PORT';/g' /etc/nginx/conf.d/default.conf

CMD ["nginx", "-g", "daemon off;"]