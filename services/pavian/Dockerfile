FROM rocker/r-ver:3.6.1

RUN apt-get update && apt-get install -y \
    sudo \
    gdebi-core \
    pandoc \
    pandoc-citeproc \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    xtail \
    wget \
    git \
    libmariadb-client-lgpl-dev \
    libbz2-dev \
    liblzma-dev \
    libxml2-dev \
    libssl-dev \
    default-libmysqlclient-dev \
    libpq-dev


# Download and install shiny server
RUN wget --no-verbose https://download3.rstudio.org/ubuntu-14.04/x86_64/VERSION -O "version.txt" && \
    VERSION=$(cat version.txt)  && \
    wget --no-verbose "https://download3.rstudio.org/ubuntu-14.04/x86_64/shiny-server-$VERSION-amd64.deb" -O ss-latest.deb && \
    gdebi -n ss-latest.deb && \
    rm -f version.txt ss-latest.deb && \
    . /etc/environment

RUN apt-get install -y libssh2-1 libssh2-1-dev

RUN R -e "install.packages('devtools');"

RUN R -e "install.packages('BiocManager');BiocManager::install('Rsamtools')"

ARG BRANCH
# Execute install_git twice, first time to download & install packages which will likely not change much, cache of this can then be used in the following step, when the up-to-date repo is pulled.
RUN R -e "devtools::install_gitlab('git@gitlab.naktuinbouw.net:bioinformatics/pavian.git', ref='${BRANCH}')"

# Use arg with current datetime to make sure everything after here is re-build without cache (so git pull is done again)
#ARG CACHE_DATE=2020-12-1
## format: YYYY-MM-DD
#RUN R -e "devtools::install_gitlab('git@gitlab.naktuinbouw.net:bioinformatics/pavian.git', ref='${BRANCH}')"

RUN R -e "sapply(list.files(system.file('shinyapp',package='pavian'),full.names=TRUE),file.copy,to='/srv/shiny-server/',recursive=TRUE)"

# Overwrite app.R from gitlab with this app.R
COPY ./app.R /srv/shiny-server/app.R

RUN chown shiny:shiny /var/lib/shiny-server
RUN mkdir -p /var/log/shiny-server
RUN chown shiny:shiny /var/log/shiny-server

ARG PAVIAN_PORT
EXPOSE $PAVIAN_PORT

CMD R -e "options(shiny.port=${PAVIAN_PORT}, shiny.host='0.0.0.0'); shiny::runApp('/srv/shiny-server/app.R')"

#COPY ./entrypoint.sh /usr/bin/shiny-server.sh
#
#CMD ["/usr/bin/shiny-server.sh"]