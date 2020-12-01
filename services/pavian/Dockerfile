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


# add credentials on build
ARG SSH_PRIVATE_KEY
ARG SSH_PUBLIC_KEY
RUN mkdir /root/.ssh/
RUN echo "${SSH_PRIVATE_KEY}" > /root/.ssh/id_rsa
RUN echo "${SSH_PUBLIC_KEY}" > /root/.ssh/id_rsa.pub

#Gitlab requires a private key with strict permission settings
RUN chmod 600 /root/.ssh/id_rsa
RUN chmod 600 /root/.ssh/id_rsa.pub

# make sure your domain is accepted
RUN touch /root/.ssh/known_hosts
RUN ssh-keyscan gitlab.naktuinbouw.net >> /root/.ssh/known_hosts

ARG BRANCH
# Execute install_git twice, first time to download & install packages which will likely not change much, cache of this can then be used in the following step, when the up-to-date repo is pulled.
RUN R -e "creds <- git2r::cred_ssh_key('/root/.ssh/id_rsa.pub', '/root/.ssh/id_rsa');devtools::install_git('git@gitlab.naktuinbouw.net:bioinformatics/pavian.git', ref='${BRANCH}', credentials=creds)"

# Use arg with current datetime to make sure everything after here is re-build without cache (so git pull is done again)
ARG CACHE_DATE=2020-12-1
## format: YYYY-MM-DD
RUN R -e "creds <- git2r::cred_ssh_key('/root/.ssh/id_rsa.pub', '/root/.ssh/id_rsa');devtools::install_git('git@gitlab.naktuinbouw.net:bioinformatics/pavian.git', ref='${BRANCH}', credentials=creds)"


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