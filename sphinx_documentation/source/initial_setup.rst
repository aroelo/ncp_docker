Initial set-up
=================

The results from the classification pipeline can be visualised with our custom pavian version, which has some extra functionalities compared to the original pavian.

| The custom pavian can be easily set-up with docker-compose. The pavian docker repository is available at:  https://github.com/aroelo/ncp_docker

Branches: master vs blast_v5
---------------------------------

| Currently, there are two branches available. The blast_v5 branch works with the latest blast database and should be used for any new set-ups.
| The difference lies in the naming of the blast database that is used with the blastdbcmd.

| For the latest blast databases there is no '_v5' suffix when the database is at version 5.
| This isn't the case for the blast databases that are currently used at Naktuinbouw, so here the 'master' branch should be used.
| The suffix has to be hardcoded in both of the docker-compose files in the BLASTDB_NT environment variable, because of how the blastdbcmd works.

| If the directory path where e.g. the nt blast database is stored is `/6_db_ssd/blastdb/nt`, then the full path that the blastdbcmd uses is `/6_db_ssd/blastdb/nt/nt_v5`.
| For the latest blast database release it would just be `/6_db_ssd/blastdb/nt/nt`.

.. note::
    The github repository is mirrored from the internal gitlab repository at http://gitlab.naktuinbouw.net/bioinformatics/pavian_docker/

    All sensitive information has been removed from git history with https://rtyley.github.io/bfg-repo-cleaner/, see http://gitlab.naktuinbouw.net/bioinformatics/pavian_docker/issues/7

Requirements
-------------

- docker-compose version 1.25.1 or later
- open/accessible ports for (defaults):
    - pavian (3838)
    - jbrowse (5004)
    - flask app (5000)
    - gunicorn (5000)
    - nginx (1337)
- output directory on host should be writable for all users

External files and databases
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
There are some external files necessary, these can be found in this repository in the 'data' directory':

- taxa.sqlite
- taxa.sqlite.traverse.pkl

| Two databases are needed, these are likely already created for the classification pipeline.
| If this is not the case, they can be created with the 'download_databases.sh' script that is available in the classification repository.
| The databases are:

- nt database
- wgs database

| See :ref:`reference_databases_setup` for instructions on running this script if you haven't done this yet.
| For the pavian visualisation it is `not` necessary to include the '-wgs_fasta' option.

.. note::
    We should include the external files & script in the docker git repository or as a submodule.
    Alternatively we could create a new git repository that only contains shared files/script with the classification pipeline & pavian docker.

.. note::
    Perhaps it would be nice to add some example files to this repository, so they can be used to test functionality in pavian

Creating the .env file
------------------------

Inside this repository there is a 'env' file which should be modified with the correct values for paths/ports etc.

.. note::
    Even though some ports seem to be custom, they are still hardcoded in some parts of the code.
    This is the case for:

        1. jbrowse (5004)
        2. pavian (3838)
        3. nginx (1337)
        4. flask (5000)

        1. Should be modified in `services/jbrowse/Dockerfile` with an arg based on env value
        2. Should be modified in `services/pavian/Dockerfile` with an arg based on env value
        3. Should be modified in `docker-compose.prod.yml` with env value
        4. Should be modifief in `docker-compose.yml` and `docker-compose.prod.yml` with env value

After changing the values in the env file, make sure to rename it to '.env', so the values are automatically read-in during docker-compose.

Docker-compose commands
----------------------------

Building the docker containers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Make sure that you have entered the directory where the docker-compose files are present. With the docker-compose file several services can be build at once.
For a full explanation of the docker-compose files (dev & prod) and the services see :ref:`docker_compose` and :ref:`docker_services`.

To build the dev environment:

`docker-compose build --build-arg CACHE_DATE=$(date +%s)`

To build the prod environment:

`docker-compose -f docker-compose.prod.yml build --build-arg CACHE_DATE=$(date +%s)`

The CACHE_DATE building arg is present to ensure that the external pavian repository is cloned from github when re-building the docker container.

On the initial build this will be done anyway, however when re-building at a later stage, docker-compose will use cache. If the build-arg is not provided it will use the current cached version of the pavian container, even if there is a newer version available from github.

A re-build of the docker-container(s) is preferred whenever there is either an update in this repository or in the external pavian repository.

Running the docker containers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Dev environment:

`docker-compose up`

Prod environment:

`docker-compose -f docker-compose.prod.yml up`


Add the '-d' parameter to run the containers 'detached' in the background.

If there are any issues with the databases, it is recommend to shut down any running containers (see below) and retry.

Shutting down the docker containers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Dev environment:

`docker-compose down -v`

Prod environment:

`docker-compose -f docker-compose.prod.yml down -v`

