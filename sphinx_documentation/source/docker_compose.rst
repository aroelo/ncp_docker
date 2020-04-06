.. _docker_compose:

Docker-compose files
======================

The docker-compose file was created at version 3.7 and contains build & run information of several different services (see :ref:`docker_services`).

Many of the parameters are obtained from ENV values that are defined in the '.env' file.

The key differences between the dev and prod as of this moment are:
    - The BRANCH arg in the 'pavian' service that determines which pavian branch is cloned from github.
    - The extra 'nginx' service in the prod environment.
    - Consequently, some different ports that are used with the 'nginx' service mentioned above.
    - Which postgres database is used in the 'web' service.

