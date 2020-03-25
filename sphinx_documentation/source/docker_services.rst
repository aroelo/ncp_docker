.. _docker_services:

Docker services
===================

There are several services build/run by docker-compose, these are all specified in the docker-compose file.
Each service has it's own Dockerfile and corresponds to a single container.

The different services for our pavian visualisation are explained below.

jbrowse
----------
| The jbrowse service hosts a jbrowse server that can be used to visualise alignments.
| The data location is hardcoded to /jbrowse/data inside the container.
| In the docker-compose file you can see that this volume links to the PAVIAN_OUT_H env variable from the host.

nginx
-------


pavian
-------

web
----

