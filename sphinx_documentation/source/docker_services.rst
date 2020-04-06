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
Nginx is only used in the production environment. It acts as a web server and reverse proxy and interacts with gunicorn.

pavian
-------
| Pavian is the actual app that is used for visualisation. This app is created in R with shiny.
| It has one custom app.R file that is used to start the app and define custom settings such as database addresses/logins.

web
----
| The web service contains the flask app. The flask app processes the additional features that we have implemented to visualise on the taxid level.
| When clicking a link in one of the yellow pop-ups in the 'comparison' overview in pavian, the request will be sent through to the flask app.


