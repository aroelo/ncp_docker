FROM nginx:1.17.4-alpine

ARG DOCKER_GID
ARG DOCKER_UID

RUN groupadd -r -g $DOCKER_UID app && useradd -r -u $DOCKER_GID -g app app

RUN rm /etc/nginx/conf.d/default.conf
COPY nginx.conf /etc/nginx/conf.d/nginx.conf.template
ARG POSTGRES_USER_ENV
ARG FLASK_PORT_ENV
#substitute postgres user & port
RUN /bin/sh -c "envsubst '\$POSTGRES_USER_ENV \$FLASK_PORT_ENV' < /etc/nginx/conf.d/nginx.conf.template > /etc/nginx/conf.d/nginx.conf"