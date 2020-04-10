FROM nginx:1.17.4-alpine

RUN rm /etc/nginx/conf.d/default.conf
COPY nginx.conf /etc/nginx/conf.d/nginx.conf.template
ARG POSTGRES_USER_ENV
ARG FLASK_PORT_ENV
#substitute postgres user & port
RUN /bin/sh -c "envsubst '\$POSTGRES_USER_ENV \$FLASK_PORT_ENV' < /etc/nginx/conf.d/nginx.conf.template > /etc/nginx/conf.d/nginx.conf"