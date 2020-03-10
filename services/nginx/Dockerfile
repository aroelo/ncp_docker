FROM nginx:1.17.4-alpine

RUN rm /etc/nginx/conf.d/default.conf
COPY nginx.conf /etc/nginx/conf.d/nginx.conf.template
ARG POSTGRES_USER_ENV
RUN /bin/sh -c "envsubst '\$POSTGRES_USER_ENV' < /etc/nginx/conf.d/nginx.conf.template > /etc/nginx/conf.d/nginx.conf"