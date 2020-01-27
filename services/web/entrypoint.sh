#!/bin/bash

if [ "$DATABASE" = "postgres" ]
then
    echo "Waiting for postgres..."

    while ! nc -z $SQL_HOST $SQL_PORT; do
      sleep 0.1
    done

    echo "PostgreSQL started"
    python manage.py create_db
    python manage.py seed_pavian_db
fi

if [ "$FLASK_ENV" = "development" ]
then
    echo "Creating the database tables..."
#    python manage.py drop_db
    python manage.py create_db
    python manage.py seed_pavian_db
    echo "Tables created"
fi

exec "$@"