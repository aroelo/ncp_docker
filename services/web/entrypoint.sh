#!/bin/bash

echo "Waiting for postgres..."

while ! nc -z $SQL_HOST $SQL_PORT; do
  sleep 0.1
done

echo "PostgreSQL started"
python manage.py create_db
echo "DB created"
python manage.py seed_pavian_db
echo "DB seeded"
exec "$@"