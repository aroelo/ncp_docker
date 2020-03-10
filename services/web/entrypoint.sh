#!/bin/bash

echo "Creating the database tables..."
#    python manage.py drop_db
python manage.py create_db
python manage.py seed_pavian_db
echo "Tables created"

exec "$@"