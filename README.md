Following this link:  
https://testdriven.io/blog/dockerizing-flask-with-postgres-gunicorn-and-nginx/

Docker daemon is hosted on port 2376

Useful cmds in Naktest-01, ~/PycharmProjects/flask-on-docker using classification-pipeline-1.0.3:  
**For dev part**  
Build  
`docker-compose build --build-arg SSH_PRIVATE_KEY="$(cat /home/aroe/.ssh/id_rsa)" --build-arg SSH_PUBLIC_KEY="$(cat /home/aroe/.ssh/id_rsa.pub)" --build-arg CACHE_DATE=$(date +%s)`  
Up  
`docker-compose up -d`  
Seed pavian datatable  
`docker-compose exec web python manage.py seed_pavian_db`  
Down  
`docker-compose down -v`  
Check logs  
`docker-compose logs -f`

