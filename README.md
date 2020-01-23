**Pavian is hosted on ***REMOVED***:3838**

Following this link to set up dev and production environment:  
https://testdriven.io/blog/dockerizing-flask-with-postgres-gunicorn-and-nginx/

Docker daemon is hosted on port 2376

Useful cmds in Naktest-01, `/5_workspace/docker/pavian_docker` using conda bio environment (just for docker-compose to work):  
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


**For prod part**  
 Build  
 `docker-compose -f docker-compose.prod.yml build --build-arg SSH_PRIVATE_KEY="$(cat /home/aroe/.ssh/id_rsa)" --build-arg SSH_PUBLIC_KEY="$(cat /home/aroe/.ssh/id_rsa.pub)" --build-arg CACHE_DATE=$(date +%s)`  
 Up  
 `docker-compose -f docker-compose.prod.yml up -d`  
 Down  
 `docker-compose -f docker-compose.prod.yml down -v`  
 Check logs  
`docker-compose -f docker-compose.prod.yml logs -f`  
 
 
 **Using this set-up somewhere else**   
 First, git clone repo from  http://gitlab.naktuinbouw.net/bioinformatics/pavian_docker.git  
 
 Then change:  
 *  The directory/file locations in the .env to their correct corresponding locations
 *  The host IP in .env.dev and/or .env.prod
 
 
 If there are any issues with the databases, it is recommend to do docker-compose down -v and then docker-compose up.  
 
**Debugging**  
To debug the script in e.g. pycharm that creates all output files per taxid, run `services/web/project/__init__.py`  with the --human argument and all other required arguments.  
