**Pavian is hosted on port 3838**

Following this link to set up dev and production environment:  
https://testdriven.io/blog/dockerizing-flask-with-postgres-gunicorn-and-nginx/

**Using this set-up somewhere else**   
 *  Use docker-compose version 1.25.1 or later
 *  Git clone repo from  https://github.com/aroelo/ncp_docker.git  
 *  Use the blast_v5 branch
 *  Make sure all ports are open/accessible:
    *  3838 - pavian (hardcoded in Dockerfile)
    *  1337 - nginx (hardcoded)
    *  $JBROWSE_PORT_ENV
    *  $FLASK_PORT_ENV

 
 Then change:  
 *  The settings in the env to their correct settings
 *  Rename 'env' file to '.env'
 
 If there are any issues with the databases, it is recommend to do docker-compose down -v and then docker-compose up.  

Useful cmds in Naktest-01, `/5_workspace/docker/pavian_docker` using conda bio environment (just for docker-compose to work):  
**For dev part**  
Build  
`docker-compose build --build-arg CACHE_DATE=$(date +%s)`  
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
 `docker-compose -f docker-compose.prod.yml build --build-arg CACHE_DATE=$(date +%s)`  
 Up  
 `docker-compose -f docker-compose.prod.yml up -d`  
 Down  
 `docker-compose -f docker-compose.prod.yml down -v`  
 Check logs  
`docker-compose -f docker-compose.prod.yml logs -f`  
 
 
**Debugging**  
To debug the script in e.g. pycharm that creates all output files per taxid, run `services/web/project/run.py`  with the -human argument and all other required arguments.  
Choose '-human main' to just run the main command, creating output files, or choose '-human app' to start up a flask app locally outside the docker environment.
For the imports to work you have to mark the `services/web` directory as sources root in pycharm.


**Syncing pavian input sql database**  
Sync is done daily with a cronjob from aroe.  
(-T fixes the input device is not a TTY .. see https://github.com/docker/compose/issues/5696)  
`30 7 * * * cd /5_workspace/docker/pavian_docker/ && docker-compose exec -T web python manage.py seed_pavian_db > /home/aroe/logs/seed_pavian_db.log 2>&1`