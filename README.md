Following this link to set up dev and production environment:  
https://testdriven.io/blog/dockerizing-flask-with-postgres-gunicorn-and-nginx/

**Using this set-up somewhere else**   
 *  Use docker-compose version 1.25.1 or later
 *  Git clone repo from  https://github.com/aroelo/ncp_docker.git  
 *  Use the blast_v5 branch
 *  Then change:  
    *  The settings in the env to their correct settings
    *  Rename 'env' file to '.env'
 *  Make sure all ports are open/accessible:
    *  $PAVIAN_PORT_ENV
    *  $NGINX_PORT_ENV
    *  $JBROWSE_PORT_ENV
    *  $FLASK_PORT_ENV


 If there are any issues with the databases, it is recommend to do docker-compose down -v and then docker-compose up.  

Useful cmds in Naktest-01, `/5_workspace/docker/pavian_docker` using conda bio environment (just for docker-compose to work):  
**For dev part - default port 3839**  
Build  
`docker-compose -p dev build --build-arg CACHE_DATE=$(date +%s)`  
Up  
`docker-compose -p dev up -d`  
Down  
`docker-compose -p dev down -v`  
Check logs  
`docker-compose -p dev logs -f`


**For prod part - default port 3838**  
 Build  
 `docker-compose build --build-arg CACHE_DATE=$(date +%s)`  
 Up  
 `docker-compose up -d`  
 Down  
 `docker-compose down -v`  
 Check logs  
`docker-compose logs -f`  
 
 
**Debugging**  
To debug the script in e.g. pycharm that creates all output files per taxid, run `services/web/project/run.py`  with the -human argument and all other required arguments.  
Choose '-human main' to just run the main command, creating output files, or choose '-human app' to start up a flask app locally outside the docker environment.


**Syncing pavian input sql database**  
Sync is done daily with a cronjob from aroe.  
(-T fixes the input device is not a TTY .. see https://github.com/docker/compose/issues/5696)  
`30 7 * * * cd /5_workspace/docker/pavian_docker/ && docker-compose exec -T web python manage.py seed_pavian_db > /home/aroe/logs/seed_pavian_db.log 2>&1`