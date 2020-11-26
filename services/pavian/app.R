if (!require(Rsamtools)) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("Rsamtools")
  library(Rsamtools)
}

if (!dir.exists(rappdirs::user_config_dir("pavian", expand = FALSE))) {
  dir.create(rappdirs::user_config_dir("pavian", expand = FALSE),
             recursive = TRUE)
}

### Option specifications
## Shiny options
# Specify the maximum web request size, which serves as a size limit for file uploads. If unset, the maximum request size defaults to 5MB
# see https://shiny.rstudio.com/reference/shiny/latest/shiny-options.html for global shiny options
options(shiny.maxRequestSize = 50 * 1024 ^ 2) # set to 50 MB
## do not set shiny.maxRequestSize here, because it overrides user options!

## DT options
# see https://datatables.net/reference/option/
options(
  DT.options = list(
    pageLength = 15,
    stateSave = TRUE,
    searchHighlight = TRUE,
    #scrollX = TRUE,
    dom = 'Bfrtip',
    ## Define the table control elements to appear
    #  B - Buttons
    #  f - filtering input
    #  r - processing display element
    #  t - The table!
    #  i - Table information summary
    #  p - pagination control
    lengthMenu = list(c(15, 25, 50, 100), c('15', '25', '50', '100')),
    search = list(regex = TRUE, caseInsensitive = TRUE)
  )
)

# Uses environment variables from .env.dev & env.dev.db or .env.prod & .env.prod.db, determined by env_file(s) in pavian service
pavian::runApp(server_dir = Sys.getenv("PAVIAN_IN"),
               flask_host = Sys.getenv("HOST_IP"), flask_port = Sys.getenv("NGINX_PORT"),
               db_type = "Postgresql", db_name = Sys.getenv("POSTGRES_DB"), db_host = Sys.getenv('SQL_HOST'),
               db_user = Sys.getenv('POSTGRES_USER'), db_passwd = Sys.getenv('POSTGRES_PASSWORD'))
# HOST_IP should be:
# paste("nginx",  Sys.getenv("HOST_DOMAIN") )
# but we don't have an "nginx-dev.naktuinbouw.cloud" and I'm not sure we need one?
