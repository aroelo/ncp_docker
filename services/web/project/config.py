import os


basedir = os.path.abspath(os.path.dirname(__file__))


class Config(object):
    SQLALCHEMY_DATABASE_URI = os.getenv("DATABASE_URL", "sqlite://")
    SQLALCHEMY_TRACK_MODIFICATIONS = False
    # STATIC_FOLDER = f"{os.getenv('APP_FOLDER')}/project/static"
    # MEDIA_FOLDER = f"{os.getenv('APP_FOLDER')}/project/media"

    HOST_IP = f"{os.getenv('HOST_IP')}"
    JBROWSE_PORT = f"{os.getenv('JBROWSE_PORT')}"
    PAVIAN_IN = f"{os.getenv('PAVIAN_IN')}"
    PAVIAN_OUT = f"{os.getenv('PAVIAN_OUT')}"
    TAXA_SQLITE = f"{os.getenv('TAXA_SQLITE')}"
    BLASTDB_NT = f"{os.getenv('BLASTDB_NT')}"
    BLASTDB_REFSEQ = f"{os.getenv('BLASTDB_REFSEQ')}"