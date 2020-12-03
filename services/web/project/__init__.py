from flask import (
    Flask,
    send_from_directory,
    render_template,
)
from flask_sqlalchemy import SQLAlchemy
import os


app = Flask(__name__)
app.config.from_object("project.config.Config")
db = SQLAlchemy(app)


@app.errorhandler(500)
def taxid_files_not_created(e):
    return render_template('500.html',
                           taxid=e.description.split('item:')[1],
                           sample=e.description.split('item:')[2],
                           traceback=e.description.split('item:')[3].replace('\n', '<br>')), 500


@app.route('/favicon.ico')
def favicon():
    return send_from_directory(os.path.join(app.root_path, 'static'),
                               'favicon.ico', mimetype='image/vnd.microsoft.icon')


# Import modules using blueprint handler variable
from project.download_pavian_data.controllers import download_pavian_data_blueprint as download_pavian_data_module
from project.viewreads.controllers import get_json_data_blueprint as get_json_module
from project.blast_reads.controllers import blast_reads_blueprint as blast_reads_module

# Register modules
app.register_blueprint(download_pavian_data_module)
app.register_blueprint(get_json_module)
app.register_blueprint(blast_reads_module)
