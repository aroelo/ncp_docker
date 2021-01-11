from project.blast_reads import NCBIWWWcustom
from flask import (
    Blueprint,
    request,
    redirect,
    abort
)
import subprocess
import os


def create_query(readids, subdir, taxid, VA):
    if VA:
        fasta_file_path = os.path.join(subdir, taxid + ".ref.fa")
    else:
        fasta_file_path = os.path.join(subdir, taxid + ".fasta")

    cmd1 = "mawk '/%s/ {print $0; getline; print $0}' %s" % (readids.strip("|"), fasta_file_path)
    ps1 = subprocess.Popen(cmd1, shell=True, executable='/bin/bash',
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    ps1.wait()
    stdout, stderr = ps1.communicate()
    stdout = stdout.decode().replace('\r', '\n')
    if not stdout:
        e = f'item:{str(taxid)}item:{os.path.basename(subdir)}item:No reads found\n'
        abort(500, e)
    return stdout


def blast_query(fasta_string):
    get_url = NCBIWWWcustom.qblast("blastn", "nt", fasta_string, format_type="HTML")
    return get_url


def get_parameters_blast():
    readids = request.args.get('readids')
    readids = str(readids)

    subdir = request.args.get('subdir')
    subdir = str(subdir)

    taxid = request.args.get('taxid')
    taxid = str(taxid)

    VA = request.args.get('VA', default=False)
    VA = bool(VA)

    return readids, subdir, taxid, VA


# Define the blueprint
blast_reads_blueprint = Blueprint('blast_reads', __name__)


@blast_reads_blueprint.route("/blast_reads")
def blast_reads():
    # get parameters, if cgi use form otherwise use argparse
    readids, subdir, taxid, VA = get_parameters_blast()
    print("Parameters: " + readids, subdir, taxid)
    fasta_string = create_query(readids, subdir, taxid, VA)
    # print("Fasta string: " + fasta_string)
    url = blast_query(fasta_string)
    # print('<meta http-equiv="refresh" content="0;URL=%s" />' % url)
    return redirect(url, code=302)
