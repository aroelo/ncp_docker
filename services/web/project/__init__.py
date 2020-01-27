import os
from ete3 import NCBITaxa
import time
import tempfile
import shutil
import subprocess
import pysam
import numpy
import pyBigWig
import json
import zipfile
import pandas as pd
import argparse
from project import NCBIWWWcustom

from werkzeug import secure_filename
from flask import (
    Flask,
    jsonify,
    send_file,
    send_from_directory,
    request,
    redirect,
    url_for,
    render_template,
    abort
)
from flask_sqlalchemy import SQLAlchemy


app = Flask(__name__)
app.config.from_object("project.config.Config")
db = SQLAlchemy(app)


class User(db.Model):
    __tablename__ = "users"

    id = db.Column(db.Integer, primary_key=True)
    email = db.Column(db.String(128), unique=True, nullable=False)
    active = db.Column(db.Boolean(), default=True, nullable=False)

    def __init__(self, email):
        self.email = email


class PavianInput(db.Model):
    """
CREATE TABLE pavian_data (
    file varchar(100),
    run varchar(100),
    sample varchar(100),
    nt varchar(10),
    date DATE,
    organism_taxid INT,
    organism_name varchar(100)
    support TINYINT(1)
    );
    """
    __tablename__ = "pavian_data"

    id = db.Column(db.Integer, primary_key=True)
    file = db.Column(db.String(100))
    run = db.Column(db.String(100))
    sample = db.Column(db.String(100))
    nt = db.Column(db.String(10))
    date = db.Column(db.Date)
    organism_taxid = db.Column(db.Integer)
    organism_name = db.Column(db.String(100))
    support = db.Column(db.Boolean)

    def __init__(self, file, run, sample, nt, date, organism_taxid, organism_name, support):
        self.file = file
        self.run = run
        self.sample = sample
        self.nt = nt
        self.date = date
        self.organism_taxid = organism_taxid
        self.organism_name = organism_name
        self.support = support


@app.route("/static/<path:filename>")
def staticfiles(filename):
    return send_from_directory(app.config["STATIC_FOLDER"], filename)


@app.route("/media/<path:filename>")
def mediafiles(filename):
    return send_from_directory(app.config["MEDIA_FOLDER"], filename)


@app.route("/upload", methods=["GET", "POST"])
def upload_file():
    if request.method == "POST":
        file = request.files["file"]
        filename = secure_filename(file.filename)
        file.save(os.path.join(app.config["MEDIA_FOLDER"], filename))
    return f"""
    <!doctype html>
    <title>upload new File</title>
    <form action="" method=post enctype=multipart/form-data>
      <p><input type=file name=file><input type=submit value=Upload>
    </form>
    """


@app.route("/")
def hello_world():
    return jsonify(hello=app.config["STRING"])


def get_parameters(human):
    if not human:
        taxid = request.args.get('taxid')
        taxid = int(taxid)

        pavian_file = request.args.get('sample')
        pavian_file = str(pavian_file)

        action = request.args.get('action')
        action = str(action)

    else:
        taxid = args.taxid
        pavian_file = args.pavian_file
        action = args.action

    return taxid, pavian_file, action


def get_child_taxa(taxid):
    """get child taxids using """
    ncbi = NCBITaxa(dbfile=app.config['TAXA_SQLITE'])
    child_taxids = ncbi.get_descendant_taxa(int(taxid),
                                            intermediate_nodes=True)
    return child_taxids


def create_taxid_tmp_file(taxids):
    taxid_tmp_file = tempfile.NamedTemporaryFile(prefix='taxids',
                                                 suffix='.txt',
                                                 dir='/tmp',
                                                 delete=False,
                                                 mode='wt')
    taxid_str = '\n'.join([str(s) for s in taxids])
    taxid_tmp_file.write('%s' % taxid_str)
    taxid_tmp_file.close()
    return taxid_tmp_file


def run_cmd(cmd, log_out):
    """
    :param cmd: Command in string format or list
    :param log_out:
    :return:
    """
    if isinstance(cmd, list):
        cmd = [" ".join(cmd)]

    ps = subprocess.Popen("time(" + cmd + ")", shell=True, executable='/bin/bash',
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = ps.communicate()
    stdout = stdout.decode().replace('\r', '\n')
    stderr = stderr.decode().replace('\r', '\n')
    print(cmd, stdout, stderr)
    if stdout:
        print('stdout:\n%s\n' % stdout)
        log_out.write('stdout:\n%s\n' % stdout)
    if stderr:
        log_out.write('stderr:\n%s\n' % stderr)
    log_out.write('finished:\t%s\n\n' % cmd)
    ps.wait()


def make_capped_coverage_bam(bam_in_path, bam_out_path, coverage_cap):
    # In case bam file is empty, check_sq = false
    bam_in = pysam.AlignmentFile(bam_in_path, "rb", check_sq=False)
    bam_out = pysam.AlignmentFile(bam_out_path, "wb", template=bam_in)
    bam_headers_dict = dict(zip(bam_in.header.references, bam_in.header.lengths))

    for contig in bam_headers_dict.keys():
        ref_length = bam_headers_dict[contig]
        coverage_array = numpy.zeros(shape=ref_length, dtype=int)
        for read in bam_in.fetch(contig):
            coverage_array[read.reference_start:read.reference_end] += 1

            if numpy.amax(coverage_array) < coverage_cap:
                bam_out.write(read)
            else:
                coverage_array[read.reference_start:read.reference_end] -= 1

    bam_in.close()
    bam_out.close()
    return None


def make_bigwig(bigwig_in_path, bigwig_out_path, header_count_path):
    """make bigwig from bam"""
    header_count = open(header_count_path, "r")
    bw = pyBigWig.open(bigwig_in_path)
    bw_out = pyBigWig.open(bigwig_out_path, "w")
    header = []
    entries = []

    for line in header_count:
        contig = line.split(' ')[0]
        length = bw.chroms(contig)
        if not length:
            continue
        header.append((contig, length))
        values = bw.values(contig, 0, length)
        entry = (contig, values)
        entries.append(entry)

    if not header:
        return None

    bw_out.addHeader(header)

    for entry in entries:
        bw_out.addEntries(entry[0], 1, values=entry[1], span=1, step=1)

    return None


def render_datatable_html(df_pickle, taxid_file, taxids, html_path, sub_dir_path):
    taxids_file = open(taxid_file, "r")
    taxid_list = taxids_file.read().splitlines()
    taxid_list = ','.join(taxid_list)

    html_file = open(html_path, "w")
    # TODO use render_template function from flask instead
    script_dir = os.path.dirname(os.path.realpath(__file__))
    reads_datatables = open(os.path.join(script_dir, "reads_datatable_blast_flask.html")).read().format(df_pickle,
                                                                                                        taxid_list,
                                                                                                        sub_dir_path,
                                                                                                        str(taxids[0]))

    html_file.write(reads_datatables)
    html_file.close()
    return None


def delete_tmp_files(*args):
    for tmp_file in args:
        if os.path.exists(tmp_file):
            os.remove(tmp_file)


def make_json(prefix, input_dir):
    """make Jbrowse style JSON listing tracks"""
    # get list of files
    file_list = os.listdir(input_dir)
    # set reference sequence
    tracklist = {'formatVersion': 1,
                 'refSeqs': '%s.ref.fa.fai' % prefix,
                 'tracks': []}
    # add reference sequence track to tracklist.json
    tracklist['tracks'].append({"category": "Reference sequence",
                                "key": "Reference sequence",
                                "label": "Reference sequence",
                                "type": "SequenceTrack",
                                "storeClass": "JBrowse/Store/SeqFeature/IndexedFasta",
                                "urlTemplate": "%s.ref.fa" % prefix, \
                                "refSeqOrder": "False"})
    # add bigwig track to trackList.json
    tracklist['tracks'].append({"category": "Sequence data",
                                "key": "Coverage",
                                "label": "Coverage",
                                "type": "JBrowse/View/Track/Wiggle/XYPlot",
                                "storeClass": "JBrowse/Store/SeqFeature/BigWig",
                                "autoscale": "local",
                                "urlTemplate": "%s.sorted.bw" % prefix
                                })
    # add BAM Sequence Coverage to trackList.json
    tracklist['tracks'].append({"category": "Sequence data",
                                "key": "Sequence reads (SNPs/Coverage)",
                                "label": "Sequence reads (SNPs/Coverage)",
                                "type": "JBrowse/View/Track/SNPCoverage",
                                "storeClass": "JBrowse/Store/SeqFeature/BAM",
                                "urlTemplate": "%s.sorted.capped.bam" % prefix,
                                "cacheMismatches": "True",
                                "chunkSizeLimit": "5000000"
                                })
    # add BAM Sequence Alignments to trackList.json
    tracklist['tracks'].append({"category": "Sequence data",
                                "key": "Sequence reads (Alignment)",
                                "label": "Sequence reads (Alignment)",
                                "type": "JBrowse/View/Track/Alignments2",
                                "storeClass": "JBrowse/Store/SeqFeature/BAM",
                                "urlTemplate": "%s.sorted.capped.bam" % prefix,
                                # add bigwig histogram option
                                "cacheMismatches": "True",
                                "chunkSizeLimit": "5000000"
                                })

    json_path = os.path.join(input_dir, 'trackList.json')
    with open(json_path, 'wt') as output_handle:
        json_raw_str = json.dumps(tracklist, indent=4)
        output_handle.write(json_raw_str)
    return 'trackList.json'


def visualize_jbrowse(taxid, sub_dir_path):
    trackList = make_json(taxid, sub_dir_path)
    JbrowseBaseUrl = "http://" + str(app.config['HOST_IP']) + ':' + str(app.config['JBROWSE_PORT'])
    # return_url = JbrowseBaseUrl + '?data=/data/' + sub_dir_path[sub_dir_path.index('pavianfiles'):]
    return_url = JbrowseBaseUrl + '?data=/data/pavianfiles/' + os.path.basename(sub_dir_path)
    handle = open(os.path.join(sub_dir_path, [f for f in os.listdir(sub_dir_path)
                                              if f.endswith('.header_count.txt')][0]))
    loc = handle.readline().split(' ')[0]
    return return_url


def make_output(sub_dir_path, taxids, bam_in_path, bigwig_path, df_reads_path):
    log_out = open(sub_dir_path + "/" + str(taxids[0]) + '.log', 'wt')
    # Write taxids to tmpfile
    taxid_tmp_file = create_taxid_tmp_file(taxids)

    # Create header_count, bam, fasta and fastq
    header_count_path = sub_dir_path + "/" + str(taxids[0]) + ".header_count.txt"
    header_count_tmp = sub_dir_path + "/" + str(taxids[0]) + ".header_count.tmp"
    bam_out_path = sub_dir_path + "/" + str(taxids[0]) + ".sorted.bam"
    sam_out_path = sub_dir_path + "/" + str(taxids[0]) + ".sorted.sam"
    fasta_out = sub_dir_path + "/" + str(taxids[0]) + ".fasta"

    #          {print ">"$1"\\n"$10 > "%s"; if(assigned_taxid==taxid[3]) print $0, header_count[$3]++}}}
    cmd = ["samtools view -@20 -F2820 %s" % bam_in_path,
           """| mawk 'NR==FNR{taxids[$1]} {split($3,taxid,"|")}
          {for (i=12;i<=NF;i++) if ($i ~/^ti:Z:/){assigned_taxid=substr($i,6); if(assigned_taxid in taxids)
          {print ">"$1"\\n"$10 > "%s"; print $0, header_count[$3]++}}}
          END {for (header in header_count) print header, header_count[header] > "%s"}' %s /dev/stdin > %s; """
           % (fasta_out, header_count_path, taxid_tmp_file.name, sam_out_path),
           "samtools view -H -@20 -F2820 %s" % bam_in_path,
           """| mawk 'NR==FNR{header_count[$1]} {if ($1=="@SQ"&&substr($2,4) in header_count) {print $0}
          else {if ($1!="@SQ"&&substr($1,1,1)=="@") print $0}}' %s /dev/stdin | cat - %s |
          samtools view -@10 -Sb - > %s; samtools index -@10 %s"""
           % (header_count_path, sam_out_path, bam_out_path, bam_out_path)]

    # cmd = '|'.join(cmd).replace('\n', '')
    cmd = ''.join(cmd).replace('\n', '')
    run_cmd(cmd, log_out)

    # Sort header_count
    nt_header_path = sub_dir_path + "/" + str(taxids[0]) + ".nt_header.txt"
    wgs_header_path = sub_dir_path + "/" + str(taxids[0]) + ".wgs_header.txt"
    cmd = "sort -nrk 2,2 %s | mawk '{if ($0 ~ /^C/){print \">\" $1 > \"%s\"; " \
          "if (c<15){c+=1; C_count[c]=$0}} " \
          "else if($0 ~ /^W|^>W/){print \">\" $1 > \"%s\";" \
          "if (w<15){w+=1; W_count[w]=$0}}} " \
          "END{for (i=1; i<=w; i++) print W_count[i]; " \
          "for (i=1; i<=c; i++) print C_count[i]}' > %s " \
          "&& mv %s %s" % (
              header_count_path, nt_header_path, wgs_header_path, header_count_tmp, header_count_tmp, header_count_path)

    run_cmd(cmd, log_out)

    # Get nt sequences
    nt_out_path = sub_dir_path + "/" + str(taxids[0]) + ".nt.fa"
    if os.path.exists(nt_header_path):
        blastdb_nt = app.config['BLASTDB_NT']
        cmd = "mawk 'BEGIN{FS=\"\\|\"} {print $4}' %s | blastdbcmd -db %s -entry_batch - " \
              "-target_only -outfmt '%%s' | paste -d '\\n' %s - > %s" \
              % (nt_header_path, blastdb_nt, nt_header_path, nt_out_path)
        run_cmd(cmd, log_out)

    # Get wgs sequences
    wgs_out_path = sub_dir_path + "/" + str(taxids[0]) + ".wgs.fa"
    if os.path.exists(wgs_header_path):
        # TODO use v5 database when present on /6_db_ssd
        # cmd = "mawk 'BEGIN{FS=\"\\|\"} {print $4}' %s | blastdbcmd -db /6_db_ssd/blast_DBv5/refseq_genomic_v5/refseq_genomic.fa " \
        #       "-entry_batch - -target_only -outfmt '%%s' | paste -d '\\n' %s - > %s" \
        #       % (wgs_header_path, wgs_header_path, wgs_out_path)

        blastdb_refseq = app.config['BLASTDB_REFSEQ']
        cmd = "mawk 'BEGIN{FS=\"\\|\"} {print $4}' %s | blastdbcmd -db %s " \
              "-entry_batch - -target_only -outfmt '%%s' | paste -d '\\n' %s - > %s" \
              % (wgs_header_path, blastdb_refseq, wgs_header_path, wgs_out_path)
        run_cmd(cmd, log_out)

    # Merge references together
    ref_out_path = sub_dir_path + "/" + str(taxids[0]) + ".ref.fa"
    cmd = "cat %s %s > %s" % (wgs_out_path, nt_out_path, ref_out_path)
    run_cmd(cmd, log_out)

    # Index reference
    cmd = "samtools faidx %s" % ref_out_path
    run_cmd(cmd, log_out)

    # Create bam file that has a cap on coverage so it can be visualised in Jbrowse.
    capped_bam_out_path = sub_dir_path + "/" + str(taxids[0]) + ".sorted.capped.bam"
    make_capped_coverage_bam(bam_out_path, capped_bam_out_path, coverage_cap=100)

    # TODO, can this be somehow piped from stdin?
    # Index capped bam file
    cmd = "samtools index %s" % capped_bam_out_path
    run_cmd(cmd, log_out)

    # Create bigwig track
    bigwig_out_path = bam_out_path.replace('.bam', '.bw')
    make_bigwig(bigwig_path, bigwig_out_path, header_count_path)

    # Get read scores per taxid, run python script with taxid and df pickle.
    html_path = sub_dir_path + "/" + str(taxids[0]) + ".read_scores.html"
    render_datatable_html(df_reads_path, taxid_tmp_file.name, taxids, html_path, sub_dir_path)

    running_log_file = os.path.join(sub_dir_path, str(taxids[0]) + "running.log")
    delete_tmp_files(sam_out_path, nt_out_path, wgs_out_path, nt_header_path, wgs_header_path, running_log_file)
    return None


@app.route("/download_pavian_data")
def main(human=False):
    """main function loop.
    Sample input that needs to be available:
       * sorted bam file (SAMPLE + ".filtered_s.bam")
       * bigwig track of sorted bam file (SAMPLE + ".filtered_s.bw")
       * reference file (SAMPLE + ".ref.fa")
       * index of reference (SAMPLE + ".ref.fa.fai")
       * pavian file (SAMPLE + ".pavian")
       * df pickle of reads&scores (SAMPLE + ".readsdf.pickle)
    where SAMPLE can be either sample.nt or sample.ntwgs
    """
    # get parameters, if app use request otherwise use argparse
    taxid, pavian_file, action = get_parameters(human)

    # get all child taxids for given taxid.
    taxids = [taxid] + get_child_taxa(taxid)
    main_dir_path = app.config['PAVIAN_OUT']

    sample = os.path.basename(pavian_file).strip('.pavian')
    support_files_path = os.path.join(os.path.dirname(pavian_file), "support_files", sample)
    bam_path = os.path.join(support_files_path, sample + ".filtered_s.bam")
    bigwig_path = os.path.join(support_files_path, sample + ".filtered_s.bw")
    df_reads_path = os.path.join(support_files_path, sample + ".readsdf.pickle")
    # Make sure that input files are available
    assert os.path.exists(bam_path)
    assert os.path.exists(bigwig_path)
    assert os.path.exists(df_reads_path)

    sub_dir_path = os.path.join(main_dir_path, str(taxid) + "_" + sample)
    running_log = os.path.join(sub_dir_path, str(taxid) + "running.log")

    # make sure that output files are generated
    output_files = [os.path.join(sub_dir_path, str(taxid) + ".fasta"),
                    os.path.join(sub_dir_path, str(taxid) + ".header_count.txt"),
                    os.path.join(sub_dir_path, str(taxid) + ".log"),
                    os.path.join(sub_dir_path, str(taxid) + ".read_scores.html"),
                    # os.path.join(sub_dir_path, str(taxid) + ".read_scores.json"),
                    os.path.join(sub_dir_path, str(taxid) + ".ref.fa"),
                    os.path.join(sub_dir_path, str(taxid) + ".ref.fa.fai"),
                    os.path.join(sub_dir_path, str(taxid) + ".sorted.bam"),
                    os.path.join(sub_dir_path, str(taxid) + ".sorted.bam.bai"),
                    os.path.join(sub_dir_path, str(taxid) + ".sorted.bw"),
                    os.path.join(sub_dir_path, str(taxid) + ".sorted.capped.bam"),
                    os.path.join(sub_dir_path, str(taxid) + ".sorted.capped.bam.bai")]
    if all([os.path.isfile(f) for f in output_files]):
        pass
    # If files are being made atm, wait until finished.
    elif os.path.exists(running_log):
        while os.path.exists(running_log):
            time.sleep(2)
        pass
    else:
        try:
            os.mkdir(sub_dir_path)
        except FileExistsError:
            shutil.rmtree(sub_dir_path)

        open(running_log, "a")
        try:
            make_output(sub_dir_path, taxids, bam_path, bigwig_path, df_reads_path)
        except Exception as e:
            print(str(e))
            shutil.rmtree(sub_dir_path)
            abort(500, e)


    if action == 'jbrowse':
        # get jbrowse url
        url = visualize_jbrowse(taxid, sub_dir_path)
        print('<meta http-equiv="refresh" content="0;URL=%s" />' % url)
        return redirect(url, code=302)

    elif action == 'viewreads':
        url = str(taxid) + ".read_scores.html"
        return send_from_directory(directory=sub_dir_path, filename=url)

    elif action == 'download':
        zip_file = tempfile.NamedTemporaryFile(prefix='zip_',
                                               suffix='.zip',
                                               dir='/tmp',
                                               delete=False,
                                               mode='wt')
        zipf = zipfile.ZipFile(zip_file.name, 'w', zipfile.ZIP_DEFLATED)
        zipf.write(os.path.join(sub_dir_path, str(taxid) + ".fasta"), arcname=str(taxid) + ".fasta")
        zipf.write(os.path.join(sub_dir_path, str(taxid) + ".sorted.bam"), arcname=str(taxid) + ".sorted.bam")
        zipf.write(os.path.join(sub_dir_path, str(taxid) + ".ref.fa"), arcname=str(taxid) + ".ref.fa")
        zipf.write(os.path.join(sub_dir_path, str(taxid) + ".ref.fa.fai"), arcname=str(taxid) + ".ref.fa.fai")
        zipf.close()

        return send_file(zip_file.name,
                         mimetype='zip',
                         attachment_filename=os.path.basename(sub_dir_path) + '.zip',
                         as_attachment=True)


@app.route("/get_json_data")
def create_json_data():
    """Create json from pandas dataframe of primary and secondary alignments of reads"""
    df_pickle = request.args.get('df_pickle')
    taxid_list = request.args.get('taxid_list')
    taxid_list = taxid_list.split(',')

    df = pd.read_pickle(df_pickle)
    df_taxid = df[df['taxid0'].isin(taxid_list)]
    json_taxid = df_taxid.to_json(orient='split')
    return json_taxid


def get_parameters_blast():
    readids = request.args.get('readids')
    readids = str(readids)

    subdir = request.args.get('subdir')
    subdir = str(subdir)

    taxid = request.args.get('taxid')
    taxid = str(taxid)
    return readids, subdir, taxid


def create_query(readids, subdir, taxid):
    fasta_file_path = os.path.join(subdir, taxid + ".fasta")

    cmd1 = "mawk '/%s/ {print $0; getline; print $0}' %s" % (readids.strip("|"), fasta_file_path)
    ps1 = subprocess.Popen(cmd1, shell=True, executable='/bin/bash',
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    ps1.wait()
    stdout, stderr = ps1.communicate()
    stdout = stdout.decode().replace('\r', '\n')
    return stdout


def blast_query(fasta_string):
    get_url = NCBIWWWcustom.qblast("blastn", "nt", fasta_string, format_type="HTML")
    return get_url


@app.route("/blast_reads")
def blast_reads():
    # get parameters, if cgi use form otherwise use argparse
    readids, subdir, taxid = get_parameters_blast()
    print("Parameters: " + readids, subdir, taxid)
    fasta_string = create_query(readids, subdir, taxid)
    # print("Fasta string: " + fasta_string)
    url = blast_query(fasta_string)
    # print('<meta http-equiv="refresh" content="0;URL=%s" />' % url)
    return redirect(url, code=302)


@app.errorhandler(500)
def taxid_files_not_created(e):
    return render_template('500.html', traceback=e), 500


@app.route('/favicon.ico')
def favicon():
    return send_from_directory(os.path.join(app.root_path, 'static'),
                               'favicon.ico', mimetype='image/vnd.microsoft.icon')

def parse_args():
    parser = argparse.ArgumentParser(
        description='get command line arguments')
    configuration = parser.add_argument_group('input', 'input parameters')
    configuration.add_argument('-taxid', type=int,
                               help='taxonomy id to retrieve')
    configuration.add_argument('-pavian_file', type=str,
                               help='sample from which taxid data is extracted')
    configuration.add_argument('-action', type=str, choices={'download', 'jbrowse', 'view_reads'},
                               help='which action will be executed (either download/jbrowse/view_reads)')

    configuration.add_argument('-host', type=str,
                               default='***REMOVED***',
                               help='which address will be used to host the app on')
    configuration.add_argument('-jbrowse_port', type=str,
                               default=5004,
                               help='port on which jbrowse is hosted')
    configuration.add_argument('-out_dir_path', type=str,
                               default="/8_expand_ssd/pavian_output",
                               help='output directory to store files')
    configuration.add_argument('-taxa_sqlite', type=str,
                               default="/5_workspace/repos/Nanopore-classification-pipeline/data/taxa.sqlite",
                               help='taxa sqlite file used to look up taxids')
    configuration.add_argument('-blastdb_nt', type=str,
                               default="/6_db_ssd/blast_DBv5/nt_v5/nt_v5",
                               help='nt blastdb to extract reference sequences from')
    configuration.add_argument('-blastdb_refseq', type=str,
                               default="/6_db_ssd/blast_DBv4/refseq_genomic/refseq_genomic",
                               help='refseq blastdb to extract reference sequences from')

    configuration.add_argument('-human',
                               choices=['main', 'app'],
                               help='Debug with pycharm, either run "main" or run "app"')

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = parse_args()
    if args.human:
        app.config['PAVIAN_OUT'] = args.out_dir_path
        app.config['BLASTDB_NT'] = args.blastdb_nt
        app.config['BLASTDB_REFSEQ'] = args.blastdb_refseq
        app.config['HOST_IP'] = args.host
        app.config['JBROWSE_PORT'] = args.jbrowse_port
        app.config['TAXA_SQLITE'] = args.taxa_sqlite
        if args.human == 'main':
            main(human=True)
        else:
            app.run(debug=True, host=args.host, port=5000)
