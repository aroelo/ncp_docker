import os
import time
import tempfile
import shutil
import zipfile
from pathlib import Path

from project import app
from project.download_pavian_data.make_output_files import make_output
from project.jbrowse.visualise_jbrowse import visualize_jbrowse
import traceback

from flask import (
    Blueprint,
    send_file,
    send_from_directory,
    request,
    redirect,
    abort
)


def get_parameters(human, args=None):
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


# Define the blueprint
download_pavian_data_blueprint = Blueprint('download_data', __name__)


@download_pavian_data_blueprint.route("/download_pavian_data")
def main(args=None, human=False):
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
    taxid, pavian_file, action = get_parameters(human, args)

    main_dir_path = app.config['PAVIAN_OUT']
    sample = Path(pavian_file).with_suffix('')

    support_files_path = Path(pavian_file).parent / "support_files" / sample
    bam_path = support_files_path / f"{sample}.filtered_s.bam"
    bigwig_path = support_files_path / f"{sample}.filtered_s.bw"
    df_reads_path = support_files_path / f"{sample}.readsdf.pickle"
    # Make sure that input files are available
    try:
        assert bam_path.exists(), f'Bam path did not exist. Given:\n{bam_path}'
        assert bigwig_path.exists(), f'Bigwig path did not exist. Given:\n{bigwig_path}'
        assert df_reads_path.exists(), f'DF reads path did not exist. Given:\n{df_reads_path}'
    except Exception:
        print(traceback.format_exc())
        e = ('item:' + str(taxid) + 'item:' + str(sample) + 'item:' + traceback.format_exc() + '\n' +
             '\n'.join([f for f in [bam_path, bigwig_path, df_reads_path] if not f.exists()]))
        abort(500, e)

    sub_dir_path = os.path.join(main_dir_path, str(taxid) + "_" + sample)
    running_log_path = os.path.join(sub_dir_path, str(taxid) + "running.log")

    # make sure that output files are generated
    output_files = [os.path.join(sub_dir_path, str(taxid) + ".fasta"),
                    os.path.join(sub_dir_path, str(taxid) + ".header_count.txt"),
                    os.path.join(sub_dir_path, str(taxid) + ".log"),
                    os.path.join(sub_dir_path, str(taxid) + ".read_scores.html"),
                    # os.path.join(sub_dir_path, str(taxid) + ".read_scores.json"),
                    os.path.join(sub_dir_path, str(taxid) + ".ref.fa"),
                    os.path.join(sub_dir_path, str(taxid) + ".ref.fa.fai"),
                    os.path.join(sub_dir_path, str(taxid) + ".cons.fa"),
                    os.path.join(sub_dir_path, str(taxid) + ".cons.fa.fai"),
                    os.path.join(sub_dir_path, str(taxid) + ".sorted.bam"),
                    os.path.join(sub_dir_path, str(taxid) + ".sorted.bam.bai"),
                    os.path.join(sub_dir_path, str(taxid) + ".sorted.bw"),
                    os.path.join(sub_dir_path, str(taxid) + ".sorted.capped.bam"),
                    os.path.join(sub_dir_path, str(taxid) + ".sorted.capped.bam.bai")]
    if all([os.path.isfile(f) for f in output_files]):
        pass
    # If files are being made atm, wait until finished.
    elif os.path.exists(running_log_path):
        print('Already creating files, waiting until finished')
        while os.path.exists(running_log_path):
            time.sleep(2)
        pass
    else:
        try:
            os.mkdir(sub_dir_path)
        except FileExistsError:
            shutil.rmtree(sub_dir_path)
            os.mkdir(sub_dir_path)

        open(running_log_path, "a")
        try:
            taxid_tmp_file = make_output(sub_dir_path, taxid, bam_path, bigwig_path, df_reads_path)
        except Exception:
            print(traceback.format_exc())
            shutil.rmtree(sub_dir_path)  # fixme: Consider moving this to an "error-spot" so we can review problems
            e = 'item:' + str(taxid) + 'item:' + str(sample) + 'item:' + traceback.format_exc()
            abort(500, e)

    if action == 'jbrowse':
        # get jbrowse url
        try:
            url = visualize_jbrowse(taxid, sub_dir_path)
        except Exception:
            print(traceback.format_exc())
            e = 'item:' + str(taxid) + 'item:' + str(sample) + 'item:' + traceback.format_exc()
            abort(500, e)
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
