from project import app
from download_pavian_data import controllers
import argparse


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
                               required=True,
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
            controllers.main(args, human=True)
        else:
            app.run(debug=True, host=args.host, port=5000)
