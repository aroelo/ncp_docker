#from services.web.project.download_pavian_data.make_output_files import make_output
import sys

def parse_URL(URL:str):
    from urllib.parse import urlparse
    # example URL:
    #     http://webproxy.naktuinbouw.cloud/download_pavian_data?taxid=42680&sample=/pavian/in/20211102_20211027_FAN40338_Diagnostiek_15.INS_21_14366_Tomaat.ntwgs.pavian&action=jbrowse
    res = urlparse(URL)
    d = {_.split('=')[0]:_.split('=')[1] for _ in res.query.split('&')}
    return d['taxid'], d['sample'].split('/')[-1].replace('.pavian','')

def test_output(URL:str):
    # proj_name = '20210323_20210317_FAN24381_Diagnostiek10.INS_2028466_2_Cannabis.ntwgs'
    # d10_can2 = f'/6_db_ssd/pavian_output/2080458_{proj_name}'
    # taxid=2080458
    td, proj_name = parse_URL(URL)

    make_output(
        sub_dir_path=d10_can2,
        taxid=td,
        bam_in_path=f'{d20_can2}/{proj_name}.filtered_s.bam',
        bigwig_path=f'{d20_can2}/{proj_name}.filtered_s.bw',
        df_reads_path=f'{d20_can2}/{proj_name}.readsdf.pickle'
    )

if __name__ == '__main__':
    print(parse_URL(sys.argv[1]))
