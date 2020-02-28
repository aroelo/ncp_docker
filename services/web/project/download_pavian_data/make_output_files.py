import os
from ete3 import NCBITaxa
import tempfile
import subprocess
import pysam
import numpy
import pyBigWig
from project import app
from flask import render_template


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


def make_capped_coverage_bam(bam_in_path, bam_out_path, header_count_path, coverage_cap):
    header_count = open(header_count_path, 'r')
    ref_list = []
    for header in header_count:
        ref_list.append(header.split(' ')[0])

    # In case bam file is empty, check_sq = false
    bam_in = pysam.AlignmentFile(bam_in_path, "rb", check_sq=False)
    bam_out = pysam.AlignmentFile(bam_out_path, "wb", template=bam_in)
    bam_headers_dict = dict(zip(bam_in.header.references, bam_in.header.lengths))
    skip = False
    for contig in bam_headers_dict.keys():
        if contig not in ref_list:
            continue
        ref_length = bam_headers_dict[contig]
        bin_size = 1
        if ref_length >= 100000:
            bin_size = 50
        elif ref_length >= 500000:
            bin_size = 100
        coverage_array = numpy.zeros(shape=-(-ref_length//bin_size), dtype=int)
        notallowed_indices = []

        read_count = 0
        for read in bam_in.fetch(contig):
            binned_start = read.reference_start//bin_size
            binned_end = read.reference_end//bin_size
            if read_count < coverage_cap-1:
                read_count += 1
                bam_out.write(read)
                coverage_array[binned_start:binned_end] += 1
                continue

            if notallowed_indices:
                for idx in notallowed_indices[0]:
                    if binned_start <= idx <= binned_end:
                        skip = True
                        break
            if skip:
                skip = False
                continue

            coverage_array[binned_start:binned_end] += 1
            if numpy.amax(coverage_array) < coverage_cap:
                bam_out.write(read)
            else:
                notallowed_indices = numpy.nonzero(coverage_array >= coverage_cap)
                coverage_array[binned_start:binned_end] -= 1
    bam_in.close()
    bam_out.close()
    return None


def make_bigwig(bigwig_in_path, bigwig_out_path, header_count_path):
    """make bigwig from bam"""
    bw = pyBigWig.open(bigwig_in_path)
    bw_out = pyBigWig.open(bigwig_out_path, "w")
    header = []

    with open(header_count_path, "r") as header_count:
        for line in header_count:
            contig = line.split(' ')[0]
            length = bw.chroms(contig)
            if not length:
                continue
            header.append((contig, length))

    if not header:
        return None

    bw_out.addHeader(header)

    with open(header_count_path, "r") as header_count:
        for line in header_count:
            contig = line.split(' ')[0]
            length = bw.chroms(contig)
            if not length:
                continue
            values = bw.values(contig, 0, length)
            entry = (contig, values)
            bw_out.addEntries(entry[0], 1, values=entry[1], span=1, step=1, validate=False)

    bw_out.close()

    return None


def render_datatable_html(df_pickle, taxid_file, taxids, html_path, sub_dir_path):
    taxids_file = open(taxid_file, "r")
    taxid_list = taxids_file.read().splitlines()
    taxid_list = ','.join(taxid_list)

    html_file = open(html_path, "w")
    # script_dir = os.path.dirname(os.path.realpath(__file__))
    # reads_datatables = open(os.path.join(script_dir, "reads_datatable_blast_flask.html")).read().format(df_pickle,
    #                                                                                                     taxid_list,
    #                                                                                                     sub_dir_path,
    #                                                                                                     str(taxids[0]))

    # Context is necessary, because template is rendered outside a request.
    with app.app_context():
        reads_datatables = render_template('reads_datatable_blast.html', df_pickle=df_pickle, taxid_list=taxid_list, sub_dir_path=sub_dir_path, taxid=str(taxids[0]))
    html_file.write(reads_datatables)
    html_file.close()
    return None


def delete_tmp_files(*args):
    for tmp_file in args:
        if os.path.exists(tmp_file):
            os.remove(tmp_file)


def make_output(sub_dir_path, taxid, bam_in_path, bigwig_path, df_reads_path):
    # get all child taxids for given taxid.
    taxids = [taxid] + get_child_taxa(taxid)

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

    # Sort header_count and cap at 15 per nt/wgs
    nt_header_path = sub_dir_path + "/" + str(taxids[0]) + ".nt_header.txt"
    wgs_header_path = sub_dir_path + "/" + str(taxids[0]) + ".wgs_header.txt"
    cmd = "sort -nrk 2,2 %s | mawk '{if ($0 ~ /^C/){ " \
          "if (c<30){c+=1; C_count[c]=$0}} " \
          "else if($0 ~ /^W|^>W/){" \
          "if (w<30){w+=1; W_count[w]=$0}}} " \
          "END{total=w+c; if(total>30){total=30}; " \
          "if (w<15){c2=total-w} else{c2=c}; if (c<15){w2=total-c} else {w2=w}; " \
          "for (i=1; i<=w2; i++){print W_count[i]; split(W_count[i],x,\" \"); print \">\" x[1] > \"%s\"}; " \
          "for (i=1; i<=c2; i++){print C_count[i]; split(C_count[i],x,\" \"); print \">\" x[1]  > \"%s\"}}' > %s " \
          "&& mv %s %s" % (
              header_count_path, wgs_header_path, nt_header_path, header_count_tmp, header_count_tmp, header_count_path)
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
    print('creating capped bam..')
    make_capped_coverage_bam(bam_out_path, capped_bam_out_path, header_count_path, coverage_cap=100)
    print('finished')

    # TODO, can this be somehow piped from stdin?
    # Index capped bam file
    cmd = "samtools index %s" % capped_bam_out_path
    run_cmd(cmd, log_out)

    bigwig_out_path = bam_out_path.replace('.bam', '.bw')
    chromsizes_path = sub_dir_path + "/" + str(taxids[0]) + ".chromsizes.txt"
    bedgraph_path = sub_dir_path + "/" + str(taxids[0]) + ".bedGraph"

    consensus_path = sub_dir_path + "/" + str(taxids[0]) + ".cons.fa"
    # Check to create consensus & method to create bigwig
    # If longest reference is bigger than 1 million bp, don't create consensus
    cmd = f"sort -nrk2,2 {ref_out_path+'.fai'} | cut -f2 | head -n1"
    ps = subprocess.Popen(cmd, shell=True, executable='/bin/bash',
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = ps.communicate()
    stdout = stdout.decode().replace('\r', '\n')
    #TODO adjust tracklist.Json or make clear to users in other way that consensus isn't made, because file is too big.
    if int(stdout) > 1000000:
        # create empty cons file
        cmd = f"touch {consensus_path}; touch {consensus_path+'.fai'}"
        run_cmd(cmd, log_out)

        # Create bigwig track
        cmd = f"samtools view -@20 -H {bam_out_path} | mawk 'BEGIN{{FS=\"\tSN:|\tLN:\";OFS=\"\\t\"}} {{if ($1==\"@SQ\") print $2, $3}}' > {chromsizes_path}; bedtools genomecov -ibam {bam_out_path} -bg -split | LC_COLLATE=C sort -k1,1 -k2,2n > {bedgraph_path}; bedGraphToBigWig {bedgraph_path} {chromsizes_path} {bigwig_out_path}"
        run_cmd(cmd, log_out)
        pass
    else:
        # Create consensus sequence
        mpileup2fa = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'mpileup2fa.py')
        cmd = f"samtools mpileup {capped_bam_out_path} -Q0 -a | python {mpileup2fa} -o {consensus_path}"
        run_cmd(cmd, log_out)

        # Index consensus
        cmd = f"samtools faidx {consensus_path}"
        run_cmd(cmd, log_out)

        # Use length of longest sequence (reference vs consensus)
        index_tmp_path = sub_dir_path + "/" + str(taxids[0]) + ".ref.fa.fai.tmp"
        cmd = f"mawk 'BEGIN{{OFS=\"\\t\"}} {{if(NR==FNR){{_[$1]=$2;next}} {{if (_[$1]>$2) print $1,_[$1],$3,$4,$5; else print $1,$2,$3,$4,$5}}}}' {consensus_path+'.fai'} {ref_out_path+'.fai'} > {index_tmp_path} & mv {index_tmp_path} {ref_out_path+'.fai'}"
        run_cmd(cmd, log_out)

        # Create bigwig track
        print('creating bigwig..')
        make_bigwig(bigwig_path, bigwig_out_path, header_count_path)
        print('finished')

    # # Get read scores per taxid, run python script with taxid and df pickle.
    html_path = sub_dir_path + "/" + str(taxids[0]) + ".read_scores.html"
    print('creating datatable html..')
    render_datatable_html(df_reads_path, taxid_tmp_file.name, taxids, html_path, sub_dir_path)
    print('finished')

    running_log_file = os.path.join(sub_dir_path, str(taxids[0]) + "running.log")
    print('deleting tmp files..')
    delete_tmp_files(sam_out_path, nt_out_path, wgs_out_path, nt_header_path, wgs_header_path, chromsizes_path, bedgraph_path, running_log_file)
    print('finished')

    return taxid_tmp_file.name