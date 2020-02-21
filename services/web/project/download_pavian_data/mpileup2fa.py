import sys
from argparse import ArgumentParser


def argument_parser():
    parser = ArgumentParser()
    # use stdin if it's full
    if not sys.stdin.isatty():
        input_stream = sys.stdin

    # otherwise, read the given filename
    else:
        input_stream = None
        parser.add_argument(
            '-i',
            help='output from mpileup')

    parser.add_argument(
        '-o',
        help='output file')
    return input_stream, parser.parse_args()


def mpileup_to_fasta(mpileup, output):
    prev_contig = None
    consensus_seq = ''
    fasta = ''
    for line in mpileup:
        contig, coord, ref_base, no_reads, read_bases, base_quals = line.strip('\n').split('\t')
        # convert low quality bases to uppercase, so they are all the same
        read_bases = read_bases.upper()
        # add header & consensus to fasta if new contig
        if contig != prev_contig and prev_contig is not None:
            fasta += '>' + prev_contig + '\n' + consensus_seq + '\n'
            consensus_seq = ''

        if no_reads == '0':
            consensus_seq += 'N'
        counter_dict, insert_dict = counter(read_bases, ref_base, int(no_reads))
        base = counter_dict[0][0]
        # if deletion, don't add anything to consensus, else add most common base
        if base != '*':
            consensus_seq += base

        # if insertion add respective inserted bases
        insert = insert_dict[0][0]
        if insert.startswith('+'):
            consensus_seq += insert[2:]

        prev_contig = contig

    # for last line add header & consensus to fasta
    fasta += '>' + prev_contig + '\n' + consensus_seq + '\n'

    # write/print output
    if output:
        output_file = open(output, 'w')
        output_file.write(fasta)
    else:
        print(fasta)


def counter(read_bases, ref_base, no_reads):
    """
    Counts occurrence per base or indel in read_bases, stores in dictionary and returns sorted dict.
    See http://www.htslib.org/doc/samtools-mpileup.html for pileup format.
    :param read_bases: string with pileup info
    :return: sorted dictionaries with how many occurrences per base or indel
    """
    # add ref_base at start, if counts are the same for multiple bases, ref_base comes first
    counter_list = [ref_base, 'A', 'C', 'G', 'T', '*', 'N']
    counter_dict = {i: 0 for i in counter_list}
    insert_dict = {}
    indel_count = 0
    it = enumerate(read_bases)
    for base in it:
        # if read base is deletion or insertion
        if base[1] == '-' or base[1] == '+':
            indel_count += 1
            # index in read_base string
            indel_idx = base[0]
            # get length of deletion/insertion
            indel_len = ''
            for x in read_bases[indel_idx + 1:]:
                if x.isdigit():
                    indel_len += x
                else:
                    break
            # if insertion, get how many & which bases are inserted. deletions can be determined per base with '*'
            if base[1] == '+':
                insert_bases = read_bases[indel_idx + 2: indel_idx + 2 + int(indel_len)]
                try:
                    insert_dict[base[1]+indel_len + insert_bases] += 1
                except KeyError:
                    insert_dict[base[1]+indel_len + insert_bases] = 1
            # skip the next characters in read_base, depending on length of deletion/insertion
            [next(it, None) for _ in range(int(indel_len)+1)]

        else:
            # match with forward or reverse
            if base[1] == '.' or base[1] == ',':
                counter_dict[ref_base] += 1
            # reference skip with forward or reverse
            elif base[1] == '>' or base[1] == '<':
                counter_dict['N'] += 1
            # first position of read, followed by mapq symbol
            elif base[1] == '^':
                [next(it, None) for _ in range(1)]
            # last position of read
            elif base[1] == '$':
                pass
            # this can be deletions '*' or any nucleotide bases
            else:
                counter_dict[base[1]] += 1
    counter_dict = sorted(counter_dict.items(), key=lambda item: item[1], reverse=True)
    # check how many reads do not have a deletion/insertion. This is necessary to determine if there are more
    # insertions than 'other' assignments.
    insert_dict['no_indel'] = no_reads - indel_count
    insert_dict = sorted(insert_dict.items(), key=lambda item: item[1], reverse=True)
    return counter_dict, insert_dict


def main():
    """main function loop"""
    input_stream, args = argument_parser()
    if not input_stream:
        input_stream = open(args.i, 'r')
    mpileup_to_fasta(input_stream, args.o)


if __name__ == '__main__':
    main()
