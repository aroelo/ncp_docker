Added functionality to pavian
================================

In the comparison tab of pavian there is some added functionality that will show up when hovering over a cell in the comparison table.

The added functionality is per taxid and should only be used at the 'species' level for the best results.

Download
---------
Download output files for the given taxid/sample combination.
Clicking this button downloads a package of multiple files that are zipped together:

- Fasta file with all reads in the sample that map to the taxid (.fasta)
- Bam file with all reads in the sample that map to the taxid (.sorted.bam)
- Fasta file with references of the taxid to which reads of the sample are mapped (.ref.fa)
- Index file of above fasta (.ref.fa.fai)

Jbrowse
-----------
Visualise alignments of reads for the given taxid/sample combination.
This will open Jbrowse externally to display several tracks.
As of now, there are 5 tracks available:

1. Alignment
2. SNPs/Coverage
3. Reference sequence
4. Consensus sequence
5. Coverage

1. Alignment.
This shows the alignment between reads and the selected reference.
The alignment track can maximally show 100 reads at any position, so a 'capped' bam file is used for this track as input.
The 'capped' bam file is created in such a way that it only contains reads so that at no position there is a coverage > 100.
Meanwhile, this subset of reads will still give a good idea of what the alignment looks like across the reference.

2. SNPs/Coverage.
The SNPs/coverage track is based on the alignment track (1.) and you can see here that the coverage is never higher than 100.


3. Reference sequence.
The sequence of the given reference.

4. Consensus sequence.
A consensus sequence based on the reads that are displayed in the alignment track (1.).
This is created by running samtools mpilup with the capped bam file, for each position, the most prevalent nucleotide or indel is then determined.
Considering the insertions and deletions in the consensus sequence, it might not line up perfectly with the reference sequence.

5. Coverage.
A coverage track displaying the actual coverage. It's based on all reads that are mapped to the given reference and isn't capped at a coverage of 100.

View reads and blast
----------------------
This will open a datatable externally that contains all reads and their alignment scores.
Every row represents one read, in case there are secondary mapping for the reads their alignment scores are displayed as well together with an 'error tag' if applicable.
This 'error tag' has been assigned during the classification pipeline with the filter_bam.py script.