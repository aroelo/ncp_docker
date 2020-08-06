import json
import os

from project import app


def visualize_jbrowse(taxid, sub_dir_path):
    trackList = make_json(taxid, sub_dir_path)
    JbrowseBaseUrl = f'http://jbrowse{app.config["HOST_DOMAIN"]}'
    # return_url = JbrowseBaseUrl + '?data=/data/' + sub_dir_path[sub_dir_path.index('pavianfiles'):]
    return_url = JbrowseBaseUrl + '?data=/data/pavianfiles/' + os.path.basename(sub_dir_path)
    handle = open(os.path.join(sub_dir_path, [f for f in os.listdir(sub_dir_path)
                                              if f.endswith('.header_count.txt')][0]))
    loc = handle.readline().split(' ')[0]
    return return_url


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
    tracklist['tracks'].append({"category": "Consensus sequence",
                                "key": "Consensus sequence",
                                "label": "Consensus sequence",
                                "type": "SequenceTrack",
                                "storeClass": "JBrowse/Store/SeqFeature/IndexedFasta",
                                "urlTemplate": "%s.cons.fa" % prefix, \
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
    # add GFF3 file to trackList.json
    tracklist['tracks'].append({"category": "Annotation",
                                "key": "Genbank annotation",
                                "label": "Genbank annotation",
                                "type": "JBrowse/View/Track/CanvasFeatures",
                                "urlTemplate": "%s.gff3.gz" % prefix,
                                "style": {
                                    "_defaultHistScale": 4,
                                    "_defaultLabelScale": 30,
                                    "_defaultDescriptionScale": 120,
                                    # Comma-separated list of case-insensitive feature tags to use
                                    # for showing the feature's label.
                                    # The first one found will be used. Default 'name,id'.
                                    "label": "produt,id",
                                    # style→description	Comma-separated list of case-insensitive
                                    # feature tags to check for the feature's long description.
                                    # The first one found will be used. Default 'note,description'.
                                    # If blank no description is used.
                                    "description": "note, description"
                                },
                                })

    json_path = os.path.join(input_dir, 'trackList.json')
    with open(json_path, 'wt') as output_handle:
        json_raw_str = json.dumps(tracklist, indent=4)
        output_handle.write(json_raw_str)
    return 'trackList.json'
