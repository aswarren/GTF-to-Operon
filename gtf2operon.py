#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os, sys
from optparse import OptionParser
from BCBio import GFF

##This script is written around some basic assumptions regarding the format of Cufflinks output. (mostly around the 'transcript' tag)
##If cufflinks output changes substantially this script may break.
##It expects the GTF/GFF files involved to be annotations of prokaryotic organisms
##The Cufflinks GTF file representing prokarytic organisms exposes a common assumption that a transcript belongs to a gene
##To use GFF3 supported conventions http://www.sequenceontology.org/gff3.shtml we should change 'transcript' to 'operon' but for now transcript will remain!


def main():
    DoNotNest =set(['transcript','inferred_parent'])
    parser = OptionParser()
    (options,args) = parser.parse_args()
    if len(args) < 2:
        sys.stderr.write("Usage: gtf2operon.py denovo.gtf someannotation.gtf/gff\n")
        sys.exit()

    #GTF files coming from Cufflinks through BCBio will have inferred parents. Go down to the subfeature level to get transcripts
    #set transcripts to be features here. screw the rest.

    in_handle=open(args[0])
    records=[]
    for rec in GFF.parse(in_handle):
        records.append(rec)
    in_handle.close()

    sequence_lookup={}
   
    for rec in records:
        sequence_lookup[rec.id]=rec
        transcripts=[]
        for feat in rec.features:
            for sub in feat.sub_features:
                if sub.type == 'transcript':
                    try:
                        sub.qualifiers['ID'] = sub.qualifiers.pop('gene_id')
                        sub.qualifiers.pop('Parent')
                    except KeyError:
                        sys.stderr.write("Faulty assumption for "+str(sub))
                    transcripts.append(sub)
        rec.features=transcripts
        rec.features.sort(key=lambda x: x.location.start, reverse=False)

    #attempt to fit existing annotations to transcript structure
    #if there is no strandedness information this is somewhat arbitrary
    #the script will compare which direction has more features fitting
    #into the transcript and then make a bulk assignment
    in_handle = open(args[1])
    annotations=[]
    for rec in GFF.parse(in_handle):
        annotations.append(rec)
        #flatten list of features (supports only 2 levels deep)
        sf = [x for y in annotations[-1].features for x in y.sub_features]
        annotations[-1].features=annotations[-1].features+sf
        for x in annotations[-1].features: x.sub_features=[]
        annotations[-1].features.sort(key=lambda x: x.location.start, reverse=False)
    in_handle.close()
    
    for seq in annotations:
        if seq.id in sequence_lookup:
            a=0
            mark=0
            for transcript in sequence_lookup[seq.id].features:
                strand_buckets={None:[],1:[],-1:[]}
                while seq.features[a].location.start < transcript.location.start and a < len(seq.features):
                    a+=1
                mark = a
                while mark < len(seq.features) and seq.features[mark].location.start < transcript.location.end:
                    if transcript.location.__contains__(seq.features[mark].location.start) and transcript.location.__contains__(seq.features[mark].location.end):
                        if not seq.features[mark].type in DoNotNest:
                            strand_buckets[seq.features[mark].strand].append(seq.features[mark])
                    mark+=1
                to_add=[]
                if transcript.strand:
                    to_add = strand_buckets[transcript.strand]
                #else:
                #    win_len=0
                #    for k in strand_buckets:
                #        if len(strand_buckets[k]) > win_len:
                #            win_len = len(strand_buckets[k])
                #            to_add =strand_buckets[k]
                for t in to_add:
                    t.qualifiers['Parent']=[transcript.id]
                    transcript.sub_features.append(t)
    GFF.write(records, sys.stdout) 
    
    
if __name__ == "__main__":
    main()
