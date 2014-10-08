''' This file converts a gtf file into protobuf format.

See Gene defintion in proto/jeweler.proto.

'''
import collections
import re
import os
import copy
import pyfasta

import sys
# Add the ptdraft folder path to the sys.path list
sys.path.append('../python')

import rnasigs_pb2 as jp2
import z
import zleveldb

gid_match = re.compile('gene_id "([^"]*)"')
gname_match = re.compile('gene_name "([^"]*)"')
tid_match = re.compile('transcript_id "([^"]*)"')
tname_match = re.compile('transcript_name "([^"]*)"')
biotype_match = re.compile('gene_biotype "([^"]*)"')

''' Parse the required information out from the extra field.

Return None if gid, gname, tid, tname, and biotype entries do not exist in 'extra'
This is a temporory solution, need a more thorough test
'''
def parse_extra(extra):

    gid = None
    gname = None
    tid = None
    tname = None
    biotype = None

    if(gid_match.search(extra)):
        gid = gid_match.search(extra).group(1)

    if(gname_match.search(extra)):
        gname = gname_match.search(extra).group(1)

    if(tid_match.search(extra)):
        tid = tid_match.search(extra).group(1)

    if(tname_match.search(extra)):
        tname = tname_match.search(extra).group(1)

    if(biotype_match.search(extra)):
        biotype = biotype_match.search(extra).group(1)

    # fix the bug in release 77: some genes / transcrips do not have a name
    if(gid != None and gname == None):
        gname = gid

    if(tid != None and tname == None):
        tname = tid

    return gid, gname, tid, tname, biotype

def identity(exons):
    return exons

class GTF2ProbobufConverter:
    def __init__(self, gtf_file, fasta_file):
        self.gtf_file = gtf_file

        def map_key(key):
            return key.split(' ')[0]

        self.fa = pyfasta.Fasta(fasta_file)
        self.fasta_file = fasta_file
        self.mapkeys = dict()
        for k in self.fa.keys():
            self.mapkeys[map_key(k)] = k

    def convert(self, out, protein_coding_only=True):

        os.system("rm -rf " + out)
        db = zleveldb.ZLevelDB(out)

        fd = open(self.gtf_file, "r")
        g = None
        last_gene = ""
        exon_meta = collections.defaultdict(list)

        for l in fd:
            if l.startswith("#"):
                continue
            data = l.strip().split('\t')
            seqname, source, feature, start, end, score, strand, frame, extra = data
            gid, gname, tid, tname, biotype = parse_extra(extra)

            if feature != "exon":
                continue
            if protein_coding_only and biotype != "protein_coding":
                continue
            if not gid:
                continue
            if gid != last_gene:
                if g:
                    self.clean_gene(g, exon_meta)
                    if g.id:
                        db.put(g.id, g)
                g = jp2.Gene()
                g.id, g.name, g.seqname, g.strand = gid, gname, seqname, strand
                last_gene = gid
                exon_meta = collections.defaultdict(list)

            if len(g.transcripts) == 0 or tid != g.transcripts[-1].id:
                t = g.transcripts.add()
                t.id = tid
                t.name = tname

            z.assertion(strand != g.strand,
                        "Strands are not matched within a gene %s." % gname)

            ## start and end are inclusive in gtf file
            ## start and end are 1-based.
            ## set end to be exclusive
            exon_meta[t.id].append([int(start), int(end) + 1])

        if gid and g:
            db.put(gid, g)
        fd.close()

    def clean_gene(self, g, exon_meta):
        def ordered(l):
            return l

        if g.strand == "+":
            order = ordered
        else:
            order = reversed

        for t in g.transcripts:
            for start, end in order(exon_meta[t.id]):
                exon = t.exons.add()
                exon.start, exon.end = start, end

                if g.seqname in self.mapkeys:
                    g.has_seq = True
                    ## fasta files are 0-based
                    chrseq = self.fa[self.mapkeys[g.seqname]]
                    exon.seq = str(chrseq[(exon.start - 1):(exon.end - 1)])
                    t.seq += exon.seq
                else:
                    g.has_seq = False
                    print("cannot find sequence " + g.seqname + " in " +
                          self.fasta_file)

            rel_start = 0
            for e in t.exons:
                for v in e.variants:
                    nv = t.variants.add()
                    nv.abs_pos = v.abs_pos
                    nv.ref_chr = v.ref_chr
                    nv.all_chr = v.all_chr
                    nv.rel_pos = rel_start + v.rel_pos
                rel_start += e.end - e.start



def iter_gene_db(db, show_progress=True):
    db = zleveldb.ZLevelDB(db)
    data = jp2.Gene()
    for gidx, (gid, data) in enumerate(db.iter(jp2.Gene)):
        if show_progress and gidx % 1000 == 0:
            print gidx
        if data.has_seq:
            yield data

def generate_transcript_fasta(db, fasta_file):
    '''This function dump all transcript sequences into a fasta file
    without the variants'''
    fd = file(fasta_file, "w+")
    for gene in iter_gene_db(db):
        for t in gene.transcripts:
            fd.write(">" + t.id + "\n")
            fd.write("" + t.seq + "\n")
    fd.close()

def generate_gene_fasta(db, gene_file):
    '''This function dump all gene sequences into a fasta file
    without the variants'''
    fd = file(gene_file, "w+")
    for gene in iter_gene_db(db):
        fd.write(">" + gene.id)
        for t in gene.transcripts:
            fd.write("|" + t.id)
        fd.write("\n")
        fd.write("|".join([t.seq for t in gene.transcripts]))
        fd.write("\n")
    fd.close()


def generate(root, gtf_file, gtf_protein_only_file, fa_file):
    ## This function is called at the folder root
    protobuf_file = gtf_file + ".jeweler.protobuf"
    protobuf_file_protein_only = gtf_file + ".jeweler.protobuf.protein_only"
    transcript_fasta = "transcripts.fa"
    gene_fasta = "genes.fa"
    transcript_fasta_pc = "transcripts.protein_coding.fa"
    gene_fasta_pc = "genes.protein_coding.fa"
    os.system("pwd")
    g2pc = GTF2ProbobufConverter(gtf_protein_only_file, fa_file)
    g2pc.convert(protobuf_file_protein_only, protein_coding_only=True)
    generate_transcript_fasta(protobuf_file_protein_only, transcript_fasta_pc)
    generate_gene_fasta(protobuf_file_protein_only, gene_fasta_pc)

    g2pc = GTF2ProbobufConverter(gtf_file, fa_file)
    g2pc.convert(protobuf_file, protein_coding_only=False)
    generate_transcript_fasta(protobuf_file, transcript_fasta)
    generate_gene_fasta(protobuf_file, gene_fasta)

    return root + transcript_fasta, root + gene_fasta, root + transcript_fasta_pc, root + gene_fasta_pc
