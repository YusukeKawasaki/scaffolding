#coding: utf-8

import subprocess
import csv
import copy
import xlrd
import inspect
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from blasttype import *
from simple_fasta import * 
from contig import *
from sqs import *

QUERY_FASTA_NAME = "ans298_170124_0.fasta"
SUBJECT_FASTA_NAME = QUERY_FASTA_NAME
BLAST_NAME = "blast_%s_%s"%(QUERY_FASTA_NAME, SUBJECT_FASTA_NAME)
OUTPUT_FASTA_NAME = "new" + QUERY_FASTA_NAME

cmd = 'makeblastdb -in %s -dbtype nucl -hash_index'%SUBJECT_FASTA_NAME
subprocess.call( cmd, shell=True  )
cmd = 'blastn -db %s -query %s -out %s -outfmt "6 std qlen slen"'%(SUBJECT_FASTA_NAME, 
                                                                    QUERY_FASTA_NAME, BLAST_NAME)
subprocess.call(cmd, shell=True)

#blastの型、始点型が１、終点型が２、内包型が３、被内包型が４
START_LINK = 1
END_LINK = 2
CONTAIN = 3
CONTAINED = 4
OTHER = 0

BLAST_DATA = Blast_data(BLAST_NAME)
QUERY_FASTA_DATA = SeqIO.parse(QUERY_FASTA_NAME, "fasta")
SUBJECT_FASTA_DATA = SeqIO.parse(SUBJECT_FASTA_NAME,"fasta")
QUERY_FASTA = Simple_multi_fasta(QUERY_FASTA_DATA)
SUBJECT_FASTA = Simple_multi_fasta(SUBJECT_FASTA_DATA)

print 'data_roading done'

NNN_CHECK_FASTA = Simple_multi_fasta()
NNN_CHECK_FASTA.add_fasta(QUERY_FASTA)

rowMax = BLAST_DATA.nrows#最終行を整数値で取得
BITSCORE_BASELINE = 1000
well_blast_list=[]#ビットスコア一定以上の行を抽出する

#well_blast_listの作成
for row in range(rowMax):
    blast = Blast(BLAST_DATA, row)
    if blast.bitscore > BITSCORE_BASELINE:
        well_blast_list.append(blast)

print 'well_blast_list(bitscore) done'
#decideでblastの型を確定
for blast in well_blast_list:
    blast.decide()

need_blasts = []
for blast in well_blast_list:
    if blast.match_num > 10000 and not blast.is_equal():
        need_blasts.append(blast)

for blast in need_blasts:
    for a_blast in need_blasts:
        if a_blast.Is_same(blast):
            if blast.sstart < 10 or blast.send < 10 or \
            blast.sstart > blast.slen - 10 or blast.send > blast.slen - 10:
                need_blasts.remove(blast)
            else:
                need_blasts.remove(a_blast)

for blast in need_blasts:
    if blast.qstart < 10:
        query_fasta = QUERY_FASTA.search_seq(blast.qname)
        query_fasta.seq = query_fasta.seq[blast.qend - query_fasta.replace_dev:]
        query_fasta.replace_dev = blast.qend - 1
    elif blast.qend > blast.qlen - 10:
        query_fasta = QUERY_FASTA.search_seq(blast.qname)
        query_fasta.seq = query_fasta.seq[:blast.qstart - query_fasta.replace_dev]

s_seq_list = [x.make_seqrecord() for x in QUERY_FASTA.fasta]
SeqIO.write(s_seq_list, OUTPUT_FASTA_NAME, "fasta")