#!/usr/bin/env python
#coding: utf-8

from datetime import datetime
import copy
import subprocess
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse
import sys

from blasttype import *
from simple_fasta import * 
from contig import *
from sqs import *

def command(argv):

    global toLog
    global toRemoveN
    global SUBJECT_FASTA_NAME
    global QUERY_FASTA_NAME
    global OUTPUT_FASTA_NAME

    toLog = False
    toRemoveN = False
    OUTPUT_FASTA_NAME = datetime.now().strftime('%m%d%H%M.fasta')
    parser = argparse.ArgumentParser()

    #required arguments
    parser.add_argument("-p", help="PacBio fasta name", required=True)
    parser.add_argument("-m", help="MiSeq fasta name", required=True)

    #not required arguments
    parser.add_argument("-o", help="Output file name(default = time.now)", required=False)

    #Boolean arguments
    parser.add_argument("-n", help="Removes Miseq N-gap", action = "store_true", required=False)
    parser.add_argument("-l", help="Creates a log file", action = "store_true", required=False)

    args = parser.parse_args(argv)

    toLog = args.l
    toRemoveN = args.n
    SUBJECT_FASTA_NAME = args.p
    QUERY_FASTA_NAME = args.m
    if args.o is not None:
        OUTPUT_FASTA_NAME = args.o

command(sys.argv[1:])

nnn_removed_name = "n_removed_%s_%s"%(QUERY_FASTA_NAME, SUBJECT_FASTA_NAME) #nnnを埋めたファイルを作る際はそのファイル名
BLAST_NAME = "blast_%s_%s"%(QUERY_FASTA_NAME, SUBJECT_FASTA_NAME) #ここにタブ形式で出力されたBLASTファイルが入る
S_BLAST_NAME = "blast_%s_%s"%(SUBJECT_FASTA_NAME, SUBJECT_FASTA_NAME) #サブジェクト同士のBLAST結果がここに入る

#################
# main- routine #
#################

#blast
Do_blast(SUBJECT_FASTA_NAME, QUERY_FASTA_NAME, BLAST_NAME)

if toLog:
    log_text_name = "%s_log.txt"%OUTPUT_FASTA_NAME #ログテキスト名
    log_text = open(log_text_name, "w")
else:
    log_text = None

#fastaの読み込み
QUERY_FASTA_DATA = SeqIO.parse(QUERY_FASTA_NAME, "fasta")
SUBJECT_FASTA_DATA = SeqIO.parse(SUBJECT_FASTA_NAME,"fasta")
QUERY_FASTA = Simple_multi_fasta(QUERY_FASTA_DATA)
SUBJECT_FASTA = Simple_multi_fasta(SUBJECT_FASTA_DATA)

print 'data_roading done'

def make_blast_list(BLAST_NAME):
    u"""マルチ形無し排除まで"""
    BLAST_DATA = Blast_data(BLAST_NAME)
    rowMax = BLAST_DATA.nrows#最終行を整数値で取得
    BITSCORE_BASELINE = 0
    well_blast_list=[]#入れ物

    #BITSCOREで選別
    for row in range(rowMax): #[0,1,2,....,rowMax -1]
        blast = Blast(BLAST_DATA, row)
        if blast.bitscore > BITSCORE_BASELINE:
            well_blast_list.append(blast)

    print 'well_blast_list(bitscore) done'

    #あるblastに対してマルチ型にあたるblast群を確定
    for blast in well_blast_list:
        blast.MULTI_LIST_DECIDE(well_blast_list)

    well_balst_list_include_multi = well_blast_list[:]
    well_blast_list_temp = well_blast_list[:]

    print 'well_blast_list(multi_decide) done'

    #マルチ型をMulti_blastクラスに変換する P1
    for blast in well_blast_list:
        if blast.multi_list:
            if blast == blast.multi_list[0]: #ひとつを代表させて
                for other_blast in blast.multi_list:
                    try:
                        well_blast_list_temp.remove(other_blast) #一旦全部削除し
                    except:
                        pass
                blast = Multi_blast(BLAST_DATA, blast, well_blast_list)
                if blast.match_num * 0.7 < blast.sum_length: #スカスカのものは採用しない
                    well_blast_list_temp.append(blast) #使えるものをひとつだけ加える

    well_blast_list = well_blast_list_temp

    print 'well_blast_list(multi_class) done'

    #一応重複を排除
    well_blast_list = list(set(well_blast_list))

    #decideでblastの型を確定
    for blast in well_blast_list:
        blast.decide()

    print 'blast.decide done'

    #型がないものを削除
    well_blast_list = [blast for blast in well_blast_list if blast.type]

    #あるcontigの組み合わせに対して(START_LINK, END_LINK)に関してはBLAST結果はひとつだけで良い
    #bitscoreが最も高いものだけ採用して他はすべて削除する
    temp = []
    for blast in well_blast_list:
        for a_blast in well_blast_list:
            if blast.qname == a_blast.qname and\
            blast.sname == a_blast.sname and\
            blast.type == (START_LINK, END_LINK) and\
            blast.row != a_blast.row and\
            blast.bitscore < a_blast.bitscore:
                temp.append(blast)
    temp = list(set(temp))
    for blast in temp:
        well_blast_list.remove(blast)

    return well_blast_list

well_blast_list = make_blast_list(BLAST_NAME)

#BLASTから得たクエリ、サブジェクトのデータを作成
qname_set = set(map(lambda x:x.qname,well_blast_list))
sname_set = set(map(lambda x:x.sname,well_blast_list))

query_list = [Query(qname, QUERY_FASTA, well_blast_list) for qname in qname_set]
subject_list = [Subject(sname, SUBJECT_FASTA, well_blast_list) for sname in sname_set]

if toRemoveN:
    #n_replaceを全てのMiSeqに対して実行
    for query in query_list:
        query.n_replace(dev = 30, query_list = query_list, subject_list = subject_list)

    #書き出し
    n_replaced_list = [x.fasta.make_seqrecord() for x in query_list]
    SeqIO.write(n_replaced_list, nnn_removed_name, "fasta")

    #BLASTかけなおし
    N_BLAST_NAME = "tes.csv"
    Do_blast(SUBJECT_FASTA_NAME, nnn_removed_name, N_BLAST_NAME)

    well_blast_list = make_blast_list(N_BLAST_NAME)

    #FASTAの再定義(SUBJECTも一応)
    QUERY_FASTA_DATA = SeqIO.parse(nnn_removed_name, "fasta")
    SUBJECT_FASTA_DATA = SeqIO.parse(SUBJECT_FASTA_NAME,"fasta")
    QUERY_FASTA = Simple_multi_fasta(QUERY_FASTA_DATA)
    SUBJECT_FASTA = Simple_multi_fasta(SUBJECT_FASTA_DATA)

    #contigの再定義
    qname_set = set(map(lambda x:x.qname,well_blast_list))
    sname_set = set(map(lambda x:x.sname,well_blast_list))

    query_list = [Query(qname,QUERY_FASTA,well_blast_list) for qname in qname_set]
    subject_list = [Subject(sname, SUBJECT_FASTA, well_blast_list) for sname in sname_set]

#SQSデータの作成
sqs_list_raw = [Sqs(query, well_blast_list) for query in query_list]
sqs_list = [sqs for sqs in sqs_list_raw if sqs.type]

print 'sqs making done'

#しりとり開始
shiritori = Sqs_chains('new', sqs_list)
shiritori.add_blast(well_blast_list)

print 'sqs_chain done'

#使っていないqueryを特定
query_list_set = set(query_list)
not_using_querys = query_list_set - shiritori.USING_QUERYS(query_list,subject_list)

#使っていないqueryでchainを作成
for query in not_using_querys:
    if query.start_link_list:
        shiritori.add_chain([max(query.start_link_list, key = lambda x:x.bitscore)])
    elif query.end_link_list:
        shiritori.add_chain([max(query.end_link_list, key = lambda x:x.bitscore)])
    elif query.contain_list:
        shiritori.add_chain([max(query.contain_list, key = lambda x:x.bitscore)])

print "sqs_chain(not using queries) done"

#内包してるものの置き換え
for contig in query_list:
    contig.contain_replace(query_list = query_list, subject_list = subject_list, length = 5000) #5000塩基以上のもので置き換え

#置き換えているもののfasta化
#s_seq_list = [x.fasta.make_seqrecord() for x in query_list]
#SeqIO.write(s_seq_list,"okikae_260_s_5000.fasta", "fasta")

#print "replace done"

#fastaファイルとして出力
complete_data = [sqs_chain.make_fasta_data(query_list, subject_list, log_text) for sqs_chain in shiritori.sqs_chains]
SeqIO.write(complete_data, OUTPUT_FASTA_NAME, "fasta")