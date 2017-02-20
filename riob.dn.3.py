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
    parser.add_argument("-o", help="Output file name(default = time)", required=False)

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

#blastの型、始点型が１、終点型が２、内包型が３、被内包型が４
START_LINK = 1
END_LINK = 2
CONTAIN = 3
CONTAINED = 4
OTHER = 0

#クエリかサブジェクトか
QUERY = 1
SUBJECT = 2

#プラスかマイナスか
PLUS = 1
MINUS = 2
direction_mark = {PLUS:'+',MINUS:'-'}

#Sqs用
START = 0
END = 1

def Do_blast(SUBJECT_FASTA_NAME, QUERY_FASTA_NAME, BLAST_NAME):
    cmd = 'makeblastdb -in %s -dbtype nucl -hash_index'%SUBJECT_FASTA_NAME #BLAST用のデータベースを作るコマンド
    subprocess.call( cmd, shell=True  ) #コマンドをターミナルで出力する
    cmd = 'blastn -db %s -query %s -out %s -outfmt "6 std qlen slen"'%(SUBJECT_FASTA_NAME, 
                                                                        QUERY_FASTA_NAME, BLAST_NAME) #BLASTデータを作成するコマンド
    subprocess.call(cmd, shell=True) #コマンドをターミナルで出力する

class Blast_data:
    "以上の操作で得られたBLASTのTABファイルを(x, y)の形で取得するためのクラス"
    def __init__(self, tab_text_name):
        f = open(tab_text_name)
        self.lines = f.readlines() #一行ごとに読み込み
        f.close()
        self.nrows = len(self.lines) #行数取得
        self.list = [line.split() for line in self.lines] #TAB位置で分割して配列化することで

    def cell(self, row, col):
        return self.list[row][col] #(x, y)の形で取得できる
class Simple_multi_fasta:
    u"""SeqIOで読んだマルチfastaファイルを指定すると
    配列と名前だけを持ったシンプルな配列情報を持つオブジェクト(=simple_fasta)を生成
    self.search_seq(コンティグ名)で、所持しているsimple_fastaオブジェクトを検索可能"""
    def __init__(self, fasta_data = None):
        if fasta_data:
            self.fasta = self.seq_loading(fasta_data) #fastaデータを読み込みSimple_fastaのリストの形で返す
        else:
            self.fasta = []
    
    def seq_loading(self, fasta_data):
        return [Simple_fasta(data) for data in fasta_data]

    def search_seq(self,name):
        u"""contig名を入れるとfastaデータを返す
        得られたfastaからはfasta.seqで配列を取得可能"""
        for fasta in self.fasta:
            if fasta.name == name:
                return fasta

    def add_fasta(self,fasta):
        u"""蛇足。fastaを追加するだけ"""
        if isinstance(self, Simple_fasta):
            self.fasta.append(fasta)
        elif isinstance(self, Simple_multi_fasta):
            self.fasta += fasta.fasta
        else:
            print 'error! fasta is not usable for add_fasta func'

class Simple_fasta:
    u"""Simple_multi_fastaが持つオブジェクトでbiopythonモジュールのseq_recordと対応
    self.seq == 配列データ
    self.make_seqrecord()でseq_recordに変換"""
    def __init__(self, fasta_data):
        u"""fasta_data には seq_recordを入れる
        必要な情報をseq_recordからコピー"""
        self.name = fasta_data.name
        self.seq = fasta_data.seq
        self.id = fasta_data.id
        self.description = fasta_data.description

    def make_seqrecord(self):
        u"""Simple_fastaが自身をseqrecordに変換する
        seqrecordに戻さないとSeqIOでの書き出しができない"""
        seq_r = SeqRecord(self.seq)
        seq_r.name = self.name
        seq_r.id = self.id
        seq_r.description = self.description
        return seq_r

class Blast:
    u"""blastのデータ、指定された行の情報を持つ
    row:(行数)(半角スペース)の形で検索可能(search_object関数)"""
    def __init__(self, BLAST_DATA, row):
        u"""BLAST_DATAオブジェクトから情報を取り出す"""
        self.row = row
        self.qname = str(BLAST_DATA.cell(row, 0))
        self.sname = str(BLAST_DATA.cell(row, 1))
        self.name = 'row:' + str(self.row + 1) + ' ' + self.qname + ', ' + self.sname
        self.match_num = int(BLAST_DATA.cell(row, 3))
        self.qstart = int(BLAST_DATA.cell(row, 6)) - 1
        self.qend = int(BLAST_DATA.cell(row, 7)) - 1
        self.sstart = int(BLAST_DATA.cell(row, 8)) - 1
        self.send = int(BLAST_DATA.cell(row, 9)) - 1
        self.bitscore = float(BLAST_DATA.cell(row, 11))
        self.qlen = int(BLAST_DATA.cell(row, 12))
        self.slen = int(BLAST_DATA.cell(row, 13))
        self.is_reverse = bool(self.sstart - self.send > 0) #クエリに対してサブジェクトが相補鎖かどうか

    def decide(self):
        u"""始点型などの判定で用いる値、真偽値を、マルチ判定を終えて初期情報が確定してから計算する。
        型(=self.type)を決定する。
        START_LINK = 1, END_LINK = 2, CONTAIN = 3, CONTAINED = 4"""
        TYPE_CONDITION_RATE = 0.2
        self.allowable_edge_error = self.match_num * TYPE_CONDITION_RATE #許容誤差
        self.query_subject_displacement = self.qstart - self.sstart #蛇足　クエリとサブジェクトのずれ
        #以下4つは一致部分から端までの距離。条件5,6で用いる値
        self.qstart_to_edge_length = self.qstart
        self.qend_to_edge_length = self.qlen - self.qend + 1
        self.sstart_to_edge_length = self.SSTART_TO_EDGE_LENGTH() #PacBioは向きによって処理が異なるので関数に回す
        self.send_to_edge_length = self.SEND_TO_EDGE_LENGTH() #同様

        self.is_qstart_exists_to_edge = \
        bool(self.allowable_edge_error > self.qstart_to_edge_length) #条件1

        self.is_qend_exists_to_edge = \
        bool(self.allowable_edge_error > self.qend_to_edge_length) #条件2

        self.is_sstart_exists_to_edge = \
        bool(self.allowable_edge_error > self.sstart_to_edge_length) #条件3

        self.is_send_exists_to_edge = \
        bool(self.allowable_edge_error > self.send_to_edge_length) #条件4

        self.is_query_long_in_start = \
        bool(self.qstart_to_edge_length > self.sstart_to_edge_length) #条件5

        self.is_query_long_in_end = \
        bool(self.qend_to_edge_length > self.send_to_edge_length) #条件6
        
        self.type = self.WHICH_TYPE() #始点型等の型。型ごとに対応する数値はdecideのコメントを参照
        self.subject_using_nucl_edge_nums = self.SUBJECT_USING_NUCL_EDGE_NUMS() #蛇足

    def Is_connectable(self,blast, contain_type = None):
        u"""blastとblastがつながるかどうか。真偽値を返す。
        blastがcontainの場合、start_link, end_linkのどちらとして扱うかを
        contain_typeに入力する"""
        if self.sname == blast.sname:
            if self.type == CONTAIN:
                type_ = contain_type #contain_type != CONTAIN に注意
            else:
                type_ = self.type
            if (type_ + blast.type + self.is_reverse + blast.is_reverse)%2 == 1: #排他的論理和マジックナンバー方式
                return True
            else:
                return False
        else:
            False

    def IS_MULTI_WITH(self, blast):
        u"""緩いマルチ型判定をおこなう
        MULTI_LIST_DECIDEで最終判定を行う前段階の関数"""
        MULTI_CONDITION_DEV = 2000
        if self.qname == blast.qname and self.sname == blast.sname\
        and self.is_reverse == blast.is_reverse\
        and abs( abs(self.qstart - blast.qstart) - abs(self.sstart - blast.sstart) )\
        < MULTI_CONDITION_DEV: #塩基数のずれが一定以下なら一応マルチ型とみなす
            return True
        else:
            return False

    def MULTI_LIST_DECIDE(self,blast_list):
        u"""自分とマルチ型の関係にあるblast(自分含む)をピックアップしたリストをself.multi_listとする
        そのようなリストが存在しない場合空のリストをself.multi_listとする"""
        l = [] #まずlに候補を入れ、その中で条件を満たすものをwell_multi_listに突っ込む
        for blast in blast_list:
            if self.IS_MULTI_WITH(blast): #is_Multi_with = ゆるいマルチ型判定を行う
                l.append(blast)
        l = list(set(l))
        l.sort(key = lambda x:x.row)
        well_multi_blasts = [self] #最終的に出力するリスト
        l.remove(self)
        Is_wrong = False
        for blast in l:
            for well_blast in well_multi_blasts:
                if well_blast.qstart < blast.qstart <  well_blast.qend or\
                well_blast.qstart < blast.qend<  well_blast.qend or\
                blast.qstart < well_blast.qstart < blast.qend:#この条件改良の余地あり これを満たすものが最終的にself.multi_listに入る
                    Is_wrong = True #Wrongなものはif文に回してappendしないようにしている
                    break
            if Is_wrong: 
                Is_wrong = False
                continue
            else:
                well_multi_blasts.append(blast) #Wrongでないものがここに来る
        well_multi_blasts.sort(key = lambda x:x.row)
        #最終的にreturnは行わずself.multi_listに格納する
        if len(well_multi_blasts) > 1:
            self.multi_list = well_multi_blasts
        else:
            self.multi_list = []

    def SAY_NAME(self):
        print self.row, '|' ,self.qname, self.sname, '|'

    def SSTART_TO_EDGE_LENGTH(self):
        u"""条件5用の値"""
        if self.is_reverse:
            return self.slen - self.sstart + 1
        else:
            return self.sstart

    def SEND_TO_EDGE_LENGTH(self):
        u"""条件6用の値"""
        if self.is_reverse:
            return self.send
        else:
            return self.slen - self.send + 1

    def WHICH_TYPE(self):
        """blastを4つの型に分類し、型に対応した1～4の整数を返す
        非分類は0を返す。型別のフローチャート参照"""
        if self.is_qstart_exists_to_edge and self.is_sstart_exists_to_edge:
            if self.is_query_long_in_start:
                self.is_qstart_exists_to_edge = False
                self.is_sstart_exists_to_edge = True
            else:
                self.is_qstart_exists_to_edge = True
                self.is_sstart_exists_to_edge = False

        if self.is_qend_exists_to_edge and self.is_send_exists_to_edge:
            if self.is_query_long_in_end:
                self.is_qend_exists_to_edge = False
                self.is_send_exists_to_edge = True
            else:
                self.is_qend_exists_to_edge = True
                self.is_send_exists_to_edge = False

        if self.is_qstart_exists_to_edge and self.is_send_exists_to_edge:
            return START_LINK
        elif self.is_sstart_exists_to_edge and self.is_qend_exists_to_edge:
            return END_LINK
        elif self.is_qstart_exists_to_edge and self.is_qend_exists_to_edge:
            return CONTAIN
        elif self.is_sstart_exists_to_edge and self.is_send_exists_to_edge:
            return CONTAINED
        else:
            return OTHER

    def SUBJECT_USING_NUCL_EDGE_NUMS(self):
        u"""クエリの末端の塩基が、サブジェクトでどの位置の塩基に対応するかを表す
        SubjectのCONTAIN_WELL_LIST関数で使用"""
        start_edge = None 
        end_edge = None
        if self.type in (START_LINK, CONTAIN):
            if not self.is_reverse:
                start_edge = self.sstart - self.qstart_to_edge_length
            else:
                start_edge =  self.sstart + self.qstart_to_edge_length
                
        if self.type in (END_LINK, CONTAIN):
            if not self.is_reverse:
                end_edge = self.send + self.qend_to_edge_length
            else:
                end_edge = self.send - self.qend_to_edge_length
        
        return (start_edge, end_edge)

    

class Multi_blast(Blast):
    """自分とマルチ関係にあるblastオブジェクトをひとつにまとめたオブジェクト
    初期値はひとつのblastデータを代表させて入れる"""
    def __init__(self,BLAST_DATA,blast,blast_list):
        Blast.__init__(self,BLAST_DATA, blast.row)
        self.multi_list = blast.multi_list #multi_listは__init__で定義していないのでここで定義
        self.qstart = min(map(lambda x:x.qstart, self.multi_list))
        self.qend = max(map(lambda x:x.qend, self.multi_list))
        self.start_of_blast = self.multi_list[ map(lambda x:x.qstart, self.multi_list).index(self.qstart) ] #マルチ型の中で一番はじめのBLAST
        self.end_of_blast = self.multi_list[ map(lambda x:x.qend, self.multi_list).index(self.qend) ] ##マルチ型の中で一番終わりのBLAST
        self.sstart = self.start_of_blast.sstart 
        self.send = self.end_of_blast.send
        self.match_num = self.qend - self.qstart + 1 #見かけ上の一致部位
        #以下、追加項目
        self.multi_null_positions = self.MULTI_NULL_POSITIONS() #マルチ同士の間に空いた隙間のpositionを示す
        self.sum_length = sum([x.match_num for x in self.multi_list]) #一致部分の塩基長の総和 スカスカだとself.match_numよりかなり小さくなる

    def MULTI_NULL_POSITIONS(self):
        """マルチ同士の間に空いた隙間のpositionを返す
        N埋めの際にNに該当する部分の配列を指す"""
        multi_list = self.multi_list[:]
        multi_null_positions = [] #returnするリスト
        multi_qstart = self.start_of_blast.qend #マルチ型間の隙間のpositionの始点 query側の塩基番号 以下同様
        multi_sstart = self.start_of_blast.send #間隙なのでstartが終点にendが始点になることに注意
        multi_qend = 0
        while multi_list:
            min_num_blast = min(multi_list, key = lambda blast: blast.qstart)
            if min_num_blast.qstart < multi_qstart: #blast同士が重複するときこのパターンはマルチ型確定の際に消しているので蛇足
                multi_list.remove(min_num_blast)
            else:
                multi_qend = min_num_blast.qstart 
                multi_send = min_num_blast.sstart
                multi_null_positions.append(((self.qname, multi_qstart, multi_qend),(self.sname, multi_sstart, multi_send)))
                #上のような形式でリストに入り、returnされる
                multi_qstart = min_num_blast.qend
                multi_sstart = min_num_blast.send
                multi_list.remove(min_num_blast)
        return multi_null_positions


class contig:
    PLUS = 1
    MINUS = 0
    u"""クエリ、サブジェクト双方の親クラス"""
    def __init__(self,name,blast_list):
        self.start_link_list = self.HAVING_BLAST(blast_list, START_LINK, self.type)
        self.end_link_list = self.HAVING_BLAST(blast_list, END_LINK, self.type)
        self.contain_list = self.HAVING_BLAST(blast_list, CONTAIN, self.type)
        self.contained_list = self.HAVING_BLAST(blast_list, CONTAINED, self.type)
        self.other_list = self.HAVING_BLAST(blast_list, OTHER, self.type)
        self.multi_null_positions = self.MULTI_NULL_POSITIONS() #自分を含むマルチ型の空白部分全て
        self.replace_dev = 0

    def HAVING_BLAST(self ,blast_list, blast_type ,contig_type):
        u"""クエリ、サブジェクト双方の親クラス
        例えばblast_type = START_LINKなら、自分を含む始点型のBLASTを全てreturnする"""
        temp = []
        for blast in blast_list:
            if contig_type == QUERY:
                search_name = blast.qname
            else:
                search_name = blast.sname
            if self.name == search_name and blast.type == blast_type:
                temp.append(blast)
        return temp

    def MULTI_NULL_POSITIONS(self):
        u""" マルチ型の間の空白の座標
        自分を含むBLASTにマルチ型のものがある場合、その全ての空白部分の座標を出力する """
        blast_lists = [self.start_link_list, self.end_link_list, self.contain_list, self.contained_list]
        positions = [] #returnするリスト
        for blast_list in blast_lists:
            blast_list.sort(key = lambda x:x.row)
            for blast in blast_list:
                if isinstance(blast,Multi_blast): #Multi_blastに該当するならば
                    position = blast.multi_null_positions #その空白部分を取得して
                    positions.extend(position) #全てリストに格納する
                    break
        return positions

    def FIND_ATCG(self,seq):
        u"""Nでない塩基配列を見つける関数
        塩基配列を最初と最後からそれぞれ探索し、初めてN以外の塩基が見つかったときそのindexを返す。
        見つからないときは-1を返す"""
        seq = str(seq).upper() #大文字に統一
        l = [seq.find('A'), seq.find('T'), seq.find('G'), seq.find('C')]
        if -1 in l: #無かったとき(結果に-1が含まれるとき)、l内の-1を全て削除する
            for i in range(l.count(-1)):
                l.remove(-1)
        try:
            index_first = min(l) #昇順で初めてATCGが見つかった塩基番号
        except:
            index_first = -1 #一つも見つからなかったときはlは空になるのでエラー。-1を返させる
        index_last = max(seq.rfind('A'), seq.rfind('T'), seq.rfind('G'), seq.rfind('C')) #降順で初めてATCGが見つかった塩基番号
        return [index_first, index_last]

    def N_POSITION_SEARCH(self):
        u"""FASTA ファイルを読み込んで、Nの位置を出力する
        結果を[(Nの開始位置, Nの終了位置),(Nの開始位置, Nの終了位置),...]のように返す
        何をしているか定かでは無いがとりあえず正しく動く"""
        fasta = self.fasta
        nnn_info_list = []
        if not fasta:
            pass
        else:
            seq = fasta.seq
            if 'NN' in seq or 'nn' in seq:
                seq = str(seq).upper()#大文字にする
                first = 0
                last = len(seq)-1
                while True:#両側から探索している
                    if seq.find('N') == -1:
                        break
                    first = first + seq.find('N')
                    last = last - len(seq) + 1 + seq.rfind('N')
                    seq = seq[seq.find('N') : seq.rfind('N') + 1]
                    nnn_info_list.append(first)
                    nnn_info_list.append(last)
                    b = self.FIND_ATCG(seq)
                    if b[0] == -1:
                        break
                    first = first + b[0]
                    last = last - len(seq) + 1 + b[1]
                    nnn_info_list.append(first - 1)
                    nnn_info_list.append(last + 1)
                    seq = seq[b[0] : b[1] + 1]
                nnn_info_list.sort()
        l = []
        for i in range(0,len(nnn_info_list),2):
            l.append((nnn_info_list[i], nnn_info_list[i + 1]))
        return l

    def n_replace(self, dev):
        u"""Nを対応するサブジェクトで置き換える
        マルチ型のデータをもとに、対応するPacBioの配列を割り出す
        クエリしか置きかえできないので注意"""
        plus_val = 0
        for n_posi in self.n_position: #n_posi = (N_start, N_end)
            frag = False
            for null_posi_set in self.multi_null_positions: #null_posi_set = (null_posi_q, null_posi_s)
                null_posi_q = null_posi_set[0] #null_posi_q = (qname, start, end) マルチ型においてqueryの空白部分
                null_posi_s = null_posi_set[1] #null_posi_s = (sname, start, end) マルチ型においてsubjectの空白部分
                if n_posi[0] - dev < null_posi_q[1] < n_posi[0] + dev and\
                n_posi[1] - dev < null_posi_q[2] < n_posi[1] + dev: # n_posiとnull_posi_qが近い値を示すとき
                    dev0 = n_posi[0] - null_posi_q[1] #差を保存して
                    dev1 = n_posi[1] - null_posi_q[2] #差を保存して
                    qname = null_posi_q[0]
                    sname = null_posi_s[0]
                    if null_posi_s[1] < null_posi_s[2]:
                        replacing_posi = [null_posi_s[1] + dev0, null_posi_s[2] + dev1] #subjectをその分ずらす
                    else:
                        replacing_posi = [null_posi_s[1] - dev0, null_posi_s[2] - dev1] #subjectをその分ずらす
                    #以上でn_posiに対応するPacBioの座標がreplacing_posiに格納された
                    frag = True
                    break
            if not frag:
                continue
            replaced_posi = n_posi
            query = search_object(qname, query_list) #qnameからクエリオブジェクトを取得する
            subject = search_object(sname, subject_list) #同様
            replaced_seq = self.fasta.seq #Miseq配列
            replacing_seq = subject.fasta.seq #PacBio配列
            s_start = replaced_posi[0] + plus_val #plus_valは置き換えで生じた塩基数のずれ
            s_end = replaced_posi[1] + plus_val 
            q_start = replacing_posi[0] 
            q_end = replacing_posi[1]
            #q => s のように上書きする
            #qをcut_seqとして切り出し、sの置き換えない部分と結合する
            if q_start < q_end:
                cut_seq = replacing_seq[q_start: q_end + 1]
                seq1 = replaced_seq[:s_start]
                seq2 = replaced_seq[s_end + 1:]
            else: #reverseの時は
                cut_seq = replacing_seq[q_end: q_start + 1]
                seq1 = replaced_seq[:s_start]
                seq2 = replaced_seq[s_end + 1:]
                cut_seq = cut_seq.reverse_complement() #相補鎖をとる必要がある

            self.fasta.seq = seq1 + cut_seq + seq2
            plus_val += len(self.fasta.seq) - len(replaced_seq) #置き換え後の長さ - 置き換え前の長さをずれとして保存しておく

    def contain_replace(self,length = 100000000):
        u"""内包しているものに置き換える"""
        plus_val = 0
        replaced_seq = self.fasta.seq
        for blast in self.contain_well_list:
            subject_posi = (blast.sstart, blast.send)
            query_posi = (blast.qstart, blast.qend)
            if isinstance(self, Query):
                subject = search_object(blast.sname, subject_list)
                replacing_seq = subject.fasta.seq
                replaced_posi = query_posi
                replacing_posi = subject_posi
            elif isinstance(self, Subject):
                query = search_object(blast.qname, query_list)
                replacing_seq = query.fasta.seq
                replaced_posi = subject_posi
                replacing_posi = query_posi
            if len(replacing_seq) < length:
                continue
            s_start = replaced_posi[0] + plus_val
            s_end = replaced_posi[1] + plus_val
            q_start = replacing_posi[0]
            q_end = replacing_posi[1]

            if s_start < s_end:
                if q_start < q_end:
                    cut_seq = replacing_seq[q_start: q_end + 1]
                    seq1 = replaced_seq[:s_start]
                    seq2 = replaced_seq[s_end + 1:]
                else:
                    cut_seq = replacing_seq[q_end: q_start + 1]
                    seq1 = replaced_seq[:s_start]
                    seq2 = replaced_seq[s_end + 1:]
                    cut_seq = cut_seq.reverse_complement()
                full_seq = seq1 + cut_seq + seq2

            else:
                if q_start < q_end:
                    cut_seq = replacing_seq[q_start: q_end + 1]
                    seq1 = replaced_seq[:s_end]
                    seq2 = replaced_seq[s_start + 1:]
                else:
                    cut_seq = replacing_seq[q_end: q_start + 1]
                    seq1 = replaced_seq[:s_end]
                    seq2 = replaced_seq[s_start + 1:]
                    cut_seq = cut_seq.reverse_complement()
                full_seq = seq1 + cut_seq + seq2
                full_seq = full_seq.reverse_complement()

            self.fasta.seq = full_seq
            plus_val += len(self.fasta.seq) - len(replaced_seq)
        self.replace_dev = len(self.fasta.seq) - self.len
        

class Query(contig):
    def __init__(self,qname,QUERY_FASTA,blast_list):
        self.name = qname
        self.type = QUERY
        self.len = self.LENGTH(blast_list)
        self.fasta = QUERY_FASTA.search_seq(self.name)
        contig.__init__(self,qname,blast_list)
        self.n_position = self.N_POSITION_SEARCH()
        self.contain_well_list = self.CONTAIN_WELL_LIST()

    def LENGTH(self,blast_list):
        for blast in blast_list:
            if blast.qname == self.name:
                i_want_this = blast
                break
        return i_want_this.qlen

    def CONTAIN_WELL_LIST(self):
        u"""被内包型に関して、コンティグ同士の重複を削除したリスト"""
        contain_list = self.contained_list[:]
        contain_well_list = []
        if contain_list:
            contain_list.sort(key = lambda blast:blast.qstart)
            contain_well_list.append(contain_list[0])
            lower_limit = contain_list[0].qend
            contain_list.remove(contain_list[0])
            for blast in contain_list:
                if blast.qstart > lower_limit - 1000:
                    contain_well_list.append(blast)
                    lower_limit =  blast.qend
                else:
                    continue
        return contain_well_list

class Subject(contig):
    def __init__(self,sname,blast_list):
        self.name = sname
        self.type = SUBJECT
        self.len = self.LENGTH(blast_list)
        self.fasta = SUBJECT_FASTA.search_seq(self.name)
        contig.__init__(self,sname,blast_list)
        self.n_position = self.N_POSITION_SEARCH()
        self.contain_well_list = self.CONTAIN_WELL_LIST()

    def LENGTH(self,blast_list):
        for blast in blast_list:
            if blast.sname == self.name:
                i_want_this = blast
                break
        return i_want_this.slen


    def CONTAIN_WELL_LIST(self):
        u"""内包型に関して、コンティグ同士の重複を削除したリスト"""
        contain_list = self.contain_list[:]
        contain_well_list = []
        if contain_list:
            contain_list.sort(key = lambda blast:min(blast.subject_using_nucl_edge_nums))
            contain_well_list.append(contain_list[0])
            lower_limit = max(contain_list[0].subject_using_nucl_edge_nums)
            for blast in contain_list:
                if min(blast.subject_using_nucl_edge_nums) > lower_limit - 1000:
                    contain_well_list.append(blast)
                    lower_limit = max(blast.subject_using_nucl_edge_nums)
                else:
                    continue
        return contain_well_list

class Fragment:
    u"""contigオブジェクト、開始位置、終了位置を持つオブジェクト
    self.seqで配列データを返す。開始位置　>　終了位置なら相補鎖を返す。"""
    def __init__(self, contig, start, end):
        self.contig = contig
        self.start = start
        self.end = end
        self.seq = self.SEQ()

    def SEQ(self):
        if self.start < self.end:
            return self.contig.fasta.seq[self.start:self.end + 1]
        else:
            return self.contig.fasta.seq[self.end: self.start + 1].reverse_complement()

    def check(self, log_text = None):
        if log_text is None:
            print "contig:",self.contig.name
            print "start:",self.start
            print "end:",self.end
            print "direction:",("+" if self.start < self.end else "-")
        else:
            log_text.write("contig:%s (%d...%d)\n"%(self.contig.name, 0, self.contig.len - 1))
            log_text.write("start:%d\n"%(self.start))
            log_text.write("end:%d\n"%(self.end))
            log_text.write("direction:%s\n\n"%("+" if self.start < self.end else "-"))


class Sqs:
    u"""SQS(サブジェクト-クエリ-サブジェクトの繋がり)のオブジェクト
    2つのサブジェクトのデータは、self.blasts, self.snamesなどにタプルとして格納されており、
    0番目が始点側、1番目が終点型側のデータ
    self.typeはサブジェクトを2つ持つときのみTrueを返す"""
    START = 0
    END = 1
    def __init__(self,query,blast_list):#queryはオブジェクトとして、サブジェクトはsnameとして扱う
        self.query = query
        self.name = self.query.name
        self.blasts = (self.DECIDE_BLAST(START_LINK),self.DECIDE_BLAST(END_LINK))
        self.snames = map(lambda x:x and x.sname or x, self.blasts)
        self.directions = map(lambda x:x and int(x.is_reverse + 1) or x,self.blasts) #MiSeqに対してBLASTが正向きなら1、逆向きなら2
        self.type = bool(self.blasts[START] and self.blasts[END]) #始点側、終点側に共にBLASTが存在すればTrue,そうでなければFalseを返します
        if self.type == True:
            for blast in blast_list:
                for a_blast in self.blasts:
                    if a_blast == blast:
                        blast.in_sqs = True #欠損ないSQSに使われているBLASTについてblast.in_sqs = True

    def DECIDE_BLAST(self,blast_type):
        u"""始点、終点に使うBLASTを確定させる関数"""
        #able_subject_list = 始点側を探索中なら始点型、終点側を探索中なら終点側のリスト
        #another_subject_list = able_subject_listの逆側のリスト
        if blast_type == START_LINK:
            able_subject_list = self.query.start_link_list[:] 
            another_subject_list = self.query.end_link_list[:]
        elif blast_type == END_LINK:
            able_subject_list = self.query.end_link_list[:]
            another_subject_list = self.query.start_link_list[:]
        elif blast_type == CONTAIN: #始点型、または終点型が無かった場合に用いる
            able_subject_list = self.query.contain_list[:]
            another_subject_list = self.query.contain_list[:]

        if len(able_subject_list) >= 1: #候補が１つ以上あるものは最も長いものを採用
            return max(able_subject_list, key = lambda x:x.slen)
        else:#候補がないときは
            if len(another_subject_list) in (1,2):#もし始点側、終点側の両方が無い等でなければ
                return self.DECIDE_BLAST(CONTAIN) #内包型を代わりに採用する
            else:
                return None

    def Is_connectable(self,sqs,my_hand,other_hand):
        if self.snames[my_hand] == sqs.snames[other_hand]:
            if (my_hand + other_hand\
            + self.directions[my_hand] + sqs.directions[other_hand])%2 == 1: #排他的論理和マジックナンバー方式
                return True
            else:
                return False
        else:
            False

    def decide_pointas(self):
        u"""pointa == 配列を繋ぎ合わせる際にどこでクエリサブジェクトを切り替えるか
        queryではreplace_devは後半 i.e. 塩基番号が大きい方のポインタに加える
        replace部分より後にあるか前にあるかがポイントとなる"""
        subjects = [search_object(sname, subject_list) for sname in self.snames] #2PacBios
        #ポインタの位置はpowerpoint資料参照
        self.q_pointas = [self.blasts[START].qend + 1, self.blasts[END].qstart  + self.query.replace_dev - 1]
        self.s_pointas = [self.blasts[START].send, self.blasts[END].sstart]
        if self.q_pointas[0] > self.q_pointas[1]: #例外処理、PacBiosが重なっているとき
            exception_dev = self.q_pointas[0] - self.q_pointas[1] #subjectがどれくらい重なっているか
            self.q_pointas[1] = self.q_pointas[0] #一方のポインタを移動するだけで対処可能
            if self.blasts[END].is_reverse: #向きによってポインタを移動する向きが異なる
                self.s_pointas[1] -= exception_dev
            else:
                self.s_pointas[1] += exception_dev
        #以下置き換え誤差処理
        if not self.blasts[START].is_reverse:
            self.s_pointas[0] += subjects[START].replace_dev
        if self.blasts[END].is_reverse:
            self.s_pointas[1] += subjects[END].replace_dev

class Sqs_chain:
    u"""SQS(サブジェクト-クエリ-サブジェクトの繋がり)を連ねたオブジェクト
    typeは環状ならcycle, 直鎖ならlineを返す
    lenは初期値から再設定されないのでseqなどの変更時は注意"""
    def __init__(self, sqs_chains, num):
        self.name = num
        self.seq = sqs_chains.seqs[num]
        self.string = sqs_chains.strings[num]
        self.type = sqs_chains.strings[num][-1]
        self.len = len(self.string) - 1

    def check(self):
        print "seq: ",
        for i in range(self.len):
            print self.string[i],
        print ""
        print "type: " + self.type

    def make_fasta_data(self, log_text):
        u"""sqs_chainからfastaデータを作成"""
        fragment_list = []
        for obj in self.seq:
            if isinstance(obj, Sqs):
                obj.decide_pointas() #pointa == blast,sqsでcontigを切り替えるポジション
        for (i,obj) in enumerate(self.seq):
            if i + 1 <= len(self.seq) - 1:
                next_obj = self.seq[i + 1]
                if i == 0:
                    if isinstance(obj,Blast):
                        query = search_object(obj.qname, query_list)
                        subject = search_object(obj.sname, subject_list)
                        next_s_pointa = next_obj.s_pointas[next_obj.chain_direction - 1]
                        if obj.type == START_LINK:
                            fragment_list.append(Fragment(query, obj.qlen - 1 + query.replace_dev, \
                            obj.qend + 1))
                            s_fragment = Fragment(subject, obj.send, next_s_pointa)
                            if not obj.is_reverse:
                                s_fragment.start += subject.replace_dev
                            fragment_list.append(s_fragment)
                        else:
                            fragment_list.append(Fragment(query, 0, obj.qstart - 1 + query.replace_dev))
                            s_fragment = Fragment(subject, obj.sstart, next_s_pointa)
                            if obj.is_reverse:
                                s_fragment.start += subject.replace_dev
                            fragment_list.append(s_fragment)
                        continue
                    else:
                        subject = search_object(obj.snames[obj.chain_direction - 1],subject_list)
                        if obj.chain_direction == PLUS:
                            end = obj.s_pointas[START]
                            temp = END
                            s_is_reverse= obj.blasts[START].is_reverse
                        else:
                            end = obj.s_pointas[END]
                            temp = START
                            s_is_reverse= bool(not obj.blasts[END].is_reverse)
                        if self.type == "line":
                            if not s_is_reverse:
                                start = 0
                            else:
                                start = subject.len - 1 + subject.replace_dev
                        else:
                            last_sqs = self.seq[-1]
                            start = last_sqs.s_pointas[temp]
                        fragment_list.append(Fragment(subject, start, end))
                query = obj.query
                subject = search_object(obj.snames[2 - obj.chain_direction],subject_list)
                now_s_pointa = obj.s_pointas[2 - obj.chain_direction]
                if isinstance(next_obj, Sqs):
                    next_s_pointa = next_obj.s_pointas[next_obj.chain_direction - 1]
                else:
                    if next_obj.type == START_LINK:
                        next_s_pointa = next_obj.send
                        if not next_obj.is_reverse:
                            next_s_pointa += subject.replace_dev
                    else:
                        next_s_pointa = next_obj.sstart
                        if next_obj.is_reverse:
                            next_s_pointa += subject.replace_dev
                if obj.chain_direction == PLUS:
                    fragment_list.append(Fragment(query, obj.q_pointas[START],obj.q_pointas[END]))
                else:
                    fragment_list.append(Fragment(query, obj.q_pointas[END],obj.q_pointas[START]))
                fragment_list.append(Fragment(subject, now_s_pointa, next_s_pointa))
            else:
                if isinstance(obj,Blast):
                    if i == 0:
                        subject = search_object(obj.sname, subject_list)
                        if obj.type == START_LINK:
                            if not obj.is_reverse:
                                fragment_list.append(Fragment(subject, 0, obj.send))
                            else:
                                fragment_list.append(Fragment(subject, obj.slen - 1, obj.send))
                        elif obj.type == END_LINK:
                            if not obj.is_reverse:
                                fragment_list.append(Fragment(subject, obj.slen - 1, obj.sstart))
                            else:
                                fragment_list.append(Fragment(subject, 0, obj.sstart))
                        else:
                            fragment_list.append(Fragment(subject, 0, obj.slen - 1))
                    query = search_object(obj.qname, query_list)
                    if obj.type == START_LINK:
                        fragment_list.append(Fragment(query, obj.qend + 1, obj.qlen - 1 + query.replace_dev))
                    else:
                        fragment_list.append(Fragment(query, obj.qstart + query.replace_dev - 1, 0))
                else:
                    query = obj.query
                    subject = search_object(obj.snames[2 - obj.chain_direction],subject_list)
                    now_s_pointa = obj.s_pointas[2 - obj.chain_direction]
                    if obj.directions[2 - obj.chain_direction] == obj.chain_direction:
                        next_end_point = subject.len - 1 + subject.replace_dev
                    else:
                        next_end_point = 1
                    if obj.chain_direction == PLUS:
                        fragment_list.append(Fragment(query, obj.q_pointas[START],obj.q_pointas[END]))
                    else:
                        fragment_list.append(Fragment(query, obj.q_pointas[END],obj.q_pointas[START]))
                    if self.type == "line":
                        fragment_list.append(Fragment(subject, now_s_pointa, next_end_point))
        full_seq = fragment_list[0].seq
        if log_text is not None:
            log_text.write(">>seq:%s\n\n"%(self.name))
        for fragment in fragment_list:
            fragment.check(log_text)
            if fragment != fragment_list[0]:
                full_seq += fragment.seq
        full_seq = SeqRecord(full_seq)
        full_seq.name = str(self.name)
        full_seq.id = str(self.name)
        full_seq.description = ''
        return full_seq


class Sqs_chains:
    u""" SQSリストからstrings(文字列),sqs_chains(sqs_chainの集まり)を生成するためのクラス
    生成時点ではsqsのみからなっているので完成させるにはself.add_blast()を実行する必要あり
    self.strings ∈ sqs_chain.string, self.seqs ∈ sqs_chain.seq """

    def __init__(self, name, sqs_list):
        self.name = name
        self.sqs_list = sqs_list
        self.strings = self.MAKE_SHIRITORIS(sqs_list)[0]
        self.seqs = self.MAKE_SHIRITORIS(sqs_list)[1]
        self.len = len(self.strings) #contig(MiSeq + PacBio)の数
        self.sqs_chains = self.MAKE_SQS_CHAIN_OBJECT()

    def MAKE_SHIRITORIS(self, sqs_list):
        u"""SQSから配列の概形を作る。
        完成させるにはself.add_blastを実行させる必要あり。"""
        def reverse_sqs_chain(sqs_chain):
            u"""sqsリストを反転させる関数"""
            for sqs in sqs_chain:
                sqs.chain_direction = 3 - sqs.chain_direction
            sqs_chain.reverse()
            return sqs_chain

        sqs_list = copy.deepcopy(sqs_list)
        shiritori_list = [] #contig名と向きが入った文字列リストをここに
        sqs_chains = [] #sqs_chain = sqsを並べたリストをここに
        while sqs_list:
            now_sqs = sqs_list.pop()
            FIRST_TYPE = START #最初に取り出したsqsの手の種類
            now_TYPE = START #今検索対象である手の種類
            sqs_chain = [now_sqs] #sqs鎖データ
            now_sqs.chain_direction = MINUS #今の鎖の向き、MINUSなのはのちにひっくり返すため
            query_subject_chain = ['+' + now_sqs.query.name] #contig名と向きが入るリスト。なくても配列には影響しない
            flag = 0
            while True:
                """now_sqs, now_TYPEは検索しているSQSと始点or終点、sqs, TYPEは候補のSQSと始点or終点"""
                for sqs in sqs_list:
                    for TYPE in [START, END]:
                        if now_sqs.Is_connectable(sqs,now_TYPE,TYPE): #接続判定
                            if now_sqs.directions[now_TYPE] == sqs.directions[TYPE]: #directions = TYPE側のPacBioの向き が等しいとき
                                sqs.chain_direction = now_sqs.chain_direction #SQSの向きはnow_sqsとsqsで等しい
                            else:
                                sqs.chain_direction = 3 - now_sqs.chain_direction #SQSの向きはnow_sqsとsqsで等しくない [1,2]なので3から引いている
                            sub_direction = [sqs.directions[TYPE],3 - sqs.directions[TYPE]][sqs.chain_direction - 1] #今ヒットしたサブジェクトの向き
                            sqs_chain.append(sqs)
                            sqs_list.remove(sqs)
                            query_subject_chain.append(direction_mark[sub_direction] + sqs.snames[TYPE]) # "+contig_24"のような形式。向き+コンティグ名
                            query_subject_chain.append(direction_mark[sqs.chain_direction] + sqs.query.name) #この2行でPacBio,Miseqの順に追加する
                            now_sqs = sqs #更新
                            now_TYPE = int(not TYPE) #始点終点更新、今ヒットした方とは逆側のPacBioを検索する
                            flag = 1 #ヒットしたことを示すフラグ
                            break
                if flag: #sqsがヒットした場合
                    flag = 0
                    continue
                #これ以降はヒットするSQSが無かった場合
                if FIRST_TYPE == START:#ひっくり返す前なら
                    if sqs_chain[0].Is_connectable(now_sqs, END, now_TYPE):#今検索対象となっているSQSがはじめのSQSに対して接続可能なら　=> 環状
                        query_subject_chain.reverse()
                        reverse_sqs_chain(sqs_chain) #最初にMINUSにしたのでひっくり返してPLUSにしておく
                        sub_direction = sqs_chain[-1].directions[END] #
                        query_subject_chain.append(direction_mark[sub_direction] + now_sqs.snames[now_TYPE])
                        query_subject_chain.append('circle') #最後に環状である旨を記載、このリストで実際に配列に影響するのはここのみ
                        shiritori_list.append(query_subject_chain) #完成したのでリストに格納
                        sqs_chains.append(sqs_chain) #同様
                        break
                    else: #そうでないならひっくり返して接続を続ける
                        sub_direction = [now_sqs.directions[now_TYPE],3 - now_sqs.directions[now_TYPE]][now_sqs.chain_direction - 1]
                        query_subject_chain.append(direction_mark[sub_direction] + now_sqs.snames[now_TYPE])
                        now_sqs = sqs_chain[0]
                        query_subject_chain.reverse() #反転
                        reverse_sqs_chain(sqs_chain) #反転、sqs_chainは処理が面倒なので関数化
                        FIRST_TYPE = END #最初の検索対象を終点側に変える
                        now_TYPE = END
                        continue
                else:#ひっくり返した後なら => 鎖状
                    sub_direction = [now_sqs.directions[now_TYPE],3 - now_sqs.directions[now_TYPE]][now_sqs.chain_direction - 1]
                    query_subject_chain.append(direction_mark[sub_direction] + now_sqs.snames[now_TYPE])
                    query_subject_chain.append('line')#最後に鎖状である旨を記載。
                    shiritori_list.append(query_subject_chain) #完成したのでリストに格納
                    sqs_chains.append(sqs_chain) #同様
                    break

        return (shiritori_list, sqs_chains) 

    def add_blast(self, blast_list):
        u"""MAKE_SHIRITORISで作ったseqデータの端に、さらにblastを繋げる関数"""
        for (i,sqs_chain) in enumerate(self.sqs_chains): #i番目のsqs_chainを取っている
            seq = sqs_chain.seq #ここでのseqはsqs_chainのSQSを並べたもの => [sqs1, sqs2, ...]
            if sqs_chain.type == 'line':
                first_blast = seq[0].blasts[seq[0].chain_direction - 1] #sqs_chainの最初に来るPacBioを含むBLAST
                last_blast = seq[-1].blasts[2 - seq[-1].chain_direction] #sqs_chainの最後に来るPacBioを含むBLAST
                more_first_blast = None #first_blastに接続するBLASTがここに
                more_last_blast = None #last_blastに接続するBLASTがここに
                for blast in blast_list:
                    try:
                        if blast.in_sqs == True: #SQSに使われているものは除外
                            continue
                    except:
                        pass
                    if blast.type > 2: #START_LINK, END_LINK以外を弾く
                        continue
                    if first_blast.Is_connectable(blast, START_LINK):
                        if more_first_blast:
                            if more_first_blast.slen > blast.slen: #slen = PacBioの全長が長いものに次々更新していく
                                continue
                        more_first_blast = blast
                    if last_blast.Is_connectable(blast, END_LINK):
                        if more_last_blast:
                            if more_last_blast.slen > blast.slen:
                                continue
                        more_last_blast = blast
                if more_first_blast:
                    sqs_chain.seq.insert(0, more_first_blast) #SQS_CHAINの手前に追加(ここが場合わけブーストの要因なので改良の余地あり)
                    sqs_chain.string.insert(0, direction_mark[3 - more_first_blast.type] + more_first_blast.qname)
                    sqs_chain.len += 1
                if more_last_blast:
                    sqs_chain.seq.append(more_last_blast) #SQS_CHAINの最後に追加(ここも場合わけブーストの要因なので改良の余地あり)
                    sqs_chain.string.insert(-1, direction_mark[more_last_blast.type] + more_last_blast.qname)
                    sqs_chain.len += 1

    def MAKE_SQS_CHAIN_OBJECT(self):
        u"""一応SQS_CHAINSから定義してしまったので、内部データとしてはSQS_CHAINオブジェクトにまとめる必要があり"""
        return [Sqs_chain(self, i) for i in range(self.len)]
    
    def USING_QUERYS(self, query_list, subject_list):
        u"""このデータ中で使用したクエリ(MiSeq)の一覧(setオブジェクト)
        SQSのMiseqだけでなく、使用したPacBioに内包されるMiSeqも含む"""
        l = []
        for sqs_chain in self.sqs_chains:
            for name in sqs_chain.string: #sqs_chain.string、思ったほど蛇足では無かったです
                if search_object(name[1:], query_list): #[1:]で0文字目を消しているのは"+""-"が含まれるため
                    query = search_object(name[1:], query_list) #query_listにあるならMiSeqなので
                    l.append(query) #そのままリストに追加
                elif search_object(name[1:], subject_list):
                    subject = search_object(name[1:], subject_list) #subject_listにあるならPacBioなので
                    contain_qname = list(set([blast.qname for blast in subject.contain_well_list])) #内包するMiSeq名を検索
                    contain_query = [search_object(qname, query_list) for qname in contain_qname] #qname => queryに変換
                    l.extend(contain_query) #まとめて追加
        l = set(l) #SETに変換することで集合演算が可能に
        return l
    
    def add_chain(self, seq_data, string_data = None):
        """sqsとblastが入ったリストを受け付ける
        string_dataがない場合はself.stringsがNoData扱いとなる（checkで参照できない）"""
        if string_data:
            self.strings.append(string_data)
        else:
            self.strings.append(["NoData", "line"])
        self.seqs.append(seq_data)
        self.len += 1
        self.sqs_chains.append(Sqs_chain(self, self.len - 1))

def search_object(obj_name, search_list):
    u"""クエリやサブジェクト、SQSをリストから検索するだけ
    blastは　row:"番号" 　contigはcontig名　sqsはクエリ名
    名前は一部でも可、完全一致を優先"""
    for i in search_list:
        if obj_name == i.name:
            return i
            break
    for i in search_list:
        if obj_name in i.name:
            return i

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

query_list = [Query(qname,QUERY_FASTA,well_blast_list) for qname in qname_set]
subject_list = [Subject(sname,well_blast_list) for sname in sname_set]

if toRemoveN:
    for query in query_list:
        query.n_replace(30)
    n_replaced_list = [x.fasta.make_seqrecord() for x in query_list]
    SeqIO.write(n_replaced_list, nnn_removed_name, "fasta")

    N_BLAST_NAME = "tes.csv"
    Do_blast(SUBJECT_FASTA_NAME, nnn_removed_name, N_BLAST_NAME)

    well_blast_list = make_blast_list(N_BLAST_NAME)

    QUERY_FASTA_DATA = SeqIO.parse(nnn_removed_name, "fasta")
    SUBJECT_FASTA_DATA = SeqIO.parse(SUBJECT_FASTA_NAME,"fasta")
    QUERY_FASTA = Simple_multi_fasta(QUERY_FASTA_DATA)
    SUBJECT_FASTA = Simple_multi_fasta(SUBJECT_FASTA_DATA)

    qname_set = set(map(lambda x:x.qname,well_blast_list))
    sname_set = set(map(lambda x:x.sname,well_blast_list))

    query_list = [Query(qname,QUERY_FASTA,well_blast_list) for qname in qname_set]
    subject_list = [Subject(sname,well_blast_list) for sname in sname_set]

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
    contig.contain_replace(5000) #5000塩基以上のもので置き換え

#置き換えているもののfasta化
#s_seq_list = [x.fasta.make_seqrecord() for x in query_list]
#SeqIO.write(s_seq_list,"okikae_260_s_5000.fasta", "fasta")

#print "replace done"

#fastaファイルとして出力
complete_data = [sqs_chain.make_fasta_data(log_text) for sqs_chain in shiritori.sqs_chains]
SeqIO.write(complete_data, OUTPUT_FASTA_NAME, "fasta")