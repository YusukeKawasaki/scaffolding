#coding: utf-8

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

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