#coding: utf-8

import subprocess

START_LINK = 1
END_LINK = 2
CONTAIN = 3
CONTAINED = 4
OTHER = 0

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
        self.qend_to_edge_length = self.qlen - self.qend - 1
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
            return self.slen - self.sstart - 1
        else:
            return self.sstart

    def SEND_TO_EDGE_LENGTH(self):
        u"""条件6用の値"""
        if self.is_reverse:
            return self.send
        else:
            return self.slen - self.send - 1

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

    def Is_same(self, blast):
        if self.row != blast.row and self.qname == blast.sname and self.sname == blast.qname \
        and self.match_num == blast.match_num and self.bitscore == blast.bitscore:
            return True
        else:
            return False

    def is_equal(self):
        if self.qname == self.sname and self.qstart == self.sstart and self.qend == self.send:
            return True
        else:
            return False

    

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