#coding: utf-8

from blasttype import *

#blastの型、始点型が１、終点型が２、内包型が３、被内包型が４
START_LINK = 1
END_LINK = 2
CONTAIN = 3
CONTAINED = 4
OTHER = 0

#クエリかサブジェクトか
QUERY = 1
SUBJECT = 2

def search_object(obj_name, search_list):
    u"""クエリやサブジェクトをリストから検索するだけ
    blastは　row:"番号" 　contigはcontig名　sqsはクエリ名
    名前は一部でも可、完全一致を優先"""
    for i in search_list:
        if obj_name == i.name:
            return i
            break
    for i in search_list:
        if obj_name in i.name:
            return i

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

    def n_replace(self, query_list, subject_list, dev = 30):
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

    def contain_replace(self, query_list, subject_list, length = 100000000):
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
    def __init__(self, sname, SUBJECT_FASTA, blast_list):
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
        self.length = abs(end - start)
        self.is_reverse = bool(start - end < 0)
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

class Sequence():
    u"""BLAST群から配列を復元するオブジェクト"""
    def __init__(self, name):
        self.name = name
        self.fragments = [] #ここにFragmentオブジェクトが入り、配列となる

    def PLUNE(self, num):
        u"""配列をnumの値分だけ削る操作
        add_blast()で用いる"""
        if not self.fragments:
            return 0
        last_frag = self.fragments[-1]
        while True:
            if last_frag.length >= num: #このときlast_fragのみを削れば良い                 
                if not last_frag.is_reverse:
                    last_frag.end -= num
                else:
                    last_frag.start -= num
                self.fragments[-1] = last_frag
                return 0
            else: #このときlast_fragはまるまる削れ、次のfragmentにまで影響が出る
                self.fragments.remove(last_frag)
                num -= last_frag.length
                last_frag = self.fragments[-1]

    def add_blast(self, blast):
        u"""blast情報を元に配列を伸長する操作"""
        if blast.type == START_LINK:
            
        elif blast.type == END_LINK: