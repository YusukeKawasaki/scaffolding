#coding: utf-8

import copy
from blasttype import *
from contig import *
from simple_fasta import *

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
        
        if len(self.query.contain_list) > 4:
            return None
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

    def decide_pointas(self, subject_list):
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

    def make_fasta_data(self, query_list, subject_list, log_text = None):
        u"""sqs_chainからfastaデータを作成"""
        fragment_list = []
        for obj in self.seq:
            if isinstance(obj, Sqs):
                obj.decide_pointas(subject_list) #pointa == blast,sqsでcontigを切り替えるポジション
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
        full_seq = fragment_list[0].seq()
        if log_text is not None:
            log_text.write(">>seq:%s\n\n"%(self.name))
        for fragment in fragment_list:
            fragment.check(log_text)
            if fragment != fragment_list[0]:
                full_seq += fragment.seq()
        full_seq = SeqRecord(full_seq)
        full_seq.name = str(self.name)
        full_seq.id = str(self.name)
        full_seq.description = ''
        return full_seq

    def make_sequence_data(self, query_list, subject_list, log_text):
        u"""fastaを返す
        sequenceクラスはfragmentの列で、add_blastを用いてblast結果を元に配列を伸長することができる"""
        sequence = Sequence(self.name)
        for obj in self.seq:
            if isinstance(obj, Blast):
                sequence.add_blast(obj, obj.type, query_list, subject_list, first_key = QUERY)
            elif isinstance(obj, Sqs):
                if self.type == "circle" and obj == self.seq[-1]:
                    if obj.chain_direction == PLUS:
                        sequence.add_blast(obj.blasts[0], START_LINK, query_list, subject_list, first_key = SUBJECT)
                        sequence.to_circular(obj.blasts[1], END_LINK)
                    else:
                        sequence.add_blast(obj.blasts[1], END_LINK, query_list, subject_list, first_key = SUBJECT)
                        sequence.to_circular(obj.blasts[0], START_LINK)
                else:
                    if obj.chain_direction == PLUS:
                        sequence.add_blast(obj.blasts[0], START_LINK, query_list, subject_list, first_key = SUBJECT)
                        sequence.add_blast(obj.blasts[1], END_LINK, query_list, subject_list)
                    else:
                        sequence.add_blast(obj.blasts[1], END_LINK, query_list, subject_list, first_key = SUBJECT)
                        sequence.add_blast(obj.blasts[0], START_LINK, query_list, subject_list)
        return sequence.output_seqrecord(log_text)

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
