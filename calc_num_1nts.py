#! /usr/bin/env python
# encoding:utf-8

import sys
from decimal import *

'''  行列変換と総行数カウント関数 '''
def changeLine(infile_name):
	with open(infile_name, 'r') as f:
		change_seq	= {}
		counter		= 0 # 行数 (grep "^@" filename | wc と同じ)
		for line in f:
			temp = line.rstrip().lstrip("@") # 処理は"@"の配列が存在するもののみ
			counter += 1
			for i in xrange(len(temp)):
				if i == 110:
					break
				t = i + 1
				if change_seq.has_key(t):
					change_seq[t] += temp[i]
				else:
					change_seq[t] = temp[i]

	return change_seq, counter


''' 塩基毎の数カウント関数  '''
def countNucleotide(seq_dict):
	A = []
	T = []
	G = []
	C = []
	for i in sorted(seq_dict.keys()):
		seq_num_A = seq_dict[i].count("A")
		seq_num_T = seq_dict[i].count("T")
		seq_num_G = seq_dict[i].count("G")
		seq_num_C = seq_dict[i].count("C")
		A.append(seq_num_A)
		T.append(seq_num_T)
		G.append(seq_num_G)
		C.append(seq_num_C)

	return A, T, G, C	


''' 1塩基目, 2塩基目, 3塩基目 ... 毎の塩基の総数(total) カウント関数 '''
def totalSeq(A, T, G, C):
	total = []
	for i in xrange(len(A)):
		num = A[i] + T[i] + G[i] + C[i]
		total.append(num)

	return total


''' comp% 塩基の比率算出関数 '''
def compRate(A, T, G, C, total):
	all_A = all_T = all_G = all_C = total_all = 0
	for i in xrange(len(A)):
		all_A		+= A[i]
		all_T		+= T[i]
		all_G		+= G[i]
		all_C		+= C[i]
		total_all	+= total[i]
	
	comp_A = comp_T = comp_G = comp_C = 0
	
	comp_A = Decimal(all_A) / total_all
	comp_T = Decimal(all_T) / total_all
	comp_G = Decimal(all_G) / total_all
	comp_C = Decimal(all_C) / total_all
	
	return comp_A, comp_T, comp_G, comp_C


''' 理論値算出関数 '''
def calcTheorynum(comp_A, comp_T, comp_G, comp_C, total):
	theory_A = []
	theory_T = []
	theory_G = []
	theory_C = []
	for i in xrange(len(total)):
		num_A = Decimal(comp_A) * total[i]
		num_T = Decimal(comp_T) * total[i]
		num_G = Decimal(comp_G) * total[i]
		num_C = Decimal(comp_C) * total[i]
		theory_A.append(num_A)
		theory_T.append(num_T)
		theory_G.append(num_G)
		theory_C.append(num_C)

	return theory_A, theory_T, theory_G, theory_C

	
''' kai値算出関数 '''
def kaiNum(A, T, G, C, theory_A, theory_T, theory_G, theory_C):
	kai_A = []
	kai_T = []
	kai_G = []
	kai_C = []
	for i in xrange(len(A)):
		num_A = round((Decimal(A[i]) - Decimal(theory_A[i])) ** 2 / Decimal(theory_A[i]), 3)
		num_T = round((Decimal(T[i]) - Decimal(theory_T[i])) ** 2 / Decimal(theory_T[i]), 3)
		num_G = round((Decimal(G[i]) - Decimal(theory_G[i])) ** 2 / Decimal(theory_G[i]), 3)
		num_C = round((Decimal(C[i]) - Decimal(theory_C[i])) ** 2 / Decimal(theory_C[i]), 3)
		kai_A.append(num_A)
		kai_T.append(num_T)
		kai_G.append(num_G)
		kai_C.append(num_C)
	return kai_A, kai_T, kai_G, kai_C

''' kai値合計算出関数 '''
def kaiTotal(kai_A, kai_T, kai_G, kai_C):
	kai_total = []
	for i in xrange(len(kai_A)):
		num = round(kai_A[i] + kai_T[i] + kai_G[i] + kai_C[i], 4)
		kai_total.append(num)
	
	return kai_total


''' メイン '''
if __name__ == "__main__":
	infile_name			= sys.argv[1:]
	for i in xrange(len(infile_name)):
		change_seq, counter = changeLine(infile_name[i])
		A, T, G, C			= countNucleotide(change_seq)
		total				= totalSeq(A, T, G, C)
		comp_A, comp_T, comp_G, comp_C			= compRate(A, T, G, C, total)
		theory_A, theory_T, theory_G, theory_C	= calcTheorynum(comp_A, comp_T, comp_G, comp_C, total)
		kai_A, kai_T, kai_G, kai_C				= kaiNum(A, T, G, C, theory_A, theory_T, theory_G, theory_C)
		kai_total			= kaiTotal(kai_A, kai_T, kai_G, kai_C)
		print kai_A[0:5]
		print kai_T[0:5]
		print kai_G[0:5]
		print kai_C[0:5]
		print "kai値合計"
		print kai_total[0:5]
		print "実測値"
		print A[0:5]
		print T[0:5]
		print G[0:5]
		print C[0:5]
		print "理論値"
		print theory_A[0:5]
		print theory_T[0:5]
		print theory_G[0:5]
		print theory_C[0:5]
		print "Total"
		print total[0:5]
		print "Comp%"
		print str(comp_A), str(comp_T), str(comp_G), str(comp_C)
