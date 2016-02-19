#! /usr/bin/env python
# encoding:utf-8

import sys
from decimal import Decimal

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
					change_seq[t] += temp[i:i+2]
				else:
					change_seq[t] = temp[i:i+2]

	return change_seq, counter


''' 塩基毎の数カウント関数  '''
def countNucleotide(seq_dict):
	AA = []
	AT = []
	AG = []
	AC = []
	TA = []
	TT = []
	TG = []
	TC = []
	GA = []
	GT = []
	GG = []
	GC = []
	CA = []
	CT = []
	CG = []
	CC = []
	for i in sorted(seq_dict.keys()):
		seq_num_AA = seq_dict[i].count("AA")
		seq_num_AT = seq_dict[i].count("AT")
		seq_num_AG = seq_dict[i].count("AG")
		seq_num_AC = seq_dict[i].count("AC")
		seq_num_TA = seq_dict[i].count("TA")
		seq_num_TT = seq_dict[i].count("TT")
		seq_num_TG = seq_dict[i].count("TG")
		seq_num_TC = seq_dict[i].count("TC")
		seq_num_GA = seq_dict[i].count("GA")
		seq_num_GT = seq_dict[i].count("GT")
		seq_num_GG = seq_dict[i].count("GG")
		seq_num_GC = seq_dict[i].count("GC")
		seq_num_CA = seq_dict[i].count("CA")
		seq_num_CT = seq_dict[i].count("CT")
		seq_num_CG = seq_dict[i].count("CG")
		seq_num_CC = seq_dict[i].count("CC")
		
		AA.append(seq_num_AA)
		AT.append(seq_num_AT)
		AG.append(seq_num_AG)
		AC.append(seq_num_AC)

		TA.append(seq_num_TA)
		TT.append(seq_num_TT)
		TG.append(seq_num_TG)
		TC.append(seq_num_TC)

		GA.append(seq_num_GA)
		GT.append(seq_num_GT)
		GG.append(seq_num_GG)
		GC.append(seq_num_GC)
		CA.append(seq_num_CA)
		CT.append(seq_num_CT)
		CG.append(seq_num_CG)
		CC.append(seq_num_CC)

	return AA, AT, AG, AC, TA, TT, TG, TC, GA, GT, GG, GC, CA, CT, CG, CC


''' 1塩基目, 2塩基目, 3塩基目 ... 毎の塩基の総数(total) カウント関数 '''
def totalSeq(AA, AT, AG, AC, TA, TT, TG, TC, GA, GT, GG, GC, CA, CT, CG, CC):
	total = []
	for i in xrange(len(AA)):
		num = AA[i] + AT[i] + AG[i] + AC[i] + TA[i] + TT[i] + TG[i] + TC[i] + GA[i] + GT[i] + GG[i] + GC[i] + CA[i] + CT[i] + CG[i] + CC[i]
		total.append(num)

	return total


''' comp% 塩基の比率算出関数 '''
def compRate(AA, AT, AG, AC, TA, TT, TG, TC, GA, GT, GG, GC, CA, CT, CG, CC, total):
	all_AA = all_AT = all_AG = all_AC = all_TA = all_TT = all_TG = all_TC = all_GA = all_GT = all_GG = all_GC = all_CA = all_CT = all_CG = all_CC = total_all = 0
	for i in xrange(len(AA)):
		all_AA		+= AA[i]
		all_AT		+= AT[i]
		all_AG		+= AG[i]
		all_AC		+= AC[i]
		all_TA		+= TA[i]
		all_TT		+= TT[i]
		all_TG		+= TG[i]
		all_TC		+= TC[i]
		all_GA		+= GA[i]
		all_GT		+= GT[i]
		all_GG		+= GG[i]
		all_GC		+= GC[i]
		all_CA		+= CA[i]
		all_CT		+= CT[i]
		all_CG		+= CG[i]
		all_CC		+= CC[i]
		total_all	+= total[i]
	
	comp_AA = comp_AT = comp_AG = comp_AC = comp_TA = comp_TT = comp_TG = comp_TC = comp_GA = comp_GT = comp_GG = comp_GC = comp_CA = comp_CT = comp_CG = comp_CC =  0
	
	comp_AA = Decimal(all_AA) / total_all
	comp_AT = Decimal(all_AT) / total_all
	comp_AG = Decimal(all_AG) / total_all
	comp_AC = Decimal(all_AC) / total_all
	comp_TA = Decimal(all_TA) / total_all
	comp_TT = Decimal(all_TA) / total_all
	comp_TG = Decimal(all_TG) / total_all
	comp_TC = Decimal(all_TC) / total_all
	comp_GA = Decimal(all_GA) / total_all
	comp_GT = Decimal(all_GT) / total_all
	comp_GG = Decimal(all_GG) / total_all
	comp_GC = Decimal(all_GC) / total_all
	comp_CA = Decimal(all_CA) / total_all
	comp_CT = Decimal(all_CT) / total_all
	comp_CG = Decimal(all_CG) / total_all
	comp_CC = Decimal(all_CC) / total_all

	return comp_AA, comp_AT, comp_AG, comp_AC, comp_TA, comp_TT, comp_TG, comp_TC, comp_GA, comp_GT, comp_GG, comp_GC, comp_CA, comp_CT, comp_CG, comp_CC


''' 理論値算出関数 '''
def calcTheorynum(comp_AA, comp_AT, comp_AG, comp_AC, comp_TA, comp_TT, comp_TG, comp_TC, comp_GA, comp_GT, comp_GG, comp_GC, comp_CA, comp_CT, comp_CG, comop_CC, total):
	theory_AA = []
	theory_AT = []
	theory_AG = []
	theory_AC = []
	theory_TA = []
	theory_TT = []
	theory_TG = []
	theory_TC = []
	theory_GA = []
	theory_GT = []
	theory_GG = []
	theory_GC = []
	theory_CA = []
	theory_CT = []
	theory_CG = []
	theory_CC = []
	for i in xrange(len(total)):
		num_AA = Decimal(comp_AA) * total[i]
		num_AT = Decimal(comp_AT) * total[i]
		num_AG = Decimal(comp_AG) * total[i]
		num_AC = Decimal(comp_AC) * total[i]
		num_TA = Decimal(comp_TA) * total[i]
		num_TT = Decimal(comp_TT) * total[i]
		num_TG = Decimal(comp_TG) * total[i]
		num_TC = Decimal(comp_TC) * total[i]
		num_GA = Decimal(comp_GA) * total[i]
		num_GT = Decimal(comp_GT) * total[i]
		num_GG = Decimal(comp_GG) * total[i]
		num_GC = Decimal(comp_GC) * total[i]
		num_CA = Decimal(comp_CA) * total[i]
		num_CT = Decimal(comp_CT) * total[i]
		num_CG = Decimal(comp_CG) * total[i]
		num_CC = Decimal(comp_CC) * total[i]
		theory_AA.append(num_AA)
		theory_AT.append(num_AT)
		theory_AG.append(num_AG)
		theory_AC.append(num_AC)
		theory_TA.append(num_TA)
		theory_TT.append(num_TT)
		theory_TG.append(num_TG)
		theory_TC.append(num_TC)
		theory_GA.append(num_GA)
		theory_GT.append(num_GT)
		theory_GG.append(num_GG)
		theory_GC.append(num_GC)
		theory_CA.append(num_CA)
		theory_CT.append(num_CT)
		theory_CG.append(num_CG)
		theory_CC.append(num_CC)

	return theory_AA, theory_AT, theory_AG, theory_AC, theory_TA, theory_TT, theory_TG, theory_TC, theory_GA, theory_GT, theory_GG, theory_GC, theory_CA, theory_CT, theory_CG, theory_CC

	
''' kai値算出関数 '''
def kaiNum(AA, AT, AG, AC, TA, TT, TG, TC, GA, GT, GG, GC, CA, CT, CG, CC, theory_AA, theory_AT, theory_AG, theory_AC, theory_TA, theory_TT, theory_TG, theory_TC, theory_GA, theory_GT, theory_GG, theory_GC, theory_CA, theory_CT, theory_CG, theory_CC):
	kai_AA = []
	kai_AT = []
	kai_AG = []
	kai_AC = []
	kai_TA = []
	kai_TT = []
	kai_TG = []
	kai_TC = []
	kai_GA = []
	kai_GT = []
	kai_GG = []
	kai_GC = []
	kai_CA = []
	kai_CT = []
	kai_CG = []
	kai_CC = []
	for i in xrange(len(AA)):
		num_AA = round((Decimal(AA[i]) - theory_AA[i]) ** 2 / theory_AA[i], 4)
		num_AT = round((Decimal(AT[i]) - theory_AT[i]) ** 2 / theory_AT[i], 4)
		num_AG = round((Decimal(AG[i]) - theory_AG[i]) ** 2 / theory_AG[i], 4)
		num_AC = round((Decimal(AC[i]) - theory_AC[i]) ** 2 / theory_AC[i], 4)
		num_TA = round((Decimal(TA[i]) - theory_TA[i]) ** 2 / theory_TA[i], 4)
		num_TT = round((Decimal(TT[i]) - theory_TT[i]) ** 2 / theory_TT[i], 4)
		num_TG = round((Decimal(TG[i]) - theory_TG[i]) ** 2 / theory_TG[i], 4)
		num_TC = round((Decimal(TC[i]) - theory_TC[i]) ** 2 / theory_TC[i], 4)
		num_GA = round((Decimal(GA[i]) - theory_GA[i]) ** 2 / theory_GA[i], 4)
		num_GT = round((Decimal(GT[i]) - theory_GT[i]) ** 2 / theory_GT[i], 4)
		num_GG = round((Decimal(GG[i]) - theory_GG[i]) ** 2 / theory_GG[i], 4)
		num_GC = round((Decimal(GC[i]) - theory_GC[i]) ** 2 / theory_GC[i], 4)
		num_CA = round((Decimal(CA[i]) - theory_CA[i]) ** 2 / theory_CA[i], 4)
		num_CT = round((Decimal(CT[i]) - theory_CT[i]) ** 2 / theory_CT[i], 4)
		num_CG = round((Decimal(CG[i]) - theory_CG[i]) ** 2 / theory_CG[i], 4)
		num_CC = round((Decimal(CC[i]) - theory_CC[i]) ** 2 / theory_CC[i], 4)
		kai_AA.append(num_AA)
		kai_AT.append(num_AT)
		kai_AG.append(num_AG)
		kai_AC.append(num_AC)
		kai_TA.append(num_TA)
		kai_TT.append(num_TT)
		kai_TG.append(num_TG)
		kai_TC.append(num_TC)
		kai_GA.append(num_GA)
		kai_GT.append(num_GT)
		kai_GG.append(num_GG)
		kai_GC.append(num_GC)
		kai_CA.append(num_CA)
		kai_CT.append(num_CT)
		kai_CG.append(num_CG)
		kai_CC.append(num_CC)
	return kai_AA, kai_AT, kai_AG, kai_AC, kai_TA, kai_TT, kai_TG, kai_TC, kai_GA, kai_GT, kai_GG, kai_GC, kai_CA, kai_CT, kai_CG, kai_CC

''' kai値合計算出関数 '''
def kaiTotal(kai_AA, kai_AT, kai_AG, kai_AC, kai_TA, kai_TT, kai_TG, kai_TC, kai_GA, kai_GT, kai_GG, kai_GC, kai_CA, kai_CT, kai_CG, kai_CC):
	kai_total = []
	for i in xrange(len(kai_AA)):
		num = round(kai_AA[i] + kai_AT[i] + kai_AG[i] + kai_AC[i] + kai_TA[i] + kai_TT[i] + kai_TG[i] + kai_TC[i] + kai_GA[i] + kai_GT[i] + kai_GG[i] + kai_GC[i] + kai_CA[i] + kai_CT[i] + kai_CG[i] + kai_CC[i], 4)
		kai_total.append(num)
	
	return kai_total


''' メイン '''
if __name__ == "__main__":
	infile_name			= sys.argv[1:]
	for i in xrange(len(infile_name)):
		change_seq, counter = changeLine(infile_name[i])
		AA, AT, AG, AC, TA, TT, TG, TC, GA, GT, GG, GC, CA, CT, CG, CC			= countNucleotide(change_seq)
		total				= totalSeq(AA, AT, AG, AC, TA, TT, TG, TC, GA, GT, GG, GC, CA, CT, CG, CC)
		comp_AA, comp_AT, comp_AG, comp_AC, comp_TA, comp_TT, comp_TG, comp_TC, comp_GA, comp_GT, comp_GG, comp_GC, comp_CA, comp_CT, comp_CG, comp_CC			= compRate(AA, AT, AG, AC, TA, TT, TG, TC, GA, GT, GG, GC, CA, CT, CG, CC, total)
		theory_AA, theory_AT, theory_AG, theory_AC, theory_TA, theory_TT, theory_TG, theory_TC, theory_GA, theory_GT, theory_GG, theory_GC, theory_CA, theory_CT, theory_CG, theory_CC	= calcTheorynum(comp_AA, comp_AT, comp_AG, comp_AC, comp_TA, comp_TT, comp_TG, comp_TC, comp_GA, comp_GT, comp_GG, comp_GC, comp_CA, comp_CT, comp_CG, comp_CC, total)
		kai_AA, kai_AT, kai_AG, kai_AC, kai_TA, kai_TT, kai_TG, kai_TC, kai_GA, kai_GT, kai_GG, kai_GC, kai_CA, kai_CT, kai_CG, kai_CC				= kaiNum(AA, AT, AG, AC, TA, TT, TG, TC, GA, GT, GG, GC, CA, CT, CG, CC, theory_AA, theory_AT, theory_AG, theory_AC, theory_TA, theory_TT, theory_TG, theory_TC, theory_GA, theory_GT, theory_GG, theory_GC, theory_CA, theory_CT, theory_CG, theory_CC)
		kai_total			= kaiTotal(kai_AA, kai_AT, kai_AG, kai_AC, kai_TA, kai_TT, kai_TG, kai_TC, kai_GA, kai_GT, kai_GG, kai_GC, kai_CA, kai_CT, kai_CG, kai_CC)
		print kai_total
