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
	
	comp_AA = Decil(all_AA) / total_all
	comp_AT = round(float(all_AT) / total_all, 2)
	comp_AG = round(float(all_AG) / total_all, 2)
	comp_AC = round(float(all_AC) / total_all, 2)
	comp_TA = round(float(all_TA) / total_all, 2)
	comp_TT = round(float(all_TA) / total_all, 2)
	comp_TG = round(float(all_TG) / total_all, 2)
	comp_TC = round(float(all_TC) / total_all, 2)
	comp_GA = round(float(all_GA) / total_all, 2)
	comp_GT = round(float(all_GT) / total_all, 2)
	comp_GG = round(float(all_GG) / total_all, 2)
	comp_GC = round(float(all_GC) / total_all, 2)
	comp_CA = round(float(all_CA) / total_all, 2)
	comp_CT = round(float(all_CT) / total_all, 2)
	comp_CG = round(float(all_CG) / total_all, 2)
	comp_CC = round(float(all_CC) / total_all, 2)

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
		num_AT = int(round(comp_AT * total[i], 0))
		num_AG = int(round(comp_AG * total[i], 0))
		num_AC = int(round(comp_AC * total[i], 0))
		num_TA = int(round(comp_TA * total[i], 0))
		num_TT = int(round(comp_TT * total[i], 0))
		num_TG = int(round(comp_TG * total[i], 0))
		num_TC = int(round(comp_TC * total[i], 0))
		num_GA = int(round(comp_GA * total[i], 0))
		num_GT = int(round(comp_GT * total[i], 0))
		num_GG = int(round(comp_GG * total[i], 0))
		num_GC = int(round(comp_GC * total[i], 0))
		num_CA = int(round(comp_CA * total[i], 0))
		num_CT = int(round(comp_CT * total[i], 0))
		num_CG = int(round(comp_CG * total[i], 0))
		num_CC = int(round(comp_CC * total[i], 0))
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
		num_AT = round((float(AT[i]) - theory_AT[i]) * (AT[i] - theory_AT[i]) / theory_AT[i], 4)
		num_AG = round((float(AG[i]) - theory_AG[i]) * (AG[i] - theory_AG[i]) / theory_AG[i], 4)
		num_AC = round((float(AC[i]) - theory_AC[i]) * (AC[i] - theory_AC[i]) / theory_AC[i], 4)
		num_TA = round((float(TA[i]) - theory_TA[i]) * (TA[i] - theory_TA[i]) / theory_TA[i], 4)
		num_TT = round((float(TT[i]) - theory_TT[i]) * (TT[i] - theory_TT[i]) / theory_TT[i], 4)
		num_TG = round((float(TG[i]) - theory_TG[i]) * (TG[i] - theory_TG[i]) / theory_TG[i], 4)
		num_TC = round((float(TC[i]) - theory_TC[i]) * (TC[i] - theory_TC[i]) / theory_TC[i], 4)
		num_GA = round((float(GA[i]) - theory_GA[i]) * (GA[i] - theory_GA[i]) / theory_GA[i], 4)
		num_GT = round((float(GT[i]) - theory_GT[i]) * (GT[i] - theory_GT[i]) / theory_GT[i], 4)
		num_GG = round((float(GG[i]) - theory_GG[i]) * (GG[i] - theory_GG[i]) / theory_GG[i], 4)
		num_GC = round((float(GC[i]) - theory_GC[i]) * (GC[i] - theory_GC[i]) / theory_GC[i], 4)
		num_CA = round((float(CA[i]) - theory_CA[i]) * (CA[i] - theory_CA[i]) / theory_CA[i], 4)
		num_CT = round((float(CT[i]) - theory_CT[i]) * (CT[i] - theory_CT[i]) / theory_CT[i], 4)
		num_CG = round((float(CG[i]) - theory_CG[i]) * (CG[i] - theory_CG[i]) / theory_CG[i], 4)
		num_CC = round((float(CC[i]) - theory_CC[i]) * (CC[i] - theory_CC[i]) / theory_CC[i], 4)
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
