#!/usr/bin/python

#python comparaBLAST.py -a Fh_com_filtro_30_bats.txt -b Fh_com_filtro_30_vir.txt -o Fh_comparado.tab
#python comparaBLAST.py -a La_com_filtro_30_bats.txt -b La_com_filtro_30_vir.txt -o La_comparado.tab
#python comparaBLAST.py -a Nm_com_filtro_30_bats.txt -b Nm_com_filtro_30_vir.txt -o Nm_comparado.tab

import argparse
import csv

ajuda = 'ComparaBLAST - comparison of sequence BLAST results against two databases\n'
ajuda=ajuda+'(c) 2021. Arthur Gruber, Giuliana Pola & Liliane S. Oliveira\n'
ajuda = ajuda + 'Usage: comparaBLAST.py -a <tabular file against database A> -b <tabular file against database B> -o <output file> \n'
ajuda = ajuda + '\nMandatory parameters:\n'
ajuda = ajuda + '-a\tTabular file against database A\n'
ajuda = ajuda + '-b\tTabular file against database B\n'
ajuda = ajuda + '\nOptional parameters:\n'
ajuda = ajuda + '-o\tOutput file (standard:"comparado.tab")\n'

parser = argparse.ArgumentParser(add_help=False)
parser.add_argument('-a')
parser.add_argument('-b')
parser.add_argument('-o')
parser.add_argument('-h', '--help', action='store_true')
args = parser.parse_args()


def ordenar(arquivo):
    seqs = dict()
    for linha in arquivo.readlines():
        colunas = linha.split('\t')
        qseqid = colunas[0]
        bitscore = float(colunas[4])
        evalue = float(colunas[6])
        if qseqid in seqs.keys():
            item = seqs[qseqid]
            if evalue < item[1]:
              #print(evalue,'<',item[1])
              seqs[qseqid] = [bitscore,evalue]
            elif evalue == item[1] and bitscore > item[0]:
                #print(evalue,'==',item[1],' and ',bitscore,'>',item[0])
                seqs[qseqid] = [bitscore,evalue]
        if not qseqid in seqs.keys():
            seqs[qseqid] = [bitscore,evalue]
    return seqs


def comparar(A, B):
    AB = dict()
    seqs = list(A.keys()) + list(B.keys())
    for qseqid in seqs:
        if qseqid in A.keys() and qseqid in B.keys():
            x = A[qseqid]
            y = B[qseqid]
            if x[1] < y[1]:
                #print(x[1],'<',y[1],' -> A')
                AB[qseqid] = [x[1], y[1], 'A']
            elif x[1] > y[1]:
                #print(x[1],'>',y[1],' -> B')
                AB[qseqid] = [x[1], y[1],'B']
            elif x[1] == y[1] and x[0] > y[0]:
                #print(x[1],'==',y[1],' and ',x[0],' > ',y[0],' -> A')
                AB[qseqid] = [x[1], y[1],'A']
            elif x[1] == y[1] and y[0] > x[0]:
                #print(x[1],'==',y[1],' and ',x[0],' < ',y[0],' -> B')
                AB[qseqid] = [x[1], y[1],'B']
            elif x[1] == y[1] and x[0] == y[0]:
                #print(x[1],'==',y[1],' and ',x[0],' == ',y[0],' -> AB')
                AB[qseqid] = [x[1], y[1],'AB']
    return AB

if args._get_args()==None:
  print(ajuda)
elif args.help == True:
    print(ajuda)
elif args.a == None or args.b == None:
    print('ERROR: Missing input file\n')
    print(ajuda)
else:
    A = ordenar(open(args.a, 'r'))
    B = ordenar(open(args.b, 'r'))
    AB = comparar(A, B)
    if args.o == None:
      saida = open('comparado.tab', mode='w')
    else:
      saida = open(args.o, mode='w')
      escritor = csv.writer(saida,
                          delimiter='\t',
                          quotechar='"',
                          quoting=csv.QUOTE_MINIMAL)
      escritor.writerow([
        'A_qseqid', 'A_evalue', 'B_qseqid', 'B_evalue', 'AB_qseqid',
        'A_evalue', 'B_evalue', 'Escolha'
    ])
      for i in range(max(len(A), len(B), len(AB))):
        linha = []
        if i < len(A):
            linha.append(list(A.keys())[i])
            linha.append(list(A[linha[0]])[0])
        else:
            linha.extend(['', ''])
        if i < len(B):
            linha.append(list(B.keys())[i])
            linha.append(list(B[linha[2]])[0])
        else:
            linha.extend(['', ''])
        if i < len(AB):
            linha.append(list(AB.keys())[i])
            linha.extend(AB[linha[4]])
        else:
            linha.extend(['', '', '', ''])
        escritor.writerow(linha)
    saida.close()