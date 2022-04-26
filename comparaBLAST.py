#!/usr/bin/python

import argparse
import csv

version='1.1.0'

menu = 'ComparaBLAST v{} - comparison of sequence BLAST results against two databases\n'.format(version)
menu = menu + '(c) 2021. Arthur Gruber, Giuliana Pola & Liliane S. Oliveira\n'
menu = menu + 'Usage: comparaBLAST.py -a <tabular file against database A> -b <tabular file against database B> -o <output file> \n'
menu = menu + '\nMandatory parameters:\n'
menu = menu + '-a\tTabular file against database A\n'
menu = menu + '-b\tTabular file against database B\n'
menu = menu + '\nOptional parameters:\n'
menu = menu + '-o\tOutput file (default:"comparado.tab")\n'

parser = argparse.ArgumentParser(add_help=False)
parser.add_argument('-a')
parser.add_argument('-b')
parser.add_argument('-o')
parser.add_argument('-h', '--help', action='store_true')
args = parser.parse_args()
processed=[]

def getfilename(path):
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)

def order(filename):
    try:
      file=open(filename, 'r')
    except:
      print("File '{}' cannot be opened!".getfilename(filename))
      quit()
    else:
      seqs = dict()
      for row in file.readlines():
          colunas = row.split('\t')
          qseqid = colunas[0]
          bitscore = float(colunas[4])
          evalue = float(colunas[6])
          if len(colunas)==8:
            stitle=colunas[7]
          if qseqid in seqs.keys():
              item = seqs[qseqid]
              if evalue < item[1]:
                  #print(evalue,'<',item[1])
                  if len(colunas)==8:
                    seqs[qseqid] = [bitscore, evalue,stitle]
                  else:
                    seqs[qseqid] = [bitscore, evalue]
              elif evalue == item[1] and bitscore > item[0]:
                  #print(evalue,'==',item[1],' and ',bitscore,'>',item[0])
                  if len(colunas)==8:
                    seqs[qseqid] = [bitscore, evalue,stitle]
                  else:
                    seqs[qseqid] = [bitscore, evalue]
          if not qseqid in seqs.keys():
              if len(colunas)==8:
                seqs[qseqid] = [bitscore, evalue,stitle]
              else:
                seqs[qseqid] = [bitscore, evalue]
      return seqs


def compare(A, B):
    AB = dict()
    seqs = list(A.keys()) + list(B.keys())
    for qseqid in seqs:
        if qseqid in A.keys() and qseqid in B.keys():
            x = A[qseqid]
            y = B[qseqid]
            if ncols==2:
              xy=[x[1],y[1]]
            elif ncols==3:
              if len(x)==3 and len(y)==3:
                xy=[x[1],x[2],y[1],y[2]]
              elif len(x)==2 and len(y)==3:
                xy=[x[1],'',y[1],y[2]]
              elif len(x)==3 and len(y)==2:
                xy=[x[1],x[2],y[1],'']
            if x[1] < y[1]:
#               print(x[1],'<',y[1],' -> A')
                xy.append('A')
                AB[qseqid] = xy
            elif x[1] > y[1]:
#               print(x[1],'>',y[1],' -> B')
                xy.append('B')
                AB[qseqid] = xy
            elif x[1] == y[1] and x[0] > y[0]:
#               print(x[1],'==',y[1],' and ',x[0],' > ',y[0],' -> A')
                xy.append('A')
                AB[qseqid] = xy
            elif x[1] == y[1] and y[0] > x[0]:
#                print(x[1],'==',y[1],' and ',x[0],' < ',y[0],' -> B')
                xy.append('B')
                AB[qseqid] = xy
            elif x[1] == y[1] and x[0] == y[0]:
#               print(x[1],'==',y[1],' and ',x[0],' == ',y[0],' -> AB')
                xy.append('AB')
                AB[qseqid] = xy
    return AB


if args.a == None and args.b == None and args.help == False and args.o == None:
    print(menu)
elif args.help == True:
    print(menu)
elif args.a == None or args.b == None:
    print('ERROR: Missing input file\n')
    print(menu)
else:
    A = order(args.a)
    B = order(args.b)
    ncols=max(len(list(A.values())[0]),len(list(B.values())[0]))
    AB = compare(A, B)
    if args.o == None:
      output = open('compared.tab', mode='w')
    else:
      output = open(args.o, mode='w')
    writer = csv.writer(output,
                          delimiter='\t')
    if ncols==2:
      writer.writerow([
          'A_qseqid', 'A_evalue', 'B_qseqid', 'B_evalue', 'AB_qseqid',
          'A_evalue', 'B_evalue', 'Escolha'
      ])
    elif ncols==3:
      writer.writerow([
          'A_qseqid', 'A_evalue', 'A_stitle', 'B_qseqid', 'B_evalue', 'B_stitle', 'AB_qseqid',
          'A_evalue',  'A_stitle', 'B_evalue', 'A_stitle', 'Escolha'
      ])
    i=0
    qseqids=sorted(list(AB.keys()))
    qseqids.extend(sorted(list(A.keys())))
    qseqids.extend(sorted(list(B.keys())))
    for qseqid in qseqids:
        row = []
        if qseqid not in processed:
         if qseqid in A.keys():
           row.append(qseqid) #A_qseqid
           row.append(list(A[qseqid])[1]) #A_evalue
           if ncols==3:
             if len(list(A[qseqid]))==3:
               row.append(list(A[qseqid])[2]) #A_stitle
             else:
               row.append('')
         else:
           if ncols==3:
             row.extend(['', '',''])
           elif ncols==2:
             row.extend(['', ''])
         if qseqid in B.keys():
           row.append(qseqid) #B_qseqid
           row.append(list(B[qseqid])[1]) #B_evalue
           if ncols==3:
             if len(list(B[qseqid]))==3:
               row.append(list(B[qseqid])[2]) #B_stitle
             else:
               row.append('')
         else:
           if ncols==3:
             row.extend(['', '',''])
           elif ncols==2:
             row.extend(['', ''])
         if qseqid in AB.keys():
           row.append(qseqid) #AB_qseqid
           row.extend(AB[qseqid]) #A_evalue, (A_stitle), B_evalue, (B_stitle), Escolha
         else:
           if ncols==3:
             row.extend(['', '','','','',''])
           elif ncols==2:
             row.extend(['', '','',''])
         writer.writerow(row)
         processed.append(qseqid)
    output.close()