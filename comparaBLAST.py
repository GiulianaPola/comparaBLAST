#!/usr/bin/python

#"/home/gpolavirus/Arthur/comparaBLAST/comparaBLAST.py"  -a Fh_com_filtro_30_bats.txt -b Fh_com_filtro_30_vir.txt -o Fh_bats_x_vir.tab
#"/home/gpolavirus/Arthur/comparaBLAST/comparaBLAST.py"  -a La_com_filtro_30_bats.txt -b La_com_filtro_30_vir.txt -o La_bats_x_vir.tab
#"/home/gpolavirus/Arthur/comparaBLAST/comparaBLAST.py"  -a Nm_com_filtro_30_bats.txt -b Nm_com_filtro_30_vir.txt -o Nm_bats_x_vir.tab
#"/home/gpolavirus/Arthur/comparaBLAST/comparaBLAST.py"  -a saida_chiroptera -b saida_virus -o chiroptera_x_virus.tab
#"/home/gpolavirus/Arthur/comparaBLAST/comparaBLAST.py" -a "/home/oocyst/itv/montagens_trinity/compara_blast/Fh_sem_filtro_Chirop.txt" -b "/home/oocyst/itv/montagens_trinity/compara_blast/Fh_sem_filtro_viruses.txt" -o Fh_chirop_x_viruses.tab


import argparse
import csv

version='1.1.3'

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

def order(filename):
    try:
      file=open(filename, 'r')
      text=file.read()
    except:
      print("File '{}' cannot be opened!".format(filename))
      quit()
    else:
      seqs = dict()
      for row in text.split('\n'):
         if not row=='':
           row=row.replace('"','').replace('  ',' ')
           columns = row.split('\t')
           qseqid = columns[0]
           bitscore = float(columns[4])
           evalue = float(columns[6])
           if len(columns)==8:
             stitle=columns[7]
           if qseqid in seqs.keys():
               item = seqs[qseqid]
               if evalue < item[1]:
                   #print(evalue,'<',item[1])
                   if len(columns)==8:
                     seqs[qseqid] = [bitscore, evalue,stitle]
                   else:
                     seqs[qseqid] = [bitscore, evalue]
               elif evalue == item[1] and bitscore > item[0]:
                   #print(evalue,'==',item[1],' and ',bitscore,'>',item[0])
                   if len(columns)==8:
                     seqs[qseqid] = [bitscore, evalue,stitle]
                   else:
                     seqs[qseqid] = [bitscore, evalue]
           if not qseqid in seqs.keys():
               if len(columns)==8:
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
            xy=[qseqid]
            if x[1] < y[1]:
#               print(x[1],'<',str(y[1]),' -> A')
                xy.append('A')
                AB[qseqid] = xy
            elif x[1] > y[1]:
#               print(x[1],'>',str(y[1]),' -> B')
                xy.append('B')
                AB[qseqid] = xy
            elif x[1] == y[1] and x[0] > y[0]:
#               print(x[1],'==',str(y[1]),' and ',x[0],' > ',y[0],' -> A')
                xy.append('A')
                AB[qseqid] = xy
            elif x[1] == y[1] and y[0] > x[0]:
#                print(x[1],'==',str(y[1]),' and ',x[0],' < ',y[0],' -> B')
                xy.append('B')
                AB[qseqid] = xy
            elif x[1] == y[1] and x[0] == y[0]:
#               print(x[1],'==',str(y[1]),' and ',x[0],' == ',y[0],' -> AB')
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
          'A_qseqid', 'A_evalue', 'B_qseqid', 'B_evalue', 'AB_qseqid', 'Selected_choice'
      ])
    elif ncols==3:
      writer.writerow([
          'A_qseqid', 'A_evalue', 'A_stitle', 'B_qseqid', 'B_evalue', 'B_stitle', 'AB_qseqid', 'Selected_choice'
      ])
    qseqids=sorted(AB.keys())+sorted(A.keys())+sorted(B.keys())
    for i in range(len(qseqids)):
        row = []
        qseqid=qseqids[i]
        if qseqid not in processed:
         if qseqid in A.keys():
           row.append(qseqid) #A_qseqid
           row.append(str(list(A[qseqid])[1])) #A_evalue
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
           row.append(str(list(B[qseqid])[1])) #B_evalue
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
           row.extend(AB[qseqid]) #AB_qseqid, Selected_choice
         else:
           row.extend(['', ''])
         writer.writerow(row)
         processed.append(qseqid)
    output.close()