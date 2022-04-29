#!/usr/bin/python

#"/home/gpolavirus/Arthur/comparaBLAST/comparaBLAST.py"  -a Fh_com_filtro_30_bats.txt -b Fh_com_filtro_30_vir.txt -o Fh_comparado.tab
#"/home/gpolavirus/Arthur/comparaBLAST/comparaBLAST.py"  -a La_com_filtro_30_bats.txt -b La_com_filtro_30_vir.txt -o La_comparado.tab
#"/home/gpolavirus/Arthur/comparaBLAST/comparaBLAST.py"  -a Nm_com_filtro_30_bats.txt -b Nm_com_filtro_30_vir.txt -o Nm_comparado.tab
#"/home/gpolavirus/Arthur/comparaBLAST/comparaBLAST.py"  -a output_chiroptera -b output_virus -o output_comparado.tab
#"/home/gpolavirus/Arthur/comparaBLAST/comparaBLAST.v.1.1.1.py" -a "/home/oocyst/itv/montagens_trinity/compara_blast/Fh_sem_filtro_Chirop.txt" -b "/home/oocyst/itv/montagens_trinity/compara_blast/Fh_sem_filtro_viruses.txt" -o teste.tab


import argparse
import csv

version='1.1.1'

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
      text=file.read()
    except:
      print("File '{}' cannot be opened!".getfilename(filename))
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
            if ncols==2:
              xy=[str(x[1]),y[1]]
            elif ncols==3:
              if len(x)==3 and len(y)==3:
                xy=[str(x[1]),x[2],str(y[1]),y[2]]
              elif len(x)==2 and len(y)==3:
                xy=[str(x[1]),'',str(y[1]),y[2]]
              elif len(x)==3 and len(y)==2:
                xy=[str(x[1]),x[2],str(y[1]),'']
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
          'A_qseqid', 'A_evalue', 'B_qseqid', 'B_evalue', 'AB_qseqid',
          'A_evalue', 'B_evalue', 'Escolha'
      ])
    elif ncols==3:
      writer.writerow([
          'A_qseqid', 'A_evalue', 'A_stitle', 'B_qseqid', 'B_evalue', 'B_stitle', 'AB_qseqid',
          'A_evalue',  'A_stitle', 'B_evalue', 'A_stitle', 'Selected_choice'
      ])
    for i in range(max(len(A),len(B),len(AB))):
        row = []
        if i<len(A):
          A_qseqid=sorted(A.keys())[i]
          row.append(A_qseqid) #A_qseqid
          row.append(str(list(A[A_qseqid])[1])) #A_evalue
          if ncols==3:
            if len(list(A[A_qseqid]))==3:
              row.append(list(A[A_qseqid])[2]) #A_stitle
            else:
              row.append('')
        else:
          if ncols==3:
            row.extend(['', '',''])
          elif ncols==2:
            row.extend(['', ''])
        if i<len(B):
          B_qseqid=sorted(B.keys())[i]
          row.append(B_qseqid) #B_qseqid
          row.append(str(list(B[B_qseqid])[1])) #B_evalue
          if ncols==3:
            if len(list(B[B_qseqid]))==3:
              row.append(list(B[B_qseqid])[2]) #B_stitle
            else:
              row.append('')
        else:
          if ncols==3:
            row.extend(['', '',''])
          elif ncols==2:
            row.extend(['', ''])
        if i<len(AB):
          AB_qseqid=sorted(AB.keys())[i]
          row.append(AB_qseqid) #AB_qseqid
          row.extend(AB[AB_qseqid]) #A_evalue, (A_stitle), B_evalue, (B_stitle), Escolha
        else:
          if ncols==3:
            row.extend(['', '','','','',''])
          elif ncols==2:
            row.extend(['', '','',''])
        writer.writerow(row)
    output.close()