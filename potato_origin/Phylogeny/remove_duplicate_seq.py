'''
Usage: $ python remove_duplicate_seq.py {inputfile} >{outputfile}
'''
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import sys
from Bio import SeqIO
seq_dict = {}
list = ["Nicotianalongiflora"]
for seq_record in SeqIO.parse(sys.argv[1], "fasta"):
    if seq_record.id in list:
        continue
    if seq_record.id not in seq_dict:
        seq_dict[seq_record.id] = seq_record.seq
    else:
        seq_dict.pop(seq_record.id)
for i in seq_dict:
    print(">"+i)
    print(seq_dict[i])