from Bio import SeqIO
from Bio.Seq import Seq	
#from Bio.Alphabet import IUPAC

import sys


db_seq_records = SeqIO.index(sys.argv[1], "fasta")
#query_seq_records = pd.read_csv(sys.argv[2])
query_seq_records =  open(sys.argv[2], "r")

#Count= 0 # read ID from the second line
for i in query_seq_records:
	ids = str(i.split()[0])
	# Count = Count + 1
	# print ids
	# if Count >= 2:
	# ids = str(i.split()[0])
	#print(db_seq_records[ids].format("fasta"))
	print (ids," ",len(db_seq_records[ids]))#.format("fasta") # print sequence lenght

