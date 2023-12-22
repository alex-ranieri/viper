import sys
from Bio import SeqIO


# sys.argv[1] -> '.Pilon.fasta' file.
# sys.argv[2] -> 'samtools depth -a ...' *.tsv file (rodar "samtools depth -a" before)



fasta_file = sys.argv[1]
depth_file = sys.argv[2]



#header = ''
new_sequence = ''

#sequence_original = ''

Lista_Depths =[]
Lista_NewPositions =[]



# Read *.FASTA
#f = open(fasta_file, 'r')
#for line in f:
#	line = line.rstrip()
#	if(line.startswith(">")):
#		header = line
#	else:
#		if(len(line) > 2):
#			sequence_original = sequence_original + line
#f.close()

record = SeqIO.read(fasta_file, "fasta")
header = ">" + record.id
sequence_original = str(record.seq)


# Read *.TSV
f = open(depth_file, 'r')

basePos = 0

for line in f:

	line = line.strip()
	temp = line.split('\t')

	depth_at_pos = int(temp[2])

	if(depth_at_pos >= 5):
		Lista_NewPositions.append(sequence_original[basePos])
	else:
		Lista_NewPositions.append('n')

	basePos += 1

f.close()



new_sequence = ''.join(Lista_NewPositions)
print(header + "\n" + new_sequence + "\n")
