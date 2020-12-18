#!/usr/bin/env python3

# adds TER cards to PDB files to be used by tleap
# TER cards are added after WAT, Cl-, and Na+ molecules, it also corrects printing of residues below 100,000 
# residues
# jose flores-canales 05/04/2013

import sys

USAGE="""
addTER.py input.pdb
"""

def main():
	if len(sys.argv)<1:
		print(USAGE)
		sys.exit()
	name_file = sys.argv[1]
	name_in = name_file.split('.pdb')[0]
	name_out=name_in+"_ter.pdb"
	infile = open(name_file,'r')
	pdbfile = open(name_out,'w')
	line =infile.readline()
	nmol=0
	while (line):
#ATOM			 1	C42 PPG X		1			 92.162 -85.409		2.210  0.00  0.00
#ATOM			 1	C42 PPG X		1			 92.162 -85.409		2.210  0.00  0.00						 
		nums=line.split()
		if len(nums)<10:
			#if nums[0] == 'END':
			#  newline="TER\n"
			#  pdbfile.write(newline)
			pdbfile.write(line)
			line=infile.readline()
			continue
		if (nums[2]=='N' or nums[2]=='N4' or (nums[2]=='H1' and nums[3]=='ACE') or  (nums[2]=='C1' and nums[3]=='CHL') or (nums[2]=='C116' and nums[3]=='PA') or (nums[2]=='C11' and nums[3]=='PC') or (nums[2]=='C12' and nums[3]=='OL') or (nums[3]=='WAT' and nums[2]=='O') or nums[3]=='Na+' or nums[3]=='Cl-'):
			nmol+=1
		newline="TER\n"
		if (nmol!=1 and ((nums[2]=='C1' and nums[3]=='CHL') or (nums[2]=='C116' and nums[3]=='PA') or (nums[2]=='C11' and nums[3]=='PC') or (nums[2]=='C12' and nums[3]=='OL') or (nums[3]=='WAT' and nums[2]=='O') or nums[3]=='Na+' or nums[3]=='Cl-')):
			pdbfile.write(newline)
			#23 - 26
		if (nmol<=10000): 
			string = '%4d' % nmol 
			newline=line[0:22] + string + line[26:]
		elif (nmol<=99999): 
			string = '{:<5}'.format(nmol)
			newline=line[0:22] + string + line[27:]
		else:
			string = '{:<6}'.format(nmol)
			newline=line[0:22] + string + line[28:]
			#print "MORE THAN 100000 ATOMS"
		pdbfile.write(newline)
		line=infile.readline()
	infile.close()
	pdbfile.close() 


if __name__=='__main__':
	main()