#!/usr/bin/env python
# after solvation of a bilayer by tleap, water molecules should be removed from the sides
# this script assumes an input PDB solvated by tleap. 
# main variables in the __main__ function
#									 name = input.pdb
#									 prefix = mod
#									 selection = "resname LA PC" or "resid 1 to 624" (select lipids)
#									 radius = 6.0 (default value is 6.0 Angstroms, but it is recomended to test)
#									 control_depth_water = 6.0 control depth insertion of water in the lipid interface (default is 6.0 Angstroms)
#									 buffer_up = -1 (default -1 no modification. If assigned a value then that would be the top water thickness )
#		 buffer_down = -1 (default -1 no modification. If assigned a value then that would be the bottom water thickness)
#
#									./remove_waters.py -p input.pdb -o mod -s selection -r radius -d control_depth_water -u buffer_up -w buffer_down 
#  Examples:
#									./remove_waters.py -p input.pdb -o mod -s "resname LA PC" 
#									./remove_waters.py -p input.pdb -o mod -s "resname LA PC" -r 6.0 -d 6.0 
#									./remove_waters.py -p input.pdb -o mod -s "resid 1 to 624" -r 6.0 -d 6.0 -u 12 -w 10
#									./remove_waters.py -p input.pdb -o mod -s "resname LA PC" -u 10 -w 10
# jose flores-canales 12/26/2014
# REV 1, jose flores-canales 01/11/2016, implemented reading of standard PDB formats
#
# TODO
# Automatic reading of molecule templates from ff14SB and Lipid14 from $AMBERHOME or force-field directory
# Automatic addition of TER cards to input files 

import math
import optparse
import sys

def parse_cmdline(cmdlineArgs):
	parser = optparse.OptionParser("Usage: python remove_waters.py "
				 "[options]")

	parser.add_option("-p", "--pdbfile", action="store", dest="pdbFile")
	parser.add_option("-o", "--prefix", action="store", dest="mod")
	parser.add_option("-s", "--selection", action="store", dest="sel")
	parser.add_option("-r", "--radius", action="store", dest="radius",default=6.0,type="float",help = "[default: %default]")
	parser.add_option("-d", "--depth", action="store", dest="depth",default=6.0,type="float",help = "[default: %default]")
	parser.add_option("-u", "--bufferup", action="store", dest="buffup",default=-1,type="float",help = "water thickness [default: %default]")
	parser.add_option("-w", "--bufferdown", action="store", dest="buffdown",default=-1,type="float",help = "water thickness [default: %default]")

	opts, args = parser.parse_args(cmdlineArgs)
	pdbFile = opts.pdbFile
	mod = opts.mod
	sel = opts.sel
	rad = opts.radius
	depth = opts.depth
	buffup = opts.buffup
	buffdown = opts.buffdown

	if (pdbFile == None) or (mod == None) or (sel == None):
		parser.print_help()
		exit()
	return pdbFile, mod,sel,rad,depth,buffup,buffdown

class MainException(Exception):
	def __init__(self, msg):
		self.msg = msg
	def __str__(self):
		return self.msg

def getX(line):
	try:
		float_n = float(line[30:38])
		return float_n
	except:
		raise MainException("atom column coordinates in PDB file are not well formatted")
	
def getY(line):
	try:
		float_n = float(line[38:46])
		return float_n
	except:
		raise MainException("atom column coordinates in PDB file are not well formatted")

def getZ(line):
	try:
		float_n = float(line[46:54])
		return float_n
	except:
		raise MainException("atom column coordinates in PDB file are not well formatted")
	
def isAtom(line):
	if (line[0:6]=='ATOM  '):
		return True
	return False

def isCard(line):
	nums = line[0:3].strip()
	cards = ['TER','END']
	if (nums in cards):
		return True
	return False

def parseSel(selection):
# very basic and rough parsing only admits resname or resid
	readresname = 0
	identifiers = []
	readcolumn = -1
	nums = selection.split()
	if ('resname' in nums):
		for nums_i in nums:
			if (nums_i == 'resname'):
				readcolumn = 3
			else:
				identifiers.append(nums_i)
		return readcolumn, identifiers	
	if ('resid' in nums):
		for nums_i in nums:
			if (nums_i=='resid'):
				readcolumn = 4
			elif (nums_i != 'to' and nums_i !='TO'):
				identifiers.append(int(nums_i))
		if (identifiers[0]> identifiers[1]):
			raise MainException("wrong order of selection")
		return readcolumn, identifiers	
	raise MainException("wrong selection")

def isSelection(line,column,identifiers):
	if (column == 3):
		if (line[17:20].strip() in identifiers):
			return True
		return False
	elif (column == 4): 
		resids = [i for i in xrange(identifiers[0],identifiers[1]+1)]
# standard columns are 22 to 26. Column 27 is taken here for residues > 9999 (consistent with tleap)
		if (int(line[22:27]) in resids):
			return True
		return False
	
def findmaxmin(name_file, selection):
	infile = open(name_file,'r')
	column, identifiers = parseSel(selection)
	X, Y, Z = [], [], []
	for line in infile:
		if isAtom(line):
			if (isSelection(line,column,identifiers)):
				X.append(getX(line))
				Y.append(getY(line))
				Z.append(getZ(line))
	max_x = max(X)
	min_x = min(X)
	max_y = max(Y)
	min_y = min(Y)
	max_z = max(Z)
	min_z = min(Z)
	print "Maximum of selection"
	print "X: ", max_x, "Y: ", max_y, "Z: ",max_z
	print "Minimum of selection"
	print "X: ", min_x, "Y: ", min_y, "Z: ",min_z
	infile.close()
	return [max_x,max_y,max_z], [min_x,min_y,min_z]

def isInsideBoxBuffer(max_coord,min_coord,tempx,tempy,tempz,radius,thick_up,thick_down,depth):
# modify this function or construct a new one for more complicated solvation boxes
	vx,vy,vz,vzbuf=[],[],[],[]
	vx=[False for i in tempx if (i > (max_coord[0]-radius) or i < (min_coord[0]+radius))]
	vy=[False for i in tempy if (i > (max_coord[1]-radius) or i < (min_coord[1]+radius))]
	vz=[False for i in tempz if (i < (max_coord[2]-depth) and i > (min_coord[2]+depth))]
	vzbuf= [False for i in tempz if (i > (max_coord[2]+thick_up) or i < (min_coord[2]-thick_down))]
	
	if (len(vx)==0 and len(vy)==0 and len(vz)==0 and len(vzbuf)==0):
		return True
	return False
	
def isInsideBox(max_coord,min_coord,tempx,tempy,tempz,radius,depth):
# modify this function or construct a new one for more complicated solvation boxes
	vx,vy,vz=[],[],[]
	vx=[False for i in tempx if (i > (max_coord[0]-radius) or i < (min_coord[0]+radius))]
	vy=[False for i in tempy if (i > (max_coord[1]-radius) or i < (min_coord[1]+radius))]
	#vxy=[False for i,j in tempx,tempy if ((i > max_coord[0] and j > max_coord[1]) or (i < min_coord[0] and j < min_coord[1]) )]
	vz=[False for i in tempz if (i < (max_coord[2]-depth) and i > (min_coord[2]+depth))]
	if (len(vx)==0 and len(vy)==0 and len(vz)==0):
		return True
	return False

def isWater(string):
	nums = string[17:20].strip()
	if (nums in ['HOH', 'WAT']):
		return True
	return False

def findmolwaters(name_file):
	infile = open(name_file,'r')
	n = 0
	for line in infile:
		if isAtom(line):
			if (isWater(line)):
				n+=1
	infile.close()
	return int(n/3)

def removeWaterBuf(name_file,selection,max_coord,min_coord,radius,thick_up,thick_down,depth):
	infile = open(name_file,'r')
	j,flag_print = 0,0
	lines=[]
	tempx, tempy,tempz,templines=[],[],[],[]
	for line in infile:
		if isAtom(line):
			if (isWater(line)):
				tempx.append(getX(line))
	tempy.append(getY(line))
	tempz.append(getZ(line))
	templines.append(line)
	j+=1
	if (j==3):
# modify or add a new function instead of isInsideBox for more complicated solvation boxes
		if (isInsideBoxBuffer(max_coord,min_coord,tempx,tempy,tempz,radius,thick_up,thick_down,depth)):
			[lines.append(i) for i in templines]
		else:
			flag_print = 1
		j = 0
		tempx, tempy,tempz,templines=[],[],[],[]
			else:
				lines.append(line)
		elif (isCard(line) and flag_print==0):
			lines.append(line)
		elif (isCard(line) and flag_print==1):
			flag_print=0
	infile.close()
	return lines

def removeWater(name_file,selection,max_coord,min_coord,radius,depth):
	infile = open(name_file,'r')
	j,flag_print = 0,0
	lines=[]
	tempx, tempy,tempz,templines=[],[],[],[]
	for line in infile:
		if isAtom(line):
			if (isWater(line)):
				tempx.append(getX(line))
	tempy.append(getY(line))
	tempz.append(getZ(line))
	templines.append(line)
	j+=1
	if (j==3):
# modify or add a new function instead of isInsideBox for more complicated solvation boxes
		if (isInsideBox(max_coord,min_coord,tempx,tempy,tempz,radius,depth)):
			[lines.append(i) for i in templines]
		else:
			flag_print = 1
		j = 0
		tempx, tempy,tempz,templines=[],[],[],[]
			else:
				lines.append(line)
		elif (isCard(line) and flag_print==0):
			lines.append(line)
		elif (isCard(line) and flag_print==1):
			flag_print=0
	infile.close()
	return lines

def addPrefix(file,tag):
	return tag+"."+file
	
def print_object(file,text):
	file_out = open(file,'w')
	file_out.writelines(text)
	file_out.close()
			
def main_func(file, sel,tag,radius,thick_up,thick_down,depth):
	waters = findmolwaters(file)
	print "Initial number of waters: ", waters
	max_sel, min_sel = findmaxmin(file,sel)
	
# max_sel and min_sel could be defined manually  
	if (thick_up == -1 and thick_down == -1):
		toprint=removeWater(file,sel,max_sel,min_sel,radius,depth)
	elif (thick_up > -1 and thick_down > -1):
		toprint=removeWaterBuf(file,sel,max_sel,min_sel,radius,thick_up,thick_down,depth)
	else:
		raise MainException("wrong values of water thickness parameters")  

	file_out = addPrefix(file,tag)
	print_object(file_out,toprint)
	waters = findmolwaters(file_out)
	print "Final number of waters: ", waters

###########################################################
# __main__ module
# variables defined here are global
###########################################################
if __name__ == "__main__":
 
	_input_pdb, _prefix, _string_sel,_r, _depth_water,_bu, _bw=parse_cmdline(sys.argv[1:])
	
	main_func(_input_pdb, _string_sel,_prefix,_r, _bu, _bw,_depth_water)
	
	exit()
