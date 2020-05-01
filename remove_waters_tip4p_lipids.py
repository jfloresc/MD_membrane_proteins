#!/usr/bin/python
# after solvation of a bilayer by tleap, water molecules should be removed from the sides
# this script assumes an input PDB solvated by tleap. 
# main variables in the __main__ function
#									 name = input.pdb
#									 prefix = mod
#									 selection = "resname LA PC" or "resid 1 to 624" (select lipids)
#									 radius = 6.0 (default value is 6.0 Angstroms, but it is recomended to test)
#									 control_depth_water = 6.0 (default value is 6.0 Angstroms, control depth insertion of water in the lipid interfacei)
#									 buffer_up = -1 (default -1). If assigned a value then that would be the top water thickness 
#									 buffer_down = -1 (default -1). If assigned a value then that would be the bottom water thickness
#									 fit_lipids = 1 or 0 (default 0). If 1 then lipids outside the water box will be removed 
#
#									./remove_waters.py -p input.pdb -o mod -s selection -r radius -d control_depth_water -u buffer_up -w buffer_down 
#  Examples:
#									./remove_waters.py -p input.pdb -o mod -s "resname LA PC" 
#									./remove_waters.py -p input.pdb -o mod -s "resname LA PC" -r 6.0 -d 6.0 
#									./remove_waters.py -p input.pdb -o mod -s "resid 1 to 624" -r 6.0 -d 6.0 -u 12 -w 10
#									./remove_waters.py -p input.pdb -o mod -s "resname LA PC" -u 10 -w 10
#									./remove_waters.py -p input.pdb -o mod -s "resname LA PC" -f 1
#									./remove_waters.py -p input.pdb -o mod -s "resid 1 to 2000" -n  12000
#									./remove_waters.py -p input.pdb -o mod -s "WAT" -n  12000
# jose flores-canales 12/26/2014
# REV 1, jose flores-canales 01/11/2016, implemented reading of standard PDB formats
# REV 2, jose flores-canales 03/21/2016, added option to remove lipids outside of the waterbox XY dimensions
# REV 3, jose flores-canales 04/15/2020, adapted for TIP4P water based models
#
# TODO
# Automatic reading of molecule templates from ff14SB and Lipid14 from $AMBERHOME or force-field directory
# Automatic addition of TER cards to input files 

import math
import optparse
import sys

def parse_cmdline(cmdlineArgs):
	usage = "Usage: python remove_waters.py [options]"
	parser = optparse.OptionParser(usage=usage)

	parser.add_option("-p", "--pdbfile", action="store", dest="pdbFile")
	parser.add_option("-o", "--prefix", action="store", dest="mod")
	parser.add_option("-s", "--selection", action="store", dest="sel")
	parser.add_option("-r", "--radius", action="store", dest="radius",default=6.0, type="float", help = "[default: %default]")
	parser.add_option("-d", "--depth", action="store", dest="depth",default=6.0, type="float", help = "[default: %default]")
	parser.add_option("-u", "--bufferup", action="store", dest="buffup",default=-1, type="float", help = "water thickness [default: %default]")
	parser.add_option("-w", "--bufferdown", action="store", dest="buffdown",default=-1, type="float", help = "water thickness [default: %default]")
	parser.add_option("-f", "--fitlipids", action="store", dest="fitlipid",default=0, type="int", help = "flag (1 or 0) to fit water box XY dimensions to lipid patch [default: %default]")
	parser.add_option("-n", "--keep", action="store", dest="keep",default=-1, type="int", help = "number of waters to keep [default: %default]")

	opts, args = parser.parse_args(cmdlineArgs)
	pdbFile = opts.pdbFile
	mod = opts.mod
	sel = opts.sel
	rad = opts.radius
	depth = opts.depth
	buffup = opts.buffup
	buffdown = opts.buffdown
	fitlipid = opts.fitlipid
	keep_nwaters = opts.keep

	if (pdbFile == None) or (mod == None) or (sel == None):
		parser.print_help()
		exit()
	return pdbFile, mod, sel, rad, depth, buffup, buffdown, fitlipid, keep_nwaters


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
	cards = ['TER', 'END']
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

def isSelection(line, selection):
	if selection == "WAT":
		if isWater(line): 
			return True
		return False
	column, identifiers = parseSel(selection)
	if (column == 3):
		if (line[17:20].strip() in identifiers):
			return True
		return False
	elif (column == 4): 
		resids = [i for i in range(identifiers[0],identifiers[1]+1)]
# standard columns are 22 to 26. Column 27 is taken here for residues > 9999 (consistent with tleap)
		if (int(line[22:27]) in resids):
			return True
		return False
	
def findmaxmin(name_file, selection):
	with open(name_file,'r') as f:
		infile = f.readlines()

	X, Y, Z = [], [], []
	for line in infile:
		if isAtom(line):
			if (isSelection(line, selection)):
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
	return [max_x,max_y,max_z], [min_x,min_y,min_z]

def isLipidInsideWaterBox(tempx,tempy,tempz,max_coord,min_coord,radius):
# modify this function or construct a new one for more complicated solvation boxes
	vx,vy,vxy=[],[],[]
	radius = 8 
	vx=[False for i in tempx if (i > (max_coord[0]+radius) or i < (min_coord[0]-radius))]
	vy=[False for i in tempy if (i > (max_coord[1]+radius) or i < (min_coord[1]-radius))]
	#vxy=[False for i,j in tempx,tempy if ((i > (max_sel_w[0] - radius) and j > (max_sel_w[1]-radius)) or (i < (min_sel_w[0]+radius) and j < (min_sel_w[1]+radius)) )]
	if (len(vx)==0 and len(vy)==0 ):
	#if len(vxy)==0 :
		return True
	return False

def isInsideBoxBuffer(max_coord, min_coord, tempx, tempy, tempz, radius, thick_up, thick_down, depth):
# modify this function or construct a new one for more complicated solvation boxes
	vx, vy, vz, vzbuf = [], [], [], []
	vx = [False for i in tempx if (i > (max_coord[0]-radius) or i < (min_coord[0]+radius))]
	vy = [False for i in tempy if (i > (max_coord[1]-radius) or i < (min_coord[1]+radius))]
	vz = [False for i in tempz if (i < (max_coord[2]-depth) and i > (min_coord[2]+depth))]
	vzbuf = [False for i in tempz if (i > (max_coord[2]+thick_up) or i < (min_coord[2]-thick_down))]
	
	if (len(vx) ==0 and len(vy) ==0 and len(vz) ==0 and len(vzbuf) == 0):
		return True
	return False
	
def isInsideBox(max_coord, min_coord, tempx, tempy, tempz, radius, depth):
# modify this function or construct a new one for more complicated solvation boxes
	vx, vy, vz = [], [], []
	vx = [False for i in tempx if (i > (max_coord[0]-radius) or i < (min_coord[0]+radius))]
	vy = [False for i in tempy if (i > (max_coord[1]-radius) or i < (min_coord[1]+radius))]
	#vxy=[False for i,j in tempx,tempy if ((i > max_coord[0] and j > max_coord[1]) or (i < min_coord[0] and j < min_coord[1]) )]
	vz = [False for i in tempz if (i < (max_coord[2]-depth) and i > (min_coord[2]+depth))]
	if (len(vx) == 0 and len(vy) == 0 and len(vz) == 0):
		return True
	return False

def isWater(string):
	nums = string[17:20].strip()
	if (nums in ['HOH', 'WAT']):
		return True
	return False

def findmolwaters(name_file):
	with open(name_file,'r') as infile:
		lines_infile = infile.readlines()

	n = 0
	for line in lines_infile:
		if isAtom(line):
			if (isWater(line)):
				n+=1
	# divided by 3 for tip3p or 4 for tip4p based models
	return int(n/4)

def removeWaterBuf(name_file, selection, max_coord, min_coord, radius, thick_up, thick_down, depth, fit_lipid_flag, max_sel_w, min_sel_w):
	with open(name_file,'r') as infile:
		lines_infile = infile.readlines()
	j,flag_print = 0,0
	line_count = 0
	lines=[]
	tempx, tempy, tempz, templines = [], [], [], []
	for line in lines_infile:
		if isAtom(line):
			if (isWater(line)):
				tempx.append(getX(line))
				tempy.append(getY(line))
				tempz.append(getZ(line))
				templines.append(line)
				j+=1
				# j = 3 for TIP3P, 4 for TIP4P
				if (j==4):
				# modify or add a new function instead of isInsideBox for more complicated solvation boxes
					if (isInsideBoxBuffer(max_coord, min_coord, tempx, tempy, tempz, radius, thick_up, thick_down, depth)):
						[lines.append(i) for i in templines]
					else:
						flag_print = 1
				j = 0
				tempx, tempy, tempz, templines = [], [], [], []
# adding a new function
			elif (isSelection(line, selection) and fit_lipid_flag == 1):
				tempx.append(getX(line))
				tempy.append(getY(line))
				tempz.append(getZ(line))
				templines.append(line)
				if (isCard(lines_infile[line_count+1])):
					if (isLipidInsideWaterBox(tempx, tempy, tempz, max_sel_w, min_sel_w, radius)):
						[lines.append(i) for i in templines]
					else:
						flag_print = 1
					tempx, tempy, tempz, templines = [], [], [], []
# end of adding a neew function
			else:
				lines.append(line)
		elif (isCard(line) and flag_print==0):
			lines.append(line)
# don't print the TER card of the removed water molecule
		elif (isCard(line) and flag_print==1):
			flag_print=0
		line_count += 1
	return lines


def removeWater(name_file, selection, max_coord, min_coord, radius, depth, fit_lipid_flag, max_sel_w, min_sel_w):
	with open(name_file,'r') as infile:
		lines_infile = infile.readlines()
	j, flag_print = 0, 0
	line_count = 0
	lines=[]
	tempx, tempy, tempz, templines = [], [], [], []
	for line in lines_infile:
		if isAtom(line):
			if (isWater(line)):
				tempx.append(getX(line))
				tempy.append(getY(line))
				tempz.append(getZ(line))
				templines.append(line)
				j+=1
				# j = 3 for TIP3P, 4 for TIP4P
				if (j==4):
				# modify or add a new function instead of isInsideBox for more complicated solvation boxes
					if (isInsideBox(max_coord, min_coord, tempx, tempy,tempz, radius, depth)):
						[lines.append(i) for i in templines]
					else:
						flag_print = 1
					j = 0
					tempx, tempy, tempz, templines = [], [], [], []
# adding a new function
			elif (isSelection(line,selection) and fit_lipid_flag == 1):
				tempx.append(getX(line))
				tempy.append(getY(line))
				tempz.append(getZ(line))
				templines.append(line)
				if (isCard(lines_infile[line_count+1])):
					if (isLipidInsideWaterBox(tempx, tempy, tempz, max_sel_w, min_sel_w, radius)):
						[lines.append(i) for i in templines]
					else:
						flag_print = 1
					tempx, tempy, tempz, templines = [], [], [], []
# end of adding a neew function
			else:
				lines.append(line)
		elif (isCard(line) and flag_print==0):
			lines.append(line)
		elif (isCard(line) and flag_print==1):
			flag_print=0
		line_count += 1
	return lines


def downsizeWater(name_file, max_sel_w, min_sel_w, diff_w, delta):
	with open(name_file,'r') as infile:
		lines_infile = infile.readlines()
	dt = delta
	while True:
		j, flag_print = 0, 0
		line_count = 0
		lines = []
		tempx, tempy, tempz, templines = [], [], [], []
		removed_w = 0
		for line in lines_infile:
			if isAtom(line):
				if (isWater(line)):
					tempx.append(getX(line))
					tempy.append(getY(line))
					tempz.append(getZ(line))
					templines.append(line)
					j+=1
					# j = 3 for TIP3P, 4 for TIP4P
					if (j==4):
					# modify or add a new function instead of isInsideBox for more complicated solvation boxes
						if isInsideBox(max_sel_w, min_sel_w, tempx, tempy, tempz, dt, dt) or removed_w >= diff_w:
							[lines.append(i) for i in templines]
						else:
							removed_w += 1
							flag_print = 1
						j = 0
						tempx, tempy, tempz, templines = [], [], [], []
				else:
					lines.append(line)
			elif (isCard(line) and flag_print==0):
				lines.append(line)
			elif (isCard(line) and flag_print==1):
				flag_print=0
			line_count += 1
		if removed_w  == diff_w:
			print "Number of waters removed", removed_w
			return lines
		dt += delta


def addPrefix(filename, tag):
	return tag + "." + filename
	
def print_object(filename, text):
	with open(filename, 'w') as file_out:
		file_out.writelines(text)


def main_func(filename, sel, tag, radius, thick_up, thick_down, depth, fit_lipid_flag, keep_nwaters):
	waters = findmolwaters(filename)
	print "Initial number of waters: ", waters
	max_sel, min_sel = findmaxmin(filename, sel)

	water_res = "WAT"
	max_sel_w, min_sel_w = findmaxmin(filename, water_res)
	
	if keep_nwaters != -1:
		print "Warning: Keep # waters option -n is being used, other options are ignored"

		diff_w = waters - keep_nwaters
		if diff_w < 0:
			print "Number of waters to keep larger than #waters in the pdb file"
			exit(0)

		delta = 0.1
		toprint = downsizeWater(filename, max_sel_w, min_sel_w, diff_w, delta)
		
# max_sel and min_sel could be defined manually  
	elif (thick_up == -1 and thick_down == -1 and keep_nwaters == -1):
		toprint = removeWater(filename, sel, max_sel, min_sel, radius, depth, fit_lipid_flag, max_sel_w, min_sel_w)
	elif (thick_up > -1 and thick_down > -1 and keep_nwaters == -1):
		toprint = removeWaterBuf(filename, sel, max_sel, min_sel, radius, thick_up, thick_down, depth, fit_lipid_flag, max_sel_w, min_sel_w)
	else:
		raise MainException("wrong values of water thickness parameters")  

	file_out = addPrefix(filename, tag)
	print_object(file_out, toprint)
	waters = findmolwaters(file_out)
	print "Final number of waters: ", waters

###########################################################
# __main__ module
# variables defined here are global
###########################################################
if __name__ == "__main__":
 
	_input_pdb, _prefix, _string_sel,_r, _depth_water,_bu, _bw, _fit_lipid, _keep_nwaters = parse_cmdline(sys.argv[1:])
	
	main_func(_input_pdb, _string_sel,_prefix,_r, _bu, _bw,_depth_water,_fit_lipid, _keep_nwaters)
	
	exit()
