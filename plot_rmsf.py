#!/usr/bin/env python
# calculate RMSF plots for MD trajectories
# written by Jose Flores-Canales, 22/02/2017

from __future__ import print_function
import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from scipy import stats

def moving_test(x,y,size): 
	p = []
	for i in xrange(0,len(x),1):
		x_temp = x[i:i+size]
		y_temp = y[i:i+size]	
		if len(x_temp) != size or len(y_temp) != size:
			continue
		else:
			#f_temp, p_temp = stats.f_oneway(x_temp,y_temp)  
			f_temp, p_temp = stats.ttest_ind(x_temp,y_temp,equal_var=True)	
			p.append(p_temp)
	n = len(p)
	x = list(xrange(17,17+n))
	p = np.array(p)
	for i in xrange(n):
		if p[i] < 0.01:
			p[i] = 0.25 
		else:
			p[i] = np.nan
	return x,p

def average(x):
	assert len(x) > 0
	#print( x)
	return float(np.sum(x)) / float(len(x))

def pearson_def(x, y):
	assert len(x) == len(y)
	n = len(x)
	assert n > 0
	avg_x = average(x)
	avg_y = average(y)
	diffprod = 0.
	xdiff2 = 0.
	ydiff2 = 0.
	for idx in range(n):
		xdiff = x[idx] - avg_x
		ydiff = y[idx] - avg_y
		diffprod += xdiff * ydiff
		xdiff2 += xdiff * xdiff
		ydiff2 += ydiff * ydiff

	return diffprod / np.sqrt(xdiff2 * ydiff2)

def rmsf(traj):
	ave_coord1, ave_sqr_coord1 = get_average_structure(traj)
	rmsf=np.sqrt(np.sum(ave_sqr_coord1-ave_coord1*ave_coord1,axis=1))
	bfactor = rmsf*rmsf*8.*np.pi*np.pi/3.
	dict_rmsf = {}
	dict_bfactor = {}
	i = 17
	for r,b in zip(rmsf,bfactor):
		dict_rmsf[i] = r
		dict_bfactor[i] = b 
		i+=1
	return dict_rmsf, rmsf, dict_bfactor, bfactor

def read_pdb(pdb,chain):
	with open(pdb,'r') as f:
		lines = f.readlines()
	bfactor = {}
	for i in xrange(17,367):
		bfactor[i] = np.nan 

	for line in lines:
		flag = line[16:17]==' ' or line[16:17]=='A'
		if (flag and line[0:6]=='ATOM  ' and line[12:16].strip() == 'CA' and line[21:22]==chain):
			bfactor[int(line[22:26].strip())] = float(line[60:66].strip())
	return bfactor, np.array([bfactor[i] for i in xrange(17,367)])

def get_average_structure(traj):
	n_frames = traj.n_frames 
	n_atoms = traj.n_atoms
	average = np.zeros((n_atoms,3))
	average2 = np.zeros((n_atoms,3))
	for frame in xrange(n_frames):
		average += traj.xyz[frame]*10. 
		average2 += traj.xyz[frame]*traj.xyz[frame]*100. 
	return average/n_frames,average2/n_frames

t1 = md.load('/net/jfloresc/kinesin_Model1/run1_T310_c0.1_v0.6_long/T0EG5.dcd',top='ADP65_ATP35_b0.2_ftrue_T310_16_restart/T0EG5.pdb')
t2 = md.load('/net/jfloresc/kinesin_Model1/ADP_like/run1_T310_c0.1_v0.6_long/T0EG5.dcd',top='ADP65_ATP35_b0.2_ftrue_T310_16_restart/T0EG5.pdb')
pdbref = md.load('ADP65_ATP35_b0.2_ftrue_T310_16_restart/T0EG5.pdb')
topology = pdbref.topology
atom_to_keep = [a.index for a in topology.atoms if a.name == 'CA']
print(len(atom_to_keep))
t1.restrict_atoms(atom_to_keep)
t2.restrict_atoms(atom_to_keep)
pdbref.restrict_atoms(atom_to_keep)
pdbref.xyz, temp = get_average_structure(t1)
t1.superpose(pdbref)
pdbref.xyz, temp = get_average_structure(t2)
t2.superpose(pdbref)
dict_rmsf_t1, rmsf_t1,dict_bfactor_t1,bfactor_t1 = rmsf(t1)
dict_rmsf_t2, rmsf_t2,dict_bfactor_t2,bfactor_t2 = rmsf(t2)
dict_b1, xray_b1 = read_pdb('3HQD.pdb','A')
dict_b2, xray_b2 = read_pdb('1II6.pdb','A')
temp_sim_b1, temp_sim_b1_no = [],[]
temp_xrayb1, temp_xrayb1_no = [],[]
temp_sim_b2, temp_sim_b2_no = [],[]
temp_xrayb2, temp_xrayb2_no = [],[]
b_xrayb1, b_xrayb2 = [],[]
sim_b1, sim_b2 = [],[]
rmsf_b1, rmsf_b2 = [],[]

factor=7.
for i in xrange(17,367):
	if (i >=271 and i <=280) or i>=350: #== 366: 
		b_xrayb2.append(np.nan)
		sim_b2.append(np.nan)
		rmsf_b2.append(dict_rmsf_t2[i])
		continue
	elif i==88 or i ==141 or i==205: 
		temp_sim_b2.append(dict_bfactor_t2[i])
		temp_xrayb2.append(dict_b2[i])

		b_xrayb2.append(dict_b2[i])
		sim_b2.append(dict_bfactor_t2[i]/factor)
		rmsf_b2.append(dict_rmsf_t2[i])
	else:
		temp_sim_b2.append(dict_bfactor_t2[i])
		temp_xrayb2.append(dict_b2[i])
		temp_sim_b2_no.append(dict_bfactor_t2[i])
		temp_xrayb2_no.append(dict_b2[i])

		b_xrayb2.append(dict_b2[i])
		sim_b2.append(dict_bfactor_t2[i]/factor)
		rmsf_b2.append(dict_rmsf_t2[i])

for i in xrange(17,367):
	temp_sim_b1.append(dict_bfactor_t1[i])
	temp_xrayb1.append(dict_b1[i])
	b_xrayb1.append(dict_b1[i])
	sim_b1.append(dict_bfactor_t1[i]/factor)
	rmsf_b1.append(dict_rmsf_t1[i])
	if i==88 or i ==141 or i==205: 
		continue
	else:
		temp_sim_b1_no.append(dict_bfactor_t1[i])
		temp_xrayb1_no.append(dict_b1[i])
###############################################################################
p1_coef = stats.pearsonr(xray_b1,bfactor_t1)
p2_coef = stats.pearsonr(xray_b2,bfactor_t2)

p1_coef1 = pearson_def(temp_xrayb1,temp_sim_b1)
p2_coef2 = pearson_def(temp_xrayb2,temp_sim_b2)
p1_coef1_no = pearson_def(temp_xrayb1_no,temp_sim_b1_no)
p2_coef2_no = pearson_def(temp_xrayb2_no,temp_sim_b2_no)
print("pearson ATP",p1_coef1, "pearson ATP no P88, H141, H205",p1_coef1_no)
print("pearson ADP",p2_coef2, "pearson ADP no P88, H141, H205",p2_coef2_no)

x_p,p_value = moving_test(rmsf_b1,rmsf_b2,3) 
for i,j in zip(x_p,p_value):
	if not np.isnan(j):
		print(i,j)
###############################################################################


plt.rcParams.update({u'legend.fontsize': 7})
plt.rcParams.update({u'font.size': 10})
fig = plt.figure(figsize=(4,6))
ax = fig.add_subplot(1, 1, 1)
ax.spines['top'].set_color('none')
ax.spines['bottom'].set_color('none')
ax.spines['left'].set_color('none')
ax.spines['right'].set_color('none')
ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
x = np.array(list(range(17,367)))

ax1 = fig.add_subplot(2, 1, 1)
line1, = ax1.plot(x,sim_b1,lw=1.5,color='k',label='Sim.')
line2, = ax1.plot(x,b_xrayb1,'b--',lw=1.5,label='Exp.')
arr=ax1.arrow(88,48,0,-7.5,head_width=5,head_length=2.,fc='k',ec='k')
arr=ax1.arrow(141,48,0,-7.5,head_width=5,head_length=2.,fc='k',ec='k')
arr=ax1.arrow(205,48,0,-7.5,head_width=5,head_length=2.,fc='k',ec='k')
plt.text(0.85, 0.75,'ATP', ha='center', va='center', transform=ax1.transAxes,size=10)
plt.text(88, 52,'P88', ha='center', va='center', size=10)
plt.text(141, 52,'H141', ha='center', va='center', size=10)
plt.text(220, 52,'H205', ha='center', va='center', size=10)
leg = ax1.legend(handles=[line1, line2])
for legobj in leg.legendHandles:
	legobj.set_linewidth(1.5)
ax1.legend(loc='best',frameon=True)
ax2 = fig.add_subplot(2, 1, 2)
line3, = ax2.plot(x,sim_b2,lw=1.5,color='k',label='Sim.')
line4, = ax2.plot(x,b_xrayb2,'b--',lw=1.5,label='Exp.')
arr=ax2.arrow(88,78,0,-12.5,head_width=5,head_length=2.5,fc='k',ec='k')
arr=ax2.arrow(141,78,0,-12.5,head_width=5,head_length=2.5,fc='k',ec='k')
arr=ax2.arrow(205,78,0,-12.5,head_width=5,head_length=2.5,fc='k',ec='k')
plt.text(0.85, 0.75,'ADP', ha='center', va='center', transform=ax2.transAxes,size=10)
plt.text(88, 84,'P88', ha='center', va='center', size=10)
plt.text(141, 84,'H141', ha='center', va='center', size=10)
plt.text(205, 84,'H205', ha='center', va='center', size=10)
ax2.legend(loc='best',frameon=True)
leg = ax2.legend(handles=[line3, line4])
for legobj in leg.legendHandles:
	legobj.set_linewidth(1.5)

ax1.tick_params(axis='both', which='major', labelsize=7)
ax1.tick_params(axis='both', which='minor', labelsize=4)
ax2.tick_params(axis='both', which='major', labelsize=7)
ax2.tick_params(axis='both', which='minor', labelsize=4)
ax.set_xlabel('Residue Index')
ax.set_ylabel(r'Bfactor ($\AA$^2)')
nameout='bfactor_atp_adp'
fig.subplots_adjust(left=0.15,bottom=0.1, top=0.95, hspace=0.1)
fig.savefig(nameout+'.png',format='png', dpi=500)
####################################################################
fig = plt.figure(figsize=(6,3))
ax = fig.add_subplot(1,1,1)
line1, = ax.plot(x,rmsf_b1,lw=1.5,color='k',label='ATP')
line2, = ax.plot(x,rmsf_b2,'b--',lw=1.5,label='ADP')
#line3, = ax.plot(x_p,p_value,'r',lw=1.5,label='p')
leg = ax.legend(handles=[line1, line2])
for legobj in leg.legendHandles:
	legobj.set_linewidth(1.5)
ax.add_patch(patches.Rectangle((22,0.15),3,0.2))
ax.add_patch(patches.Rectangle((59,0.15),3,0.2))
ax.add_patch(patches.Rectangle((90,0.15),6,0.2))# 94
ax.add_patch(patches.Rectangle((113,0.15),3,0.2))
ax.add_patch(patches.Rectangle((122,0.15),6,0.2))# 125
ax.add_patch(patches.Rectangle(( 141,0.15),7,0.2))#142, 145
ax.add_patch(patches.Rectangle(( 149,0.15),3,0.2))
ax.add_patch(patches.Rectangle((176,0.15),4,0.2))#177
ax.add_patch(patches.Rectangle(( 222,0.15),3,0.2)) 
ax.add_patch(patches.Rectangle((226,0.15),3,0.2))
ax.add_patch(patches.Rectangle(( 263,0.15),3,0.2))
ax.add_patch(patches.Rectangle(( 268,0.15),3,0.2)) 
ax.add_patch(patches.Rectangle((272,0.15),9,0.2))# 273 274 275 276 277 278
ax.add_patch(patches.Rectangle(( 284,0.15),7,0.2)) #287 288 
ax.add_patch(patches.Rectangle((297,0.15),3,0.2)) 
ax.add_patch(patches.Rectangle((301,0.15),6,0.2)) # 304 
ax.add_patch(patches.Rectangle((322,0.15),6,0.2))# 325 
ax.add_patch(patches.Rectangle((332,0.15),5,0.2)) #333 334
ax.add_patch(patches.Rectangle(( 346,0.15),3,0.2))
ax.add_patch(patches.Rectangle(( 351,0.15),9,0.2))# 352 355 357
ax.add_patch(patches.Rectangle(( 363,0.15),4,0.2))# 364
ax.legend(loc='best',frameon=True)
ax.set_xlim([16,367])
ax.set_ylim([0,7])
ax.set_xticks(np.arange(17,368,15))
for label in ax.get_xticklabels()[1::2]:
	label.set_visible(False)
ax.set_xlabel('Residue Index')
ax.set_ylabel(r'RMSF ($\AA$)')
nameout='rmsf_atp_adp'
fig.subplots_adjust(left=0.15,bottom=0.15, top=0.95, hspace=0.1)
fig.savefig(nameout+'.png',format='png', dpi=500)
