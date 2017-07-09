#!/usr/bin/env python
# loads multiple trajectories and aligns them relative to a reference pdb structure
# prints trajectories in the first 2 PC vectors
# plots 2D plots on the first 2 PC vectors
# plots variance
# plots RMSD mode weighted plots of the first 2 PC vectors 
# written by Jose Flores-Canales 08/2016


from __future__ import print_function
import mdtraj as md
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import numpy as np

def correlation_numpy(A, n_frames, n_atoms):
  #3000, 54, 3
  print("size of traj: ", A.shape)
  MA=np.mean(A, axis=0)
  print("mean(A): ", MA.shape)
  M = (A-MA) # subtract the mean (along columns)
#14000, 1050
  print("M=A-MA.T: ", M.shape)
#1050,
  M = M.reshape(n_frames,n_atoms*3)
  MT = np.memmap('mt.npy',dtype=np.float64,mode='w+',shape=(n_atoms*3,n_frames))
  MT[:] = M.T[:]
  cov_nn = np.memmap('covnn.npy',dtype=np.float64,mode='w+',shape=(n_atoms*3,n_atoms*3))
  cov_nn[:] = MT.dot(M)[:]/n_frames
  #cov_nn = np.tensordot(M,MT,axes=(-1,-1))/n_frames
  diag = np.copy(cov_nn.diagonal())
  for i in xrange(n_atoms*3):
    for j in xrange(n_atoms*3):
      cov_nn[i,j] = cov_nn[i,j]/np.sqrt(diag[i]*diag[j])
  return cov_nn

def slice_before_load(netcdf, prmtop):
  topology = md.load_topology(prmtop)
  with md.open(netcdf) as f:
    f.seek(START)
    t = f.read_as_traj(
      topology, n_frames=STOP-START, stride=STRIDE
    )
    print(t)
    return t

def get_average_structure(traj):
  n_frames = traj.n_frames 
  n_atoms = traj.n_atoms
  average = np.zeros((n_atoms,3))
  for frame in xrange(n_frames):
    average += traj.xyz[frame] 
  return average/n_frames

START = 1 
STOP = 65000
STRIDE = 10
T1 = md.load('/scratch/jfloresc/Transitions/ADP65_ATP35_b0.2_ftrue_T310_16_restart/T0EG5.dcd',top='/net/jfloresc/kinesin_Model1/traj6/T0EG5.pdb')
T2 = md.load('/scratch/jfloresc/Transitions/ADP65_ATP35_b0.2_ftrue_T310_17_restart/T0EG5.dcd',top='/net/jfloresc/kinesin_Model1/traj6/T0EG5.pdb')
T3 = md.load('/scratch/jfloresc/Transitions/ADP65_ATP35_b0.2_ftrue_T310_18_restart/T0EG5.dcd',top='/net/jfloresc/kinesin_Model1/traj6/T0EG5.pdb')
T4 = md.load('/scratch/jfloresc/Transitions/ADP65_ATP35_b0.2_ftrue_T310_19_restart/T0EG5.dcd',top='/net/jfloresc/kinesin_Model1/traj6/T0EG5.pdb')
T5 = md.load('/scratch/jfloresc/Transitions/ADP65_ATP35_b0.2_ftrue_T310_20_restart/T0EG5.dcd',top='/net/jfloresc/kinesin_Model1/traj6/T0EG5.pdb')
T6 = md.load('/scratch/jfloresc/Transitions/ADP65_ATP35_b0.2_ftrue_T310_21_restart/T0EG5.dcd',top='/net/jfloresc/kinesin_Model1/traj6/T0EG5.pdb')
T7 = md.load('/scratch/jfloresc/Transitions/ADP65_ATP35_b0.2_ftrue_T310_22_restart/T0EG5.dcd',top='/net/jfloresc/kinesin_Model1/traj6/T0EG5.pdb')
T8 = md.load('/scratch/jfloresc/Transitions/ADP65_ATP35_b0.2_ftrue_T310_23_restart/T0EG5.dcd',top='/net/jfloresc/kinesin_Model1/traj6/T0EG5.pdb')
T9 = md.load('/scratch/jfloresc/Transitions/ADP65_ATP35_b0.2_ftrue_T310_24_restart/T0EG5.dcd',top='/net/jfloresc/kinesin_Model1/traj6/T0EG5.pdb')
T10 = md.load('/scratch/jfloresc/Transitions/ADP65_ATP35_b0.2_ftrue_T310_25_restart/T0EG5.dcd',top='/net/jfloresc/kinesin_Model1/traj6/T0EG5.pdb')
T1a = T1.join(T2)
T2a = T1a.join(T3)
T3a = T2a.join(T4)
T4a = T3a.join(T5)
T5a = T4a.join(T6)
T6a = T5a.join(T7)
T7a = T6a.join(T8)
T8a = T7a.join(T9)
T9a = T8a.join(T10)

##load a reference or initial PDB structure 'Any_Reference_Structure.pdb'
#pdbref = md.load('1-average.pdb')
pdbref = md.load('/net/jfloresc/kinesin_Model1/traj6/T0EG5.pdb')
atpref = md.load('/net/jfloresc/kinesin_Model1/traj6/T0EG5.pdb')
adpref = md.load('/net/jfloresc/kinesin_Model1/ADP_like/traj1_u/T0EG5.pdb')
topology = pdbref.topology

#atom_to_keep = topology.select('name CA and (resid 184 to 218 or resid 589 to 623 or resid 994 to 1028 or resid 1399 to 1433)')
atom_to_keep = topology.select('name CA')
#atom_to_keep = topology.select_atom_indices('heavy')
#atom_to_keep = topology.select('name CA and (resid 1 to 356)')

print("number of atoms: ", len(atom_to_keep))

T9a.restrict_atoms(atom_to_keep)
T1.restrict_atoms(atom_to_keep)
T2.restrict_atoms(atom_to_keep)
T3.restrict_atoms(atom_to_keep)
T4.restrict_atoms(atom_to_keep)
T5.restrict_atoms(atom_to_keep)
T6.restrict_atoms(atom_to_keep)
T7.restrict_atoms(atom_to_keep)
T8.restrict_atoms(atom_to_keep)
T9.restrict_atoms(atom_to_keep)
T10.restrict_atoms(atom_to_keep)
pdbref.restrict_atoms(atom_to_keep)
atpref.restrict_atoms(atom_to_keep)
adpref.restrict_atoms(atom_to_keep)

##"get average coordinates"
pdbref.xyz = get_average_structure(T9a)
rmsds = md.rmsd(T9a,pdbref,0)
##"get structure closest to the average coordinates"
pdbref.xyz = T9a.xyz[np.argmin(rmsds)] 
T9a.superpose(pdbref)
#pdbref.save_pdb('average_TMD_M3.pdb')
pdbref.save_pdb('average1.pdb')

T1.superpose(pdbref)
T2.superpose(pdbref)
T3.superpose(pdbref)
T4.superpose(pdbref)
T5.superpose(pdbref)
T6.superpose(pdbref)
T7.superpose(pdbref)
T8.superpose(pdbref)
T9.superpose(pdbref)
T10.superpose(pdbref)
atpref.superpose(pdbref)
adpref.superpose(pdbref)

print("Number of frames:",T9a.n_frames)
pca1 = PCA(n_components=10)
reduced_cartesian = pca1.fit_transform(T9a.xyz.reshape(T9a.n_frames, T9a.n_atoms * 3)*10.)

centered = T9a.xyz.reshape(T9a.n_frames,T9a.n_atoms*3)-pca1.mean_
centered_xyz = centered.reshape(T9a.n_frames,T9a.n_atoms,3)
corr_matrix = correlation_numpy(centered_xyz*10., T9a.n_frames, T9a.n_atoms)
print("dimensions of corr matrix ",corr_matrix.shape)
corr_matrix_sqr = corr_matrix*corr_matrix
print("dimensions of corr matrix sqr:",corr_matrix_sqr.shape)

# Print 
print("pca1.components_.shape", pca1.components_.shape) #modes, natoms*3 
traj2_centered = T1.xyz.reshape(T1.n_frames, T1.n_atoms*3)*10. - pca1.mean_
scores_1 = np.dot(pca1.components_, traj2_centered.T).T 

traj3_centered = T2.xyz.reshape(T2.n_frames, T2.n_atoms*3)*10. - pca1.mean_
scores_2 = np.dot(pca1.components_, traj3_centered.T).T 

traj4_centered = T3.xyz.reshape(T3.n_frames, T3.n_atoms*3)*10. - pca1.mean_
scores_3 = np.dot(pca1.components_, traj4_centered.T).T 

traj5_centered = T4.xyz.reshape(T4.n_frames, T4.n_atoms*3)*10. - pca1.mean_
scores_4 = np.dot(pca1.components_, traj5_centered.T).T 

traj6_centered = T5.xyz.reshape(T5.n_frames, T5.n_atoms*3)*10. - pca1.mean_
scores_5 = np.dot(pca1.components_, traj6_centered.T).T 

traj7_centered = T6.xyz.reshape(T6.n_frames, T6.n_atoms*3)*10. - pca1.mean_
scores_6 = np.dot(pca1.components_, traj7_centered.T).T 

traj8_centered = T7.xyz.reshape(T7.n_frames, T7.n_atoms*3)*10. - pca1.mean_
scores_7 = np.dot(pca1.components_, traj8_centered.T).T 

traj9_centered = T8.xyz.reshape(T8.n_frames, T8.n_atoms*3)*10. - pca1.mean_
scores_8 = np.dot(pca1.components_, traj9_centered.T).T 

traj10_centered = T9.xyz.reshape(T9.n_frames, T9.n_atoms*3)*10. - pca1.mean_
scores_9 = np.dot(pca1.components_, traj10_centered.T).T 

traj11_centered = T10.xyz.reshape(T10.n_frames, T10.n_atoms*3)*10. - pca1.mean_
scores_10 = np.dot(pca1.components_, traj11_centered.T).T 

atp_centered = atpref.xyz.reshape(1, atpref.n_atoms*3)*10. - pca1.mean_
scores_atp = np.dot(pca1.components_, atp_centered.T).T 
adp_centered = adpref.xyz.reshape(1, adpref.n_atoms*3)*10. - pca1.mean_
scores_adp = np.dot(pca1.components_, adp_centered.T).T 
#r1 = pca1.transform(t1.xyz.reshape(t1.n_frames, t1.n_atoms * 3))
#all_r.append(r1)

#print("print the two vectors: ", reduced_cartesian)
#print(r1[0])
#print(r1[1])
#print("element 0 :", reduced_cartesian[0])
#print("element 1 :", reduced_cartesian[1])

# recons = PCscores * Eigenvectors.T + Average

#########################################
# MODIFY INDEXES of reduced_cartesian and pca1.components_
# [:,0,np.newaxis] [np.newaxis,0,:] --> 1st PC
# [:,1,np.newaxis] [np.newaxis,1,:] --> 2nd PC
recons = np.dot(reduced_cartesian[:,0,np.newaxis] , pca1.components_[np.newaxis,0,:]).reshape(T9a.n_frames,T9a.n_atoms,3) + pca1.mean_.reshape(T9a.n_atoms,3) #pdbref.xyz#.reshape(t1.n_atoms*3)
T9a.xyz=recons
T9a.save('traj_pc1_heavy.dcd')
recons = np.dot(reduced_cartesian[:,1,np.newaxis] , pca1.components_[np.newaxis,1,:]).reshape(t9a.n_frames,t1.n_atoms,3) + pdbref.xyz#.reshape(t1.n_atoms*3)
#t1.xyz=recons
#t1.save('traj_pc2_all.pdb')
#recons = np.dot(reduced_cartesian[:,2,np.newaxis] , pca1.components_[np.newaxis,2,:]).reshape(t1.n_frames,t1.n_atoms,3) + pdbref.xyz#.reshape(t1.n_atoms*3)
#t1.xyz=recons
#t1.save('traj_pc3_all.pdb')
# keep adding more modes if needed
#################################################

# First trajectory is plotted on its PCs
all_r = []
all_r.append(reduced_cartesian)
all_t = []
all_t.append(T9a.time)
print("Dimensions of the reduced space: ", reduced_cartesian.shape)
print("Dimensions of the trajectory file:  ", T9a.xyz.shape)
print("Explained variance ratio:", pca1.explained_variance_ratio_)
plt.figure(1)
plt.scatter(all_r[0][:, 0], all_r[0][:,1], marker='x', c=all_t[0]/10.)
plt.scatter(scores_atp[:, 0], scores_atp[:,1], marker='o',c= 'darkorange',s=20)#c=all_t[a])
plt.scatter(scores_adp[:, 0], scores_adp[:,1], marker='d',c= 'darkgreen',s=20)#c=all_t[a])
plt.xlabel('PC1')
plt.ylabel('PC2')
#cbar = plt.colorbar()
#cbar.set_label('Time steps')
nameout='pca_traj1_all'
plt.savefig(nameout+'.png')

# Second trajectory is plotted over the PCs from the first trajectory
all_r2 = []
all_r2.append(scores_1)
all_r2.append(scores_2)
all_r2.append(scores_3)
all_r2.append(scores_4)
all_r2.append(scores_5)
all_r2.append(scores_6)
all_r2.append(scores_7)
all_r2.append(scores_8)
all_r2.append(scores_9)
all_r2.append(scores_10)
all_t = []
all_t.append(T1.time)
all_t.append(T2.time)
all_t.append(T3.time)
all_t.append(T4.time)
all_t.append(T5.time)
all_t.append(T6.time)
all_t.append(T7.time)
all_t.append(T8.time)
all_t.append(T9.time)
all_t.append(T10.time)
colors_list = ['#67001f','#b2182b','#d6604d','#f4a582','#fddbc7','#d1e5f0','#92c5de','#4393c3','#2166ac','#053061']

print("Dimensions of the reduced space of the 2nd traj.", scores_1.shape)
print("Dimensions of the 2nd trajectory file:  ", T1.xyz.shape)
plt.rcParams.update({u'legend.fontsize': 6})
plt.rcParams.update({u'font.size': 10})

plt.figure(figsize=(4,3))
symbols=['o','v','^','<','>','8','s','p','*','x']
for a in xrange(10):
  plt.scatter(all_r2[a][:, 0], all_r2[a][:,1], marker=symbols[a],c= colors_list[a],s=10,edgecolors='none',label='T%s'%(a+1))#c=all_t[a])
plt.scatter(scores_atp[:, 0], scores_atp[:,1], marker='<',c= 'darkorange',s=40,edgecolors='none')#c=all_t[a])
plt.scatter(scores_adp[:, 0], scores_adp[:,1], marker='d',c= 'darkgreen',s=40,edgecolors='none')#c=all_t[a])
plt.tick_params(axis='both', which='major', labelsize=7)
plt.tick_params(axis='both', which='minor', labelsize=4)
plt.xlabel(r'PC1 ($\AA$)')
plt.ylabel(r'PC2 ($\AA$)')
#cbar = plt.colorbar()
#cbar.set_label('Time steps')
plt.legend(loc='best',frameon=True,scatterpoints=1)
nameout='pca_traj2_all'
plt.subplots_adjust(bottom=0.15, top=0.95)
plt.savefig(nameout+'.png',format='png',dpi=300)


loadings1 = pca1.components_[0]*pca1.explained_variance_ratio_[0]
loadings2 = pca1.components_[1]*pca1.explained_variance_ratio_[1]
loadings3 = pca1.components_[2]*pca1.explained_variance_ratio_[2]
x = np.array(list(xrange(1,len(loadings1)+1)))

fig = plt.figure(figsize=(4,3))
#line1, = plt.plot(x,loadings1,lw=1.5,color='k',label='PC1')
#line2, = plt.plot(x,loadings2,'b--',lw=1.5,label='PC2')
line3, = plt.plot(x,loadings3,'r+-',lw=1.5,label='PC3')
#leg = plt.legend(handles=[line1, line2,line3])
leg = plt.legend(handles=[line3])
for legobj in leg.legendHandles:
  legobj.set_linewidth(1.5)
plt.legend(loc='best',frameon=True)
plt.xlabel('Variables')
plt.ylabel('Loadings')
nameout='Figure_PC_loadings'
plt.subplots_adjust(left=0.2,bottom=0.15, top=0.95)
plt.savefig(nameout+'.png',format='png',dpi=300)

plt.figure(figsize=(4,3))
symbols=['o','v','^','<','>','8','s','p','*','x']
plt.scatter(loadings1, loadings2, marker='o',c= 'darkorange',s=10,edgecolors='none')#c=all_t[a])
plt.tick_params(axis='both', which='major', labelsize=7)
plt.tick_params(axis='both', which='minor', labelsize=4)
plt.xlabel('Loadings PC1')
plt.ylabel('Loadings PC2')
#cbar = plt.colorbar()
#cbar.set_label('Time steps')
plt.legend(loc='best',frameon=True,scatterpoints=1)
nameout='Scatter_loadings_PC1_2'
plt.subplots_adjust(left=0.2,bottom=0.15, top=0.95)
plt.savefig(nameout+'.png',format='png',dpi=300)

rmsd1 = np.sqrt(np.sum(pca1.components_[0].reshape(T9a.n_atoms,3)*pca1.components_[0].reshape(T9a.n_atoms,3),axis=1)*pca1.explained_variance_ratio_[0])
np.savetxt('pca_rmsd1.txt',rmsd1)
rmsd2 = np.sqrt(np.sum(pca1.components_[1].reshape(T9a.n_atoms,3)*pca1.components_[1].reshape(T9a.n_atoms,3),axis=1)*pca1.explained_variance_ratio_[1])
np.savetxt('pca_rmsd2.txt',rmsd2)
rmsd3 = np.sqrt(np.sum(pca1.components_[2].reshape(T9a.n_atoms,3)*pca1.components_[2].reshape(T9a.n_atoms,3),axis=1)*pca1.explained_variance_ratio_[2])
x = np.array(list(xrange(17,T9a.n_atoms+17)))

fig = plt.figure(figsize=(4,3))
line1, = plt.plot(x,rmsd1,lw=1.5,color='k',label='PC1')
line2, = plt.plot(x,rmsd2,'b--',lw=1.5,label='PC2')
#line3, = plt.plot(x,rmsd3,'r+-',lw=1.5,label='PC3')
leg = plt.legend(handles=[line1, line2])
#leg = plt.legend(handles=[line3])
for legobj in leg.legendHandles:
  legobj.set_linewidth(1.5)
plt.legend(loc='best',frameon=True)
plt.xlabel('Residue Index')
plt.ylabel(r'Weighted RMSD Modes ($\AA$)')
plt.xlim((16,367))
nameout='Figure_RMSD_PC_loadings'
plt.subplots_adjust(left=0.2,bottom=0.15, top=0.95)
plt.savefig(nameout+'.png',format='png',dpi=300)
