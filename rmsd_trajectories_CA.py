#!/usr/bin/env python
# calculate rmsd plots relative to two references
# written by Jose Flores-Canales, 06/2016
 
import sys
import os
import numpy as np 
import mdtraj as md
from pylab import *

USAGE="""
rmsd_calc.py refA refB traj_dcd rmsd_dat 
"""

def calc_rmsd(refA, refB, traj, L1):
  #traj.superpose(refA,0)
  rmsdA = md.rmsd(traj,refA,0,atom_indices=L1)
  #traj.superpose(refB,0)
  rmsdB = md.rmsd(traj,refB,0,atom_indices=L1)
  return rmsdA, rmsdB

def print_rmsd(nameout,rmsdA,rmsdB, flag):
  out = open(nameout+'.txt','w')
  n = rmsdA.shape[0] 
  for i in xrange(n):
    if flag == 1:
      print >> out, i+1, rmsdA[i]*10., rmsdB[i]*10.
    else:
      print >> out, i+1, rmsdA[i][0]*10., rmsdB[i][0]*10.
  out.close()

def plot_rmsd(nameout,rmsdA, rmsdB, traj):
  figure()
  #title('RMSD')
  scatter(rmsdA*10., rmsdB*10., marker='x', c=traj.time)
  cbar = colorbar()
  cbar.set_label('Time steps')
  xlabel(r'RMSD-C$\alpha$ ATP ($\AA$)')
  xlim(0, 9)
  ylabel(r'RMSD-C$\alpha$ ADP ($\AA$)')
  ylim(0, 9)
  #show()
  axes().set_aspect('equal', 'datalim')
  savefig(nameout+'.png')
  close()
 
  figure()
  t = np.array([i for i in xrange(0,rmsdA.shape[0])])
  #plot(t, rmsdA*10,'r--',t, rmsdB*10, 'b')
  xlabel(r'Time Steps') 
  ylabel(r'RMSD-C$\alpha$') 
  line1, = plot(t, rmsdA*10, label="ATP", linestyle='--', color='red')
  line2, = plot(t, rmsdB*10, label="ADP", color='blue') 
  legend(handles=[line1, line2])
  savefig(nameout+'_time.png')
  close()

def plot_distance(nameout,dA, dB, dC, traj):
  print nameout
  figure()
  t = np.array([i for i in xrange(0,dA.shape[0])])
  xlabel(r'Time Steps') 
  ylabel(r'Distance') 
  line1, = plot(t, dA*10, label="V256-V365", linestyle='--', color='red')
  line2, = plot(t, dB*10, label="L30-M228", color='blue') 
  line3, = plot(t, dC*10, label="D91-K146", color='orange') 
  legend(handles=[line1, line2, line3])
  savefig(nameout+'_time.png')
  close()

def plot_hist(pwd, dA, dB,dC):
  dA *= 10.
  dB *= 10.
  dC *= 10.
  plt.figure(1)
  bins = 100
  xmin = np.min(dA)
  xmax = np.max(dA)
  plt.hist(dA, bins)
  plt.title('Histogram of distance')
  plt.xlabel(r'Distance V256-V365C$ \alpha$ ATP ($\AA$)')
  plt.ylabel(r'Frequency')
  axes = plt.gca()
  axes.set_xlim([xmin, xmax])
  filename = 'hist_rmsd_per_atom_1'
  namepng = filename + '.png'
  #namepng = pwd + '/' + filename + '.png'
  plt.savefig(namepng)
  plt.close()

  plt.figure(2)
  bins = 100
  xmin = np.min(dB)
  xmax = np.max(dB)
  plt.hist(dB, bins)
  plt.title('Histogram of distance')
  plt.xlabel(r'Distance L30-M228C$ \alpha$ ATP ($\AA$)')
  plt.ylabel(r'Frequency')
  axes = plt.gca()
  axes.set_xlim([xmin, xmax])
  filename = 'hist_rmsd_per_atom_2'
  namepng = filename + '.png'
  plt.savefig(namepng)
  plt.close()

  plt.figure(3)
  bins = 100
  xmin = np.min(dC)
  xmax = np.max(dC)
  plt.hist(dC, bins)
  plt.title('Histogram of distance')
  plt.xlabel(r'Distance D91-K146C$ \alpha$ ATP ($\AA$)')
  plt.ylabel(r'Frequency')
  axes = plt.gca()
  axes.set_xlim([xmin, xmax])
  filename = 'hist_rmsd_per_atom_3'
  namepng = filename + '.png'
  plt.savefig(namepng)
  plt.close()

def main():
  if len(sys.argv)<4:
    print USAGE
    sys.exit()
  refA = sys.argv[1]
  refB = sys.argv[2]
  traj_dcd = sys.argv[3]
  nameout = sys.argv[4]
  nameout = nameout.split('/')[0]
  traj = md.load(traj_dcd, top=refA)
  refA_pdb = md.load(refA)
  refB_pdb = md.load(refB)
  #atom_to_keep = [a.index for a in traj.topology.atoms if a.name == 'CA']
  atom_to_keep = [a.index for a in refA_pdb.topology.atoms if a.name == 'CA']
  traj.restrict_atoms(atom_to_keep)
  atom_to_keep = [a.index for a in refA_pdb.topology.atoms if a.name == 'CA']
  refA_pdb.restrict_atoms(atom_to_keep)
  atom_to_keep = [a.index for a in refB_pdb.topology.atoms if a.name == 'CA']
  refB_pdb.restrict_atoms(atom_to_keep)
  #topology = traj.topology
  topology = refA_pdb.topology
  #L1 = topology.select('resid 207 to 214 or resid 343 to 349 and name CA')
  NL_L1 = topology.select('resid 343 to 349 and name CA') # NL
  SW1_L1 = topology.select('resid 207 to 214 and name CA') # SW1
  #L1 = topology.select('name CA and not resid 255 to 264 and not resid 343 to 350') # no missing loop and not NL 
  L1 = topology.select('name CA and not resid 255 to 264') # no missing loop
  rmsdA, rmsdB = calc_rmsd(refA_pdb, refB_pdb, traj,NL_L1)
  filename = nameout + '_NL' 
  print_rmsd(filename,rmsdA, rmsdB,1)
  plot_rmsd(filename, rmsdA, rmsdB,traj)
  rmsdA, rmsdB = calc_rmsd(refA_pdb, refB_pdb, traj,SW1_L1)
  filename = nameout + '_SW1' 
  print_rmsd(filename,rmsdA, rmsdB,1)
  plot_rmsd(filename, rmsdA, rmsdB,traj)
  rmsdA, rmsdB = calc_rmsd(refA_pdb, refB_pdb, traj,L1)
  filename = nameout 
  print_rmsd(filename,rmsdA, rmsdB,1)
  plot_rmsd(filename, rmsdA, rmsdB,traj)
  #1-350 PDB tinker model
  #17-366 PDB xray model
  #pairs1 = topology.select('resid 256 or resid 365 and name CA') # NL
  #pairs2 = topology.select('resid 30 or resid 228 and name CA') # SW1
  pairs1 = topology.select('resid 240 or resid 349 and name CA') # NL
  pairs2 = topology.select('resid 14 or resid 212 and name CA') # SW1
  pairs3 = topology.select('resid 75 or resid 130 and name CA') # SW1
  dist1 = md.compute_distances(traj, np.array([pairs1])) 
  dist2 = md.compute_distances(traj, np.array([pairs2])) 
  dist3 = md.compute_distances(traj, np.array([pairs3])) 
  nameout1 = nameout + '_dist'
  print nameout1,nameout
  print_rmsd(nameout1, dist1, dist2,0)
  plot_distance(nameout1, dist1, dist2, dist3, traj)
  plot_hist(nameout, dist1, dist2, dist3)

  dist1_A = md.compute_distances(refA_pdb, np.array([pairs1])) 
  dist2_A = md.compute_distances(refA_pdb, np.array([pairs2])) 
  dist3_A = md.compute_distances(refA_pdb, np.array([pairs3])) 
  print "Ref. A NL dist:", dist1_A[0][0]*10., "SW1 dist:", dist2_A[0][0]*10., "K146-D91:", dist3_A[0][0]*10.
  dist1_B = md.compute_distances(refB_pdb, np.array([pairs1])) 
  dist2_B = md.compute_distances(refB_pdb, np.array([pairs2])) 
  dist3_B = md.compute_distances(refB_pdb, np.array([pairs3])) 
  print "Ref. B NL dist:", dist1_B[0][0]*10., "SW1 dist:", dist2_B[0][0]*10., "K146-D91:", dist3_B[0][0]*10.

if __name__ == '__main__':
    main()

