#!/Users/martinmccullagh/anaconda/bin/python
#USAGE : 
import scipy
import sys
import os
import numpy
import math
import MDAnalysis
from MDAnalysis.analysis.align import *

# Subroutines
def get_rot_mat(vec1,vec2):
	v = numpy.cross(vec1,vec2)
	cosTheta = numpy.dot(vec1,vec2)
	vx = numpy.matrix([[0,-v[2],v[1]],[v[2],0,-v[0]],[-v[1],v[0],0]],dtype=float)
	rot = numpy.diag([1.0,1.0,1.0]) + vx + numpy.dot(vx,vx)*(1-cosTheta)/math.pow(linalg.norm(v),2)
	return rot

# read the configuration file and populate the global variables
def ParseConfigFile(cfg_file):
	global inp_pdb_file,ref_pdb_file,out_pdb_file,inp_sel,ref_sel
	f = open(cfg_file)
	atom_sel = 'blank'
	for line in f:
		# first remove comments
		if '#' in line:
			line, comment = line.split('#',1)
		if '=' in line:
			option, value = line.split('=',1)
			option = option.strip()
			value = value.strip()
			print "Option:", option, " Value:", value
			# check value
			if option.lower()=='refpdbfile':
				ref_pdb_file = value
			elif option.lower()=='inppdbfile':
				inp_pdb_file = value
			elif option.lower()=='outpdbfile':
				out_pdb_file = value
			elif option.lower()=='atomsel':
				atom_sel = value
			elif option.lower()=='refatomsel':
				ref_sel = value
			elif option.lower()=='inpatomsel':
				inp_sel = value
			else :
				print "Option:", option, " is not recognized"


	if atom_sel == 'blank':
		print "atom selection not defined"
	else:
		print "Using the same atom selection for both pdbs.  Change this by removing keyword atomsel and adding refatomsel and inpatomsel."
		inp_sel = atom_sel
		ref_sel = atom_sel

# compute the distance between two position vectors taking into account PBC
# Main Program

# read in command line argument
cfg_file = sys.argv[1]

# read cfg file
ParseConfigFile(cfg_file)

print "Initial PDB file:", inp_pdb_file
print "Input atom selection:", inp_sel
print "Reference PDB file:", ref_pdb_file
print "Reference atom selection:", ref_sel

# read in intial coordinates
coord = MDAnalysis.Universe(inp_pdb_file)

# read in reference pdb file (coordinates and number of atoms)
ref = MDAnalysis.Universe(ref_pdb_file)

# make atom selections
inp_sel = coord.selectAtoms(inp_sel)
ref_sel = ref.selectAtoms(ref_sel)

# Check that number of atoms in selection are the same
if inp_sel.numberOfAtoms() != ref_sel.numberOfAtoms():
	sys.exit("Number of atoms in selection are different.  Bombing out...")

print "Number of atoms in selection:", inp_sel.numberOfAtoms()

# subtract the average from the selected coordinates
disp_vec = ref_sel.centerOfMass() 
ref_sel.translate(-disp_vec)
coord.atoms.translate(-inp_sel.centerOfMass())

#compute rotation matrix
R, rmsd = MDAnalysis.analysis.align.rotation_matrix(inp_sel.positions,ref_sel.positions)
#print R, disp_vec

#Apply rotation and translation and then print new pdb
coord.atoms.rotate(R)
coord.atoms.translate(disp_vec)
coord.atoms.write(out_pdb_file)

