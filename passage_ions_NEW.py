#!/usr/bin/env python3

import numpy as np
import MDAnalysis as mda
import argparse

#error if less than 4 arguments are passed

parser = argparse.ArgumentParser(description="Analysis \
        of ions permeation through channels")

parser.add_argument("-s", "--structure", dest="pdb", \
        required=True, help="structure file (e.g. .pdb,.psf or .gro)")

parser.add_argument("-t", "--traj", dest = "traj", \
        required = True, help ="trajectory files, also an array can be passed")

parser.add_argument("-sel", "--selection", dest = "sel",\
        required = True, help ="ions selection for which the analysis must be done (e.g. resname NA)")

parser.add_argument("-ref", "--reference", dest = "ref", \
        required = True, help ="From this selection will be guessed the first \
        and last atoms of the reference channel (e.g. protein and backbone)")

parser.add_argument("-st", "--starttime", dest = "starttime",\
        type=float, default=0, help = "time starting (ns), default = 0")

parser.add_argument("-st", "--endtime", dest = "endtime", \
        type=float, default=-1, help = "time ending (ns), default = last frame")

parser.add_argument("-j", "--stride", dest = "stride", \
        type=int, default=1, help = "stride, default = 1")

parser.add_argument("-dx", "--width", dest = "width", \
        type=float, default=0.05, help = "bin width for the channel axis, default = 0.05")

args = parser.parse_args()

sel=args.sel
ref = args.ref

u = mda.Universe(args.pdb, args.traj)

sel_atoms = u.select_atoms(args.sel)

ref_atoms  = u.select_atoms(args.ref)

border_atoms_up = C_atoms[:14]

border_atoms_lw = C_atoms[-14:]

box_up = 90.0

box_lw = 0.0

# arrays to contain x,y,z coordinates of the center top and bottom of the CNT

center_cnt_up = np.zeros(3)
center_cnt_lw = np.zeros(3)

#vector used down below in order to follow the numerical approach of
# https://stackoverflow.com/questions/47932955/how-to-check-if-a-3d-point-is-inside-a-cylinder
q  = np.zeros(3)
p1 = np.zeros(3)
p2 = np.zeros(3)

# declare an empty list to store the labels
labelist = []

#radius of the CNT measured prev through MDAnalysis
r = 4.745

# labeling every position of the NA
# +1 above the CNT
# 0 inside the CNT
# -1 below the CNT
# since it starts from outside the CNT the first label should be 1
# btw initialized out of range so it can append immediately the first value

old_step = 2

#variable counting the passages of NA ions through channel

pass_count = 0

for n in range (0,len(NA_atoms)):

	#counter variable, should start from 0 but 1 is to have match with the vmd reading trajectory
#	i = 1

	for frm in u.trajectory[:]:

    		#find the center (baricenter definition) of the upper part of CNT
		center_cnt_up[2] = sum(border_Catoms_up.positions[:,2])/len(border_Catoms_up.positions[:,2])
        
		center_cnt_up[1] = sum(border_Catoms_up.positions[:,1])/len(border_Catoms_up.positions[:,1])

		center_cnt_up[0] = sum(border_Catoms_up.positions[:,0])/len(border_Catoms_up.positions[:,0])


    		#find the center (baricenter definition) of the lower part of CNT
		center_cnt_lw[2] = sum(border_Catoms_lw.positions[:,2])/len(border_Catoms_lw.positions[:,2])

		center_cnt_lw[1] = sum(border_Catoms_lw.positions[:,1])/len(border_Catoms_lw.positions[:,1])

		center_cnt_lw[0] = sum(border_Catoms_lw.positions[:,0])/len(border_Catoms_lw.positions[:,0])

    		#set the x,y,z coordinate for the NAR ion
		
		z_coordinate_NA = NA_atoms.positions[n,2]

		y_coordinate_NA = NA_atoms.positions[n,1]

		x_coordinate_NA = NA_atoms.positions[n,0]

		p1 = np.array([center_cnt_lw[0],center_cnt_lw[1],center_cnt_lw[2]])

		p2 = np.array([center_cnt_up[0],center_cnt_up[1],center_cnt_up[2]])
    
		#define difference between vectors of the two centers (N.B possible only if these two are numpy arrays)
		vec = p2 - p1
		#follow the link above
		const = r * np.linalg.norm(vec)

		q = np.array([x_coordinate_NA,y_coordinate_NA,z_coordinate_NA])

		if ( np.dot(q - p1, vec) >= 0 and np.dot(q - p2, vec) >= 0 and np.linalg.norm(np.cross(q - p1, vec)) <= const and old_step != 1):

			labelist.append(1)
		
			old_step = 1

		if ( np.dot(q - p1, vec) >= 0 and np.dot(q - p2, vec) <= 0 and np.linalg.norm(np.cross(q - p1, vec)) <= const and old_step != 0):

			labelist.append(0)
		
			old_step = 0


		if ( np.dot(q - p1, vec) <= 0 and np.dot(q - p2, vec) <= 0 and np.linalg.norm(np.cross(q - p1, vec)) <= const and old_step != -1):

			labelist.append(-1)
		
			old_step = -1


# append an out of range value to separate different atoms passages

	labelist.append(2)


labelarray = np.array(labelist)

# loop to verify if the passage happened through the check of the position ordered string 1,0,-1
# for up to down permeation events (our case)

for i in range (0,len(labelist)-3):
	
	if ( labelarray[i]==1 and labelarray[i+1]==0 and labelarray[i+2]==-1):
		
		pass_count+=1

print ('\n The NA ions have passed through the CNT %i times\n' % ( pass_count ))

