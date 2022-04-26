#!/usr/bin/env python3

import numpy as np
import MDAnalysis as mda
import argparse

def binner(coords,tobinarray,bins=100):

    if isinstance(bins, Iterable):
        bins = np.asarray(bins)
        bins = bins[np.argsort(bins)]
    else:
        range = (coords.min(), coords.max())
        xmin, xmax = range
        if xmin == xmax:
            xmin -= 0.5
            xmax += 0.5
        bins = np.linspace(xmin, xmax, bins+1, endpoint=True)

    idx = np.argsort(coords)
    coords = coords[idx]
    var = tobinarray[idx]
    # left: inserts at i where coords[:i] < edge
    # right: inserts at i where coords[:i] <= edge
    # r_ concatenates

    bin_idx = np.r_[coords.searchsorted(bins, side='right')]
    binned = [var[i:j] for i, j in zip(bin_idx[:-1], bin_idx[1:])]

    return binned, bins


#error if less than 4 arguments are passed

parser = argparse.ArgumentParser(description="Analysis \
        of ions permeation through channels")

parser.add_argument("-s", "--structure", dest="pdb", \
        required=True, help="<Required> structure file (e.g. .pdb,.psf or .gro)")

parser.add_argument("-t", "--traj", dest = "traj", \
        action='append', required = True, help ="<Required> traj files\
        , an array can be passed through many invocation of -t")

parser.add_argument("-sel", "--selection", dest = "sel",\
        required = True, help ="<Required> ions selection for which \
        the analysis must be done (e.g. resname NA)")

parser.add_argument("-ref", "--reference", dest = "ref", \
        required = True, help ="<Required> From this selection will \
        be guessed the first and last atoms of the reference \
        channel (e.g. protein and backbone)")

parser.add_argument("-st", "--starttime", dest = "startt",\
        type=float, default=0, help = "time starting (ns), default = 0")

parser.add_argument("-et", "--endtime", dest = "endt", \
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

# insert here a RMSD fit of the traj on the structure,
# it should solve the alignment problem!!

z_up=np.amax(ref_atoms.positions[:,2])
z_lw=np.amin(ref_atoms.positions[:,2])

# this is an huuuuge guess, I understand
# the generalization issue but you have
# to be smarter!!

lim_up = [ref_atoms.center_of_mass()[0],\
        ref_atoms.center_of_mass()[1],z_up]
lim_lw = [ref_atoms.center_of_mass()[0],\
        ref_atoms.center_of_mass()[1],z_lw]


#vector used down below in order to follow the numerical approach of
# https://stackoverflow.com/questions/47932955/how-to-check-if-a-3d-point-is-inside-a-cylinder
q  = np.zeros(3)
p1 = np.array(lim_lw)
p2 = np.array(lim_up)

#radius of the CNT measured prev through MDAnalysis
r = 5.00 

#define difference between vectors of the two centers (N.B possible only if these two are numpy arrays)
vec = p2 - p1
#follow the link above
const = r * np.linalg.norm(vec)

# declare an empty list to store the labels
labelist = []
x=[]
y=[]
z=[]

# labeling every position of the NA
# +1 above the channel
# 0 inside the channel
# -1 below the channel
# since it starts from outside the CNT the first label should be 1
# btw initialized out of range so it can append immediately the first value

old_step = 2

# here we assume the traj centered (user selection frozen)
# again above do the RMSD fit!!


for n in range (len(sel_atoms)):
    for frm in u.trajectory[args.startt:args.endt:args.stride]:

        z_sel = sel_atoms.positions[n,2]

        y_sel = sel_atoms.positions[n,1]

        x_sel = sel_atoms.positions[n,0]

        q = np.array([x_sel,y_sel,z_sel])

        if ( np.dot(q - p1, vec) >= 0 and np.dot(q - p2, vec) >= 0 and np.linalg.norm(np.cross(q - p1, vec)) <= const and old_step != 1):

            labelist.append(1)
            old_step = 1

        elif ( np.dot(q - p1, vec) >= 0 and np.dot(q - p2, vec) <= 0 and np.linalg.norm(np.cross(q - p1, vec)) <= const and old_step != 0):

           labelist.append(0)
           old_step = 0
           x.append(x_sel)
           y.append(y_sel)
           z.append(z_sel)

        elif ( np.dot(q - p1, vec) <= 0 and np.dot(q - p2, vec) <= 0 and np.linalg.norm(np.cross(q - p1, vec)) <= const and old_step != -1):

           labelist.append(-1)
           old_step = -1


# append an out of range value to separate different atoms passages

    labelist.append(2)


labelist = np.array(labelist)
x=np.array(x)
y=np.array(y)
z=np.array(z)

delta=2
centers=np.arange(z.min(),z.max()+delta,delta)
binx=binner(z,x,centers)
biny=binner(z,y,centers)
binz=binner(z,z,centers)

c1=list(map(np.mean,binx))
c2=list(map(np.mean,biny))
c3=list(map(np.mean,binz))

c4=list(map(np.std,binx))
c5=list(map(np.std,biny))
c6=list(map(np.std,binz))

np.savetxt('ions-trajectory.dat',np.c_[c1,c2,c3,c4,c5,c6])

# loop to verify if the passage happened through
# the position ordered string 1,0,-1
# for up to down permeation events (is it our case?)

pass_count = 0

if labelist[0]==0 and labelist[1]==-1:
    pass_count+=1

if labelist[1]==0 and labelist[2]==-1:
    pass_count+=1

for i in range (2,len(labelist)):

    if ( labelist[i-2]==1 and labelist[i-1]==0 and labelist[i]==-1):

        pass_count+=1

print (f'\n The {args.sel} ions have passed through the \
        {args.ref} {pass_count} times')
