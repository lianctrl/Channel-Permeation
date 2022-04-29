#!/usr/bin/env python3

import os
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import align
import argparse
from collections.abc import Iterable
import warnings

warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning) 


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

# function to verify if the passage happened through
# the position ordered string 1,0,-1
# for up to down permeation events (is it our case?)

def count_perm(labelist):

    pass_count = 0

    for i in range (2,len(labelist)):

        if ( labelist[i-2]==1 and labelist[i-1]==0 and labelist[i]==-1):

            pass_count+=1

    return pass_count

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

parser.add_argument("-com", "--center", dest = "com", \
        required = True, help ="<Required> From this selection will \
        be guessed the center of the channel (e.g. resid 100-200)")

parser.add_argument("-r", "--radius", dest = "radius", \
        type=float, default=5.0, help = "estimate of channel radius, default = 5.0 (Ang)")

parser.add_argument("-st", "--starttime", dest = "startt",\
        type=float, default=0, help = "time starting (ns), default = 0")

parser.add_argument("-et", "--endtime", dest = "endt", \
        type=float, default=-1, help = "time ending (ns), default = last frame")

parser.add_argument("-j", "--stride", dest = "stride", \
        type=int, default=1, help = "stride, default = 1")

parser.add_argument("-dx", "--width", dest = "width", \
        type=float, default=0.5, help = "bin width for the channel axis, default = 0.5 (Ang)")

feature_parser = parser.add_mutually_exclusive_group(required=False)
feature_parser.add_argument('--align', dest='toalign', action='store_true')
feature_parser.add_argument('--no-align', dest='toalign', action='store_false')
parser.set_defaults(toalign=False)

args = parser.parse_args()

sel=args.sel
ref = args.ref

u = mda.Universe(args.pdb, args.traj)

sel_atoms = u.select_atoms(args.sel)

ref_atoms = u.select_atoms(args.ref)

com_atoms = u.select_atoms(args.com)

print(f'\n The selection provided is of {len(sel_atoms)} ions \
and the trajectory contains {u.trajectory.n_frames} frames')

# Computed time is about 90 frames/second per ion

time=len(sel_atoms)*u.trajectory.n_frames/90

print(f'\n Estimated time to finish is {time/60} minutes or {time/3600} hours')

# trajectory alignment if requested,
if args.toalign:

    name=os.path.splitext(args.pdb)[0]+'-aligned.dcd'

    print(f'The trajectory {args.traj} will be aligned to {args.pdb} and saved as {name}')

    algn=mda.Universe(args.pdb)

    alignment = align.AlignTraj(u, algn, filename=name, select=args.ref)

    alignment.run()
    
    del u
    u = mda.Universe(args.pdb, name)

z_up=np.amax(ref_atoms.positions[:,2])
z_lw=np.amin(ref_atoms.positions[:,2])

# this is an huuuuge guess, I understand
# the generalization issue but you have
# to be smarter!!

lim_up = [com_atoms.centroid()[0],\
        com_atoms.centroid()[1],z_up]
lim_lw = [com_atoms.centroid()[0],\
        com_atoms.centroid()[1],z_lw]


#vector used down below in order to follow the numerical approach of
# https://stackoverflow.com/questions/47932955/how-to-check-if-a-3d-point-is-inside-a-cylinder
q  = np.zeros(3)
p1 = np.array(lim_lw)
p2 = np.array(lim_up)

#radius of the CNT measured prev through MDAnalysis
r = args.radius 

#define difference between vectors of the two centers (N.B possible only if these two are numpy arrays)
vec = p2 - p1
#follow the link above
const = r * np.linalg.norm(vec)

# declare an empty list to store the labels
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

pass_count = 0

print('\n Main loop starting, good luck and wait',flush=True)

for n in range (len(sel_atoms)):

    labelist = []
    labelist.append(1)

    for frm in u.trajectory[args.startt:args.endt:args.stride]:

        z_sel = sel_atoms.positions[n,2]

        y_sel = sel_atoms.positions[n,1]

        x_sel = sel_atoms.positions[n,0]

        q = np.array([x_sel,y_sel,z_sel])

        a = np.dot(q - p1, vec)
        b = np.dot(q - p2, vec)
        c = np.linalg.norm(np.cross(q - p1, vec))

        if ( a >= 0 and b >= 0 and c <= const and old_step != 1):

            labelist.append(1)
            old_step = 1

        elif ( a >= 0 and b <= 0 and c <= const ):

            x.append(x_sel)
            y.append(y_sel)
            z.append(z_sel)

            if old_step != 0:
                labelist.append(0)
                old_step = 0

        elif ( a <= 0 and b <= 0 and c <= const and old_step != -1):

            labelist.append(-1)
            old_step = -1


    pass_count+=count_perm(labelist)


x=np.array(x)
y=np.array(y)
z=np.array(z)

delta=args.width
ns=int((z.max()-z.min())/delta)
centers=np.linspace(z.min(),z.max(),ns+1,endpoint=True)

binx,_=binner(z,x,centers)
biny,_=binner(z,y,centers)
binz,_=binner(z,z,centers)

c1=np.asarray(list(map(np.mean,binx)))
c2=np.asarray(list(map(np.mean,biny)))
c3=np.asarray(list(map(np.mean,binz)))

c4=np.asarray(list(map(np.std,binx)))
c5=np.asarray(list(map(np.std,biny)))
c6=np.asarray(list(map(np.std,binz)))

c7=(centers[1:]+centers[:-1])/2

np.savetxt('ions-trajectory.dat',np.c_[c1,c2,c7,c4,c5,c6])

print (f'\n The {args.sel} ions have passed through the \
{args.ref} {pass_count} times')
