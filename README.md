# Channel-Permeation 
## Description
This python script helps the user in perform calculations of permeation of ions through biological channels. \
It also determines a trajectory averaged path followed by the ion through the channel.
## Library used
This script needs the following python libs
1. numpy >= 1.21  
2. MDAnalysis >= 2.0.0
3. multiprocessing for the parallelized version of the script
## How to use
It is possible to run the python scripts following the below line command as an example:\
python permeation-analysis.py -s structure-file -t trajectory1 -t trajectory2 ... -sel "selection" -ref "reference" -r channel_radius -dx delta

The **selection** refers to the ions types to which one wants to compute the permeation. \
The **reference** is the atoms of the channel one wants to select (e.g protein or protein and backbone, etc.)\
**selection** and **reference** must be strings in VMD style selection. \
The **delta** is the step to which the channel has been discretized \
It is also possible to select a part of the trajectory with the optional flags:
1. -st start at n frame
2. -et end at m frame 
3. -j skip every k frames 

Now it is also possible to align a trajectory to the reference structure with respect to ref with the **--align** flag
