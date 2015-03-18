Readme for the Parentfinder code  -- MS 10/10/10 !! 

The code finds descendants for halos in a dark matter
only simulation, where the halos have been identified
by P-Groupfinder. Matching is done based on particle
id -- memory intensive. The computer you run this
on should be able to handle O(N^3) bytes 
(where N^3 is the simulation size, e.g., 512^3). 

These are the steps in the code:


0. Finds the hierarchy level of the subhalos using octants. (hierarchy.c)

1. Match NextFofs to PrevFofs (using all particles, including
all the subhalos) based on the highest Ncommon. (findfofparents in findallparents.c)

2. Match remaining PrevGroups to NextGroups based on the 
symmetrized binding energy rank (ref: Boylan-Kolchin, MS-II paper).
There is also the option of implementing a MAXRANKLOC -- the max.
innermost particles that will be used. (findallparents in findallparents.c)

3. See if FOF+Sub in a major merger have their definitions switched in the next
snapshot. What ends up happening is both the PrevFof and PrevSub point to 
NextFof and some random PrevSub points to NextSub (a dissolution event, instead 
of a natural progression). check_fof_matches in switchfof.c fixes these

4. Massive tidal stripping sometimes causes a subhalo to point to the Fof background
at the next step even though the core of the subhalo still survives. Problem has
already been alleviated quite a bit by symmetrizing the binding energy rank and
using MAXRANKLOC in findallparents(). However, if any such subhalos still remain,
they are fixed in find_progenitor() in findprogenitor.c 

Steps 3 and 4 are only executed when comparing two consecutive snapshots. 

5. Match all remaining PrevGroups to the following snapshot using findallparents(). 
Repeat step till all groups have been found or NextGroup is 3 snapshots away from
PrevGroup. 


IDL code to run (and time) the parentfinding code is in runparentfinder.pro -- dependency
load_groups_header.pro is also there. If P-groupfinder has group velocities in the 
groups files, then the GET_GROUPVEL option in the Makefile, and GET_VEL for load_groups.pro
and load_groups_header.pro, must be enabled. 


