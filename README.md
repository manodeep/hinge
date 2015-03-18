/* Author: Manodeep Sinha <manodeep@gmail.com>
   Date: March 18, 2015.
   LICENSE: MIT
*/


Description
======================

This repo contains a set of 3 tools to create a full interaction
network from a cosmological simulation. The code is divided into
3 stages:

HaloParentFinder
-------------------
Creates the progenitor-descendant pairs between a target snapshot
and future snapshots. Refer to the README in the HaloParentFinder
sub-directory for more details. 

OrphanFixer
-------------------
Corrects the situation where the halofinder misses subhalos that
pass too close to the host halo center. Such subhalos appear to have
(artificially) disrupted; and reappear as a progenitor-less subhalo
after a few snapshots.

This step can only be run after ``HaloParentFinder`` has processed
all of the snapshots. OrphanFixer creates the exact same interaction
network as the next step, ``MergerTree``. 

MergerTree
------------------
Creates the full interaction network and outputs various interactions.
Adds stellar masses from various fits in the literature. genplotdata.c
contains majority of the routines to generate data used to make plots.

Author
=====================
HINGE was written by Manodeep Sinha. Please contact the author in
case of any issues manodeep@gmail.com

LICENSE
=====================
HINGE is released under the MIT license. Basically, do what you want
with the code including using it in commercial application (however,
in that case, please send me an email - I would like to know how
you could possibly make a commercial application application out
of this code).

Project URL
=====================

* version control (https://bitbucket.org/manodeep/hinge)

