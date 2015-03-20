/* Author: Manodeep Sinha <manodeep@gmail.com>
   Date: March 18, 2015.
   LICENSE: MIT
*/

# Contents 

* [Description](#description)
* [Installation](#installing-hinge)
* [Running the codes](#running-hinge)


# Description

This repo contains a set of 3 tools (written in C) to create a
full interaction network from a cosmological simulation. The code
is divided into 3 stages:

## HaloParentFinder

Creates the progenitor-descendant pairs between a target snapshot
and future snapshots. Refer to the README in the HaloParentFinder
sub-directory for more details. 

## OrphanFixer

Corrects the situation where the halofinder misses subhalos that
pass too close to the host halo center. Such subhalos appear to have
(artificially) disrupted; and reappear as a progenitor-less subhalo
after a few snapshots.

This step can only be run after ``HaloParentFinder`` has processed
all of the snapshots. OrphanFixer creates the exact same interaction
network as the next step, ``MergerTree``. 

## MergerTree

Creates the full interaction network and outputs various interactions.
Adds stellar masses from various fits in the literature. genplotdata.c
contains majority of the routines to generate data used to make plots.


# Installing HINGE


## Prerequisites

[GNU Scientific Library](http://www.gnu.org/software/gsl/ "GSL") is the 
only pre-requisite for HINGE. And a working C compiler, of course. 


## HINGE Compile Options

Get the source code from the repo, [here](https://bitbucket.org/manodeep/hinge/ "HINGE repository"), ::

    hg clone https://USERNAME@bitbucket.org/manodeep/hinge

Edit the ``common.mk`` file to set up your compilation options, 

* BIGSIM     -- Particle load is > INT_MAX (~ 2 Billion). All loop counters will be 64-bit integers
* LONGIDS    -- Particle IDs are 64-bit integers
* MAKE_LEAN  -- Remove memory allocations when they are no longer required (recommended)

* SUSSING_TREES -- Read in AHF halos generated with the SUSSING TREES compile option
* BGC2 					-- Read in Rockstar halos written in the bgc2 format

* WMAP5 				-- Set up WMAP5 cosmology
* WMAP3 				-- Set up WMAP3 cosmology
* WMAP1 				-- Set up WMAP1 cosmology


## Deprecated HINGE Compile Options

* FOF_ONLY			-- Generate the interaction network just based on FOF halos
* GET_GROUPVEL	-- Read in group velocities for [SubFind](http://enzo-project.org/ "SubFind is bundled with Enzo") halos. 


## Options under development

These are the following compile options that will be active in the near-future, 

* BGC2  		-- Read-in [Rockstar](https://bitbucket.org/gfcstanford/rockstar "Rockstar Repository") halo catalogs that have been written in the ``bgc2`` format
* USE_OMP   -- Enable OpenMP parallelization

# Running HINGE

## Input Formats for Halo Catalogs

HINGE currently supports standard SubFind halos, and [AHF](http://popia.ft.uam.es/AHF/ "Download AHF") halos generated
with the ``SUSSING2013`` Makefile option. 


## Common parameters

Here are the parameters that are common to all three codes, 

Parameter name                                     |  Parameter Effect
[----------------------------------------------------|-------------------------------
 MIN_SNAPSHOT_NUM			 | Minimum snapshot number (integer).
 MAX_SNAPSHOT_NUM			| Maximum snapshot number (integer).
 SNAPSHOT_DIR					| Directory where all the snapshots are stored (string).
 SNAPSHOT_BASE				| Basename for all the snapshots (string). Fully qualified snapshot names are generated using the C printf format ``sprintf(%s/%s_%03d,SNAPSHOT_DIR,SNAPSHOT_BASE,snapshot_number)``.
 GROUP_DIR						| Directory where the halo catalogs are stored (string)
 GROUP_BASE					| Basename for all the halo catalogs (string). For instance, 			SubFind writes out halo catalogs as groups_XXX.*. So, if you are using SubFind halos then, ``GROUP_BASE`` should be `groups`.
 OUTPUT_DIR					| Directory where all output is written


## Parameter File for HaloParentFinder

Please refer to the [HaloParentFinder/README](haloparentfinder/README.md) file. 

## Parameter File for OrphanFixer

Please refer to the [OrphanFixer/README](orphanfixer/README.md) file. 

## Parameter File for MergerTree

Please refer to the [MergerTree/README](mergertree/README.md) file.

# Author


HINGE was written by Manodeep Sinha. Please contact the [author](mailto:manodeep@gmail.com) in
case of any issues.

# LICENSE

HINGE is released under the MIT license. Basically, do what you want
with the code including using it in commercial application (however,
in that case, please send me an email - I would like to know how
you could possibly make a commercial application application out
of this code).

# Project URL
 
* version control (https://bitbucket.org/manodeep/hinge)