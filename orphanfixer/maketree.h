#pragma once

#include <stdio.h>
#include <stdlib.h>

#include "defs.h"
#include "io.h"
#include "loadparents.h"

struct node_data {
  int64 nodeloc;
  int64 haloid;
  float z;
  short snapshot;

  float xcen;
  float ycen;
  float zcen;

  float vxcen;
  float vycen;
  float vzcen;

  float meanvel[3]; /*all velocities are assumed to be "Gadget" snapshot
                       velocities and converted to peculiar velocities by
                       multiplying with sqrt(a) */

#ifdef SUSSING_TREES
  long SUSSING_haloID;
#endif

  double Mtot;
  double BoundFofMtot;

  int64
      ParentID; /* Group number for the parent group at a following snapshot */
  float ParentZ;
  short ParentSnapshot;

  short isFof;
  short ParentLevel; /* misnomer. means the hierarchy level */
  int64 ContainerId; /* Group number that contains this halo. For a Fof, this
                        will be the Fof itself */

#ifndef FOF_ONLY
  int64 Nsub; /* actual number of subhalos. Nsub equals to subhalos for the Fof,
                 and the number of (sub)_n-subs for subhalos */
#endif
  int64 Nchild; /* actual number of pointers that consider this halo to be a
                   parent */

  struct node_data *FofHalo; /* FOF container at this snapshot */
  struct node_data *Sibling;
  struct node_data
      *Parent; /* The halo this is going to go into in a future snapshot.
                  Depends on the parent matching algorithm - mostly at the
                  immediately next snapshot */
  struct node_data *
      BigChild; /* Progenitor at a previous snapshot. Similar to Parent (mostly
                   the immediately previous snapshot but not necessarily so ) */

  struct node_data *ContainerHalo; /* container halo; for a FOF/sub-halo,
                                      container is the Fof; a sub-subhalo will
                                      have a container as a subhalo and so on*/
};

/* functions in maketree */
void maketree(struct parent_data *allparents[], int64 *Ngroups,
              struct node_data *tree[]);
void assign_haloid(struct node_data *tree[], int64 *Ngroups);
void assign_parent(struct node_data *halo, int64 haloid,
                   struct node_data *parenthalo, int64 parentid);
