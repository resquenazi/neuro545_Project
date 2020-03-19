/*****************************************************************************/
/*                                                                           */
/*  pop_low_util.c                                                           */
/*  wyeth bair                                                               */
/*                                                                           */
/*  Lowest level of 'pop' utilites.                                          */
/*                                                                           */
/*C  (C) Copyright  Wyeth Bair 1993-2010                                     */
/*                                                                           */
/*****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "my_util.h"
#include "carray_util.h"
#include "iarray_util.h"
#include "farray_util.h"
#include "paramfile_util.h"

#include "mod.h" // Data structures
#include "ifc.h" // IFC and POP data structures
#include "paramfile.h"

/**************************************-**************************************/
/*                                                                           */
/*                          POP_LOW_GET_LAYER_FOR_NAME                       */
/*                                                                           */
/*  Return a pointer to the layer with the given name, or NULL.              */
/*                                                                           */
/*****************************************************************************/
struct pop_layer *pop_low_get_layer_for_name(laylist,nlay,name)
     struct pop_layer **laylist;    /* List of layers */
     int nlay;                      /* Number of layers */
     char name[];                   /* Layer name */
{
  int i;
  int done;
  struct pop_layer *pl;
  
  done = 0;
  i = 0;
  pl = NULL;
  while(!done){
    if (i >= nlay)
      done = 1;
    else if (strcmp(laylist[i]->name,name)==0){
      pl = laylist[i];
      done = 1;
    }else
      i += 1;
  }

  return pl;
}
