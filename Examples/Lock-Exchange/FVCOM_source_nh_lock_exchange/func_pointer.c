/*
 * $Id: func_pointer.c,v 1.1.2.3 2008/04/03 16:28:02 dstuebe Exp $
 * $Name: New_Input $
 * $Revision: 1.1.2.3 $
 */

#include <stdio.h>

#define maxfuncs 20

/* Can you pass fortran pointers through c? */

typedef struct {
  void (*func)();
} VOIDFUNC_TABLE;


VOIDFUNC_TABLE voidfuncs[maxfuncs];

static int nvoidfuncs =0;

void register_func_(void (*func)(), int *i, int *stat) {
   *stat = -1;
     if(nvoidfuncs <= maxfuncs)
       {
         /* Convert from fortran indexing */
         voidfuncs[nvoidfuncs++].func = func; /*Insert at nvoidfuncs */
         *i = nvoidfuncs; /*Return value nvoidfuncs+1 */
         *stat =0;
       }
}



void call_func_(int *i, int *stat){
  *stat = -1;
  int j = (*i)-1; /* Convert from fortran indexing */
  if ( j >= 0 && j < nvoidfuncs)
    {
      voidfuncs[j].func();
      *stat=0;
    }
}
