#include <R_ext/Rdynload.h>
#include "binning.h"

static const R_CMethodDef cMethods[] = {
    /* binning.c */
    {".binning", (DL_FUNC) & binning, 7},
{NULL, NULL, 0}
};


void R_init_methylPipe(DllInfo * info)
{
    R_registerRoutines(info, cMethods, NULL, NULL, NULL);
}

void R_unload_methylPipe(DllInfo *info)
{
   /* any clean-up when package unloaded */
    (void) info;
}
