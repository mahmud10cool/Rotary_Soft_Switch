#include "trial_with_intermediate_tanks_capi_host.h"
static trial_with_intermediate_tanks_host_DataMapInfo_T root;
static int initialized = 0;
rtwCAPI_ModelMappingInfo *getRootMappingInfo()
{
    if (initialized == 0) {
        initialized = 1;
        trial_with_intermediate_tanks_host_InitializeDataMapInfo(&(root), "trial_with_intermediate_tanks");
    }
    return &root.mmi;
}

rtwCAPI_ModelMappingInfo *mexFunction(){return(getRootMappingInfo());}
