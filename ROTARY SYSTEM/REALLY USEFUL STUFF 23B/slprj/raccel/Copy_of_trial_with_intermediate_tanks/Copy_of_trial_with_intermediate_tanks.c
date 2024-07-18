#include "Copy_of_trial_with_intermediate_tanks.h"
#include "rtwtypes.h"
#include "mwmathutil.h"
#include "Copy_of_trial_with_intermediate_tanks_private.h"
#include "rt_logging_mmi.h"
#include "Copy_of_trial_with_intermediate_tanks_capi.h"
#include "zero_crossing_types.h"
#include "Copy_of_trial_with_intermediate_tanks_dt.h"
extern void * CreateDiagnosticAsVoidPtr_wrapper ( const char * id , int nargs
, ... ) ; extern ssExecutionInfo gblExecutionInfo ; RTWExtModeInfo *
gblRTWExtModeInfo = NULL ; void raccelForceExtModeShutdown ( boolean_T
extModeStartPktReceived ) { if ( ! extModeStartPktReceived ) { boolean_T
stopRequested = false ; rtExtModeWaitForStartPkt ( gblRTWExtModeInfo , 3 , &
stopRequested ) ; } rtExtModeShutdown ( 3 ) ; }
#include "slsv_diagnostic_codegen_c_api.h"
#include "slsa_sim_engine.h"
#ifdef RSIM_WITH_SOLVER_MULTITASKING
boolean_T gbl_raccel_isMultitasking = 1 ;
#else
boolean_T gbl_raccel_isMultitasking = 0 ;
#endif
boolean_T gbl_raccel_tid01eq = 0 ; int_T gbl_raccel_NumST = 4 ; const char_T
* gbl_raccel_Version = "23.2 (R2023b) 01-Aug-2023" ; void
raccel_setup_MMIStateLog ( SimStruct * S ) {
#ifdef UseMMIDataLogging
rt_FillStateSigInfoFromMMI ( ssGetRTWLogInfo ( S ) , & ssGetErrorStatus ( S )
) ;
#else
UNUSED_PARAMETER ( S ) ;
#endif
} static DataMapInfo rt_dataMapInfo ; DataMapInfo * rt_dataMapInfoPtr = &
rt_dataMapInfo ; rtwCAPI_ModelMappingInfo * rt_modelMapInfoPtr = & (
rt_dataMapInfo . mmi ) ; int_T enableFcnCallFlag [ ] = { 1 , 1 , 1 , 1 } ;
const char * raccelLoadInputsAndAperiodicHitTimes ( SimStruct * S , const
char * inportFileName , int * matFileFormat ) { return
rt_RAccelReadInportsMatFile ( S , inportFileName , matFileFormat ) ; }
#include "simstruc.h"
#include "fixedpoint.h"
#include "slsa_sim_engine.h"
#include "simtarget/slSimTgtSLExecSimBridge.h"
#define bjcsxykqxt (-1)
B rtB ; X rtX ; DW rtDW ; PrevZCX rtPrevZCX ; static SimStruct model_S ;
SimStruct * const rtS = & model_S ; void MdlInitialize ( void ) { boolean_T
tmp ; rtX . eer5ulgpq1 = rtP . Integrator4_IC ; rtDW . aku4xs4fdc = 1 ; if (
ssIsFirstInitCond ( rtS ) ) { rtX . gvqywnxquo = 1.0E+7 ; tmp =
slIsRapidAcceleratorSimulating ( ) ; if ( tmp ) { tmp =
ssGetGlobalInitialStatesAvailable ( rtS ) ; rtDW . aku4xs4fdc = ! tmp ; }
else { rtDW . aku4xs4fdc = 1 ; } rtX . kllv1nez4n = 2.0E+7 ; } rtX .
mui3hk2kdn = rtP . Integrator13_IC ; rtDW . aigzt14viy = 1 ; if (
ssIsFirstInitCond ( rtS ) ) { tmp = slIsRapidAcceleratorSimulating ( ) ; if (
tmp ) { tmp = ssGetGlobalInitialStatesAvailable ( rtS ) ; rtDW . aigzt14viy =
! tmp ; } else { rtDW . aigzt14viy = 1 ; } rtX . orxgnzk15p = 1.0E+7 ; } rtDW
. jwuiy2atz3 = 1 ; if ( ssIsFirstInitCond ( rtS ) ) { tmp =
slIsRapidAcceleratorSimulating ( ) ; if ( tmp ) { tmp =
ssGetGlobalInitialStatesAvailable ( rtS ) ; rtDW . jwuiy2atz3 = ! tmp ; }
else { rtDW . jwuiy2atz3 = 1 ; } rtX . iikwtjyqxf = 200.0 ; } rtX .
apttuwn22n = rtP . Integrator14_IC ; rtX . avebadcay4 = rtP . Integrator7_IC
; rtX . o2kfetkxay = rtP . Integrator12_IC ; rtX . fvfvpkek2n = rtP .
Integrator10_IC ; rtX . aqqzk12fdv = rtP . Integrator8_IC ; rtX . e30fqjkckb
= rtP . Integrator9_IC ; rtX . mizgi5ztop = rtP . Integrator11_IC ; rtDW .
isc2qgvxqr = 1 ; if ( ssIsFirstInitCond ( rtS ) ) { tmp =
slIsRapidAcceleratorSimulating ( ) ; if ( tmp ) { tmp =
ssGetGlobalInitialStatesAvailable ( rtS ) ; rtDW . isc2qgvxqr = ! tmp ; }
else { rtDW . isc2qgvxqr = 1 ; } } rtX . c5bbd2cwaf = rtP . Integrator6_IC ;
rtX . loqzkyynb5 = rtP . Integrator16_IC ; rtX . fklbo31ygf = rtP .
Integrator5_IC ; rtDW . j0cswhuhs0 = rtP . Delay_InitialCondition ; rtX .
kwoozlzy2w = rtP . Integrator15_IC ; rtX . klye04tpmo [ 0 ] = 0.0 ; rtX .
ht5behjiqe [ 0 ] = 0.0 ; rtX . aa2dkshjtd [ 0 ] = 0.0 ; rtX . kco1vq4cxg [ 0
] = 0.0 ; rtX . mvzjw32hje [ 0 ] = 0.0 ; rtX . klye04tpmo [ 1 ] = 0.0 ; rtX .
ht5behjiqe [ 1 ] = 0.0 ; rtX . aa2dkshjtd [ 1 ] = 0.0 ; rtX . kco1vq4cxg [ 1
] = 0.0 ; rtX . mvzjw32hje [ 1 ] = 0.0 ; rtDW . jz3v2t25p2 = false ; rtDW .
loaiy2cacx = bjcsxykqxt ; rtDW . fxawe3htlm = false ; rtDW . kvbq5jwbvw =
bjcsxykqxt ; rtDW . eyyp3zs3le = false ; rtDW . apkag20puf = bjcsxykqxt ;
rtDW . eiohhcwuqa = false ; rtDW . glwrm0q1go = bjcsxykqxt ; rtDW .
ffd0olavlx = false ; rtDW . inlqdh4kh1 = bjcsxykqxt ; rtDW . hk25kkxqwm =
false ; rtDW . atxpvhc5n0 = bjcsxykqxt ; rtDW . nas3vtbw5t = false ; rtDW .
je5l4cyavv = bjcsxykqxt ; rtB . oguodtbfqn = rtP . Out1_Y0 ; rtDW .
awqg33it32 = false ; rtDW . ebwtwofxpd = bjcsxykqxt ; rtDW . dfvcjzwy4h =
false ; rtDW . hamy5qs5qe = bjcsxykqxt ; rtDW . m0i52xetkk = false ; rtDW .
e0qqpzfruh = bjcsxykqxt ; rtDW . irq32iym0t = false ; rtDW . audz3xgmq3 =
bjcsxykqxt ; rtDW . jkvzly4a1o = false ; rtDW . epcaebihmu = bjcsxykqxt ;
rtDW . ct4cduudvv = false ; rtDW . d34s0zbduf = bjcsxykqxt ; rtDW .
ljfwvqlys1 = false ; rtDW . p1ctkshy11 = bjcsxykqxt ; rtDW . n4ne2vzdw3 =
false ; rtDW . a4pqll03x1 = bjcsxykqxt ; rtDW . nakeaw2rna = false ; rtDW .
mpaqmqlwie = bjcsxykqxt ; } void MdlStart ( void ) { { bool
externalInputIsInDatasetFormat = false ; void * pISigstreamManager =
rt_GetISigstreamManager ( rtS ) ;
rtwISigstreamManagerGetInputIsInDatasetFormat ( pISigstreamManager , &
externalInputIsInDatasetFormat ) ; if ( externalInputIsInDatasetFormat ) { }
} { { { bool isStreamoutAlreadyRegistered = false ; { sdiSignalSourceInfoU
srcInfo ; sdiLabelU loggedName = sdiGetLabelFromChars ( "P3" ) ; sdiLabelU
origSigName = sdiGetLabelFromChars ( "P3" ) ; sdiLabelU propName =
sdiGetLabelFromChars ( "P3" ) ; sdiLabelU blockPath = sdiGetLabelFromChars (
"Copy_of_trial_with_intermediate_tanks/To Workspace" ) ; sdiLabelU blockSID =
sdiGetLabelFromChars ( "" ) ; sdiLabelU subPath = sdiGetLabelFromChars ( "" )
; sdiDims sigDims ; sdiLabelU sigName = sdiGetLabelFromChars ( "P3" ) ;
sdiAsyncRepoDataTypeHandle hDT = sdiAsyncRepoGetBuiltInDataTypeHandle (
DATA_TYPE_DOUBLE ) ; { sdiComplexity sigComplexity = REAL ;
sdiSampleTimeContinuity stCont = SAMPLE_TIME_CONTINUOUS ; int_T sigDimsArray
[ 1 ] = { 1 } ; sigDims . nDims = 1 ; sigDims . dimensions = sigDimsArray ;
srcInfo . numBlockPathElems = 1 ; srcInfo . fullBlockPath = ( sdiFullBlkPathU
) & blockPath ; srcInfo . SID = ( sdiSignalIDU ) & blockSID ; srcInfo .
subPath = subPath ; srcInfo . portIndex = 0 + 1 ; srcInfo . signalName =
sigName ; srcInfo . sigSourceUUID = 0 ; rtDW . homjuevopz . AQHandles =
sdiStartAsyncioQueueCreation ( hDT , & srcInfo , rt_dataMapInfo . mmi .
InstanceMap . fullPath , "679d4c65-e4b9-4c1e-a9cc-57cd620e7bec" ,
sigComplexity , & sigDims , DIMENSIONS_MODE_FIXED , stCont , "" ) ;
sdiCompleteAsyncioQueueCreation ( rtDW . homjuevopz . AQHandles , hDT , &
srcInfo ) ; if ( rtDW . homjuevopz . AQHandles ) {
sdiSetSignalSampleTimeString ( rtDW . homjuevopz . AQHandles , "Continuous" ,
0.0 , ssGetTFinal ( rtS ) ) ; sdiSetSignalRefRate ( rtDW . homjuevopz .
AQHandles , 0.0 ) ; sdiSetRunStartTime ( rtDW . homjuevopz . AQHandles ,
ssGetTaskTime ( rtS , 0 ) ) ; sdiAsyncRepoSetSignalExportSettings ( rtDW .
homjuevopz . AQHandles , 1 , 0 ) ; sdiAsyncRepoSetSignalExportName ( rtDW .
homjuevopz . AQHandles , loggedName , origSigName , propName ) ;
sdiAsyncRepoSetBlockPathDomain ( rtDW . homjuevopz . AQHandles ) ; }
sdiFreeLabel ( sigName ) ; sdiFreeLabel ( loggedName ) ; sdiFreeLabel (
origSigName ) ; sdiFreeLabel ( propName ) ; sdiFreeLabel ( blockPath ) ;
sdiFreeLabel ( blockSID ) ; sdiFreeLabel ( subPath ) ; } } if ( !
isStreamoutAlreadyRegistered ) { { sdiLabelU varName = sdiGetLabelFromChars (
"P1" ) ; sdiRegisterWksVariable ( rtDW . homjuevopz . AQHandles , varName ,
"timeseries" ) ; sdiFreeLabel ( varName ) ; } } } } } { { { bool
isStreamoutAlreadyRegistered = false ; { sdiSignalSourceInfoU srcInfo ;
sdiLabelU loggedName = sdiGetLabelFromChars ( "Integrator7" ) ; sdiLabelU
origSigName = sdiGetLabelFromChars ( "" ) ; sdiLabelU propName =
sdiGetLabelFromChars ( "Integrator7" ) ; sdiLabelU blockPath =
sdiGetLabelFromChars ( "Copy_of_trial_with_intermediate_tanks/To Workspace1"
) ; sdiLabelU blockSID = sdiGetLabelFromChars ( "" ) ; sdiLabelU subPath =
sdiGetLabelFromChars ( "" ) ; sdiDims sigDims ; sdiLabelU sigName =
sdiGetLabelFromChars ( "Integrator7" ) ; sdiAsyncRepoDataTypeHandle hDT =
sdiAsyncRepoGetBuiltInDataTypeHandle ( DATA_TYPE_DOUBLE ) ; { sdiComplexity
sigComplexity = REAL ; sdiSampleTimeContinuity stCont =
SAMPLE_TIME_CONTINUOUS ; int_T sigDimsArray [ 1 ] = { 1 } ; sigDims . nDims =
1 ; sigDims . dimensions = sigDimsArray ; srcInfo . numBlockPathElems = 1 ;
srcInfo . fullBlockPath = ( sdiFullBlkPathU ) & blockPath ; srcInfo . SID = (
sdiSignalIDU ) & blockSID ; srcInfo . subPath = subPath ; srcInfo . portIndex
= 0 + 1 ; srcInfo . signalName = sigName ; srcInfo . sigSourceUUID = 0 ; rtDW
. h1nrd1cc1r . AQHandles = sdiStartAsyncioQueueCreation ( hDT , & srcInfo ,
rt_dataMapInfo . mmi . InstanceMap . fullPath ,
"afa1a4dd-40b7-42cf-9a73-58da46a7060a" , sigComplexity , & sigDims ,
DIMENSIONS_MODE_FIXED , stCont , "" ) ; sdiCompleteAsyncioQueueCreation (
rtDW . h1nrd1cc1r . AQHandles , hDT , & srcInfo ) ; if ( rtDW . h1nrd1cc1r .
AQHandles ) { sdiSetSignalSampleTimeString ( rtDW . h1nrd1cc1r . AQHandles ,
"Continuous" , 0.0 , ssGetTFinal ( rtS ) ) ; sdiSetSignalRefRate ( rtDW .
h1nrd1cc1r . AQHandles , 0.0 ) ; sdiSetRunStartTime ( rtDW . h1nrd1cc1r .
AQHandles , ssGetTaskTime ( rtS , 0 ) ) ; sdiAsyncRepoSetSignalExportSettings
( rtDW . h1nrd1cc1r . AQHandles , 1 , 0 ) ; sdiAsyncRepoSetSignalExportName (
rtDW . h1nrd1cc1r . AQHandles , loggedName , origSigName , propName ) ;
sdiAsyncRepoSetBlockPathDomain ( rtDW . h1nrd1cc1r . AQHandles ) ; }
sdiFreeLabel ( sigName ) ; sdiFreeLabel ( loggedName ) ; sdiFreeLabel (
origSigName ) ; sdiFreeLabel ( propName ) ; sdiFreeLabel ( blockPath ) ;
sdiFreeLabel ( blockSID ) ; sdiFreeLabel ( subPath ) ; } } if ( !
isStreamoutAlreadyRegistered ) { { sdiLabelU varName = sdiGetLabelFromChars (
"Regen" ) ; sdiRegisterWksVariable ( rtDW . h1nrd1cc1r . AQHandles , varName
, "timeseries" ) ; sdiFreeLabel ( varName ) ; } } } } } { { { bool
isStreamoutAlreadyRegistered = false ; { sdiSignalSourceInfoU srcInfo ;
sdiLabelU loggedName = sdiGetLabelFromChars ( "Electric Torque Controller" )
; sdiLabelU origSigName = sdiGetLabelFromChars ( "" ) ; sdiLabelU propName =
sdiGetLabelFromChars ( "Electric Torque Controller" ) ; sdiLabelU blockPath =
sdiGetLabelFromChars ( "Copy_of_trial_with_intermediate_tanks/To Workspace10"
) ; sdiLabelU blockSID = sdiGetLabelFromChars ( "" ) ; sdiLabelU subPath =
sdiGetLabelFromChars ( "" ) ; sdiDims sigDims ; sdiLabelU sigName =
sdiGetLabelFromChars ( "Electric Torque Controller" ) ;
sdiAsyncRepoDataTypeHandle hDT = sdiAsyncRepoGetBuiltInDataTypeHandle (
DATA_TYPE_DOUBLE ) ; { sdiComplexity sigComplexity = REAL ;
sdiSampleTimeContinuity stCont = SAMPLE_TIME_CONTINUOUS ; int_T sigDimsArray
[ 1 ] = { 1 } ; sigDims . nDims = 1 ; sigDims . dimensions = sigDimsArray ;
srcInfo . numBlockPathElems = 1 ; srcInfo . fullBlockPath = ( sdiFullBlkPathU
) & blockPath ; srcInfo . SID = ( sdiSignalIDU ) & blockSID ; srcInfo .
subPath = subPath ; srcInfo . portIndex = 0 + 1 ; srcInfo . signalName =
sigName ; srcInfo . sigSourceUUID = 0 ; rtDW . l444m1qlra . AQHandles =
sdiStartAsyncioQueueCreation ( hDT , & srcInfo , rt_dataMapInfo . mmi .
InstanceMap . fullPath , "eb0a971b-4404-4ef0-a6e2-34f2840dface" ,
sigComplexity , & sigDims , DIMENSIONS_MODE_FIXED , stCont , "" ) ;
sdiCompleteAsyncioQueueCreation ( rtDW . l444m1qlra . AQHandles , hDT , &
srcInfo ) ; if ( rtDW . l444m1qlra . AQHandles ) {
sdiSetSignalSampleTimeString ( rtDW . l444m1qlra . AQHandles , "Continuous" ,
0.0 , ssGetTFinal ( rtS ) ) ; sdiSetSignalRefRate ( rtDW . l444m1qlra .
AQHandles , 0.0 ) ; sdiSetRunStartTime ( rtDW . l444m1qlra . AQHandles ,
ssGetTaskTime ( rtS , 0 ) ) ; sdiAsyncRepoSetSignalExportSettings ( rtDW .
l444m1qlra . AQHandles , 1 , 0 ) ; sdiAsyncRepoSetSignalExportName ( rtDW .
l444m1qlra . AQHandles , loggedName , origSigName , propName ) ;
sdiAsyncRepoSetBlockPathDomain ( rtDW . l444m1qlra . AQHandles ) ; }
sdiFreeLabel ( sigName ) ; sdiFreeLabel ( loggedName ) ; sdiFreeLabel (
origSigName ) ; sdiFreeLabel ( propName ) ; sdiFreeLabel ( blockPath ) ;
sdiFreeLabel ( blockSID ) ; sdiFreeLabel ( subPath ) ; } } if ( !
isStreamoutAlreadyRegistered ) { { sdiLabelU varName = sdiGetLabelFromChars (
"T_elec" ) ; sdiRegisterWksVariable ( rtDW . l444m1qlra . AQHandles , varName
, "timeseries" ) ; sdiFreeLabel ( varName ) ; } } } } } { { { bool
isStreamoutAlreadyRegistered = false ; { sdiSignalSourceInfoU srcInfo ;
sdiLabelU loggedName = sdiGetLabelFromChars ( "Integrator6" ) ; sdiLabelU
origSigName = sdiGetLabelFromChars ( "" ) ; sdiLabelU propName =
sdiGetLabelFromChars ( "Integrator6" ) ; sdiLabelU blockPath =
sdiGetLabelFromChars ( "Copy_of_trial_with_intermediate_tanks/To Workspace11"
) ; sdiLabelU blockSID = sdiGetLabelFromChars ( "" ) ; sdiLabelU subPath =
sdiGetLabelFromChars ( "" ) ; sdiDims sigDims ; sdiLabelU sigName =
sdiGetLabelFromChars ( "Integrator6" ) ; sdiAsyncRepoDataTypeHandle hDT =
sdiAsyncRepoGetBuiltInDataTypeHandle ( DATA_TYPE_DOUBLE ) ; { sdiComplexity
sigComplexity = REAL ; sdiSampleTimeContinuity stCont =
SAMPLE_TIME_CONTINUOUS ; int_T sigDimsArray [ 1 ] = { 1 } ; sigDims . nDims =
1 ; sigDims . dimensions = sigDimsArray ; srcInfo . numBlockPathElems = 1 ;
srcInfo . fullBlockPath = ( sdiFullBlkPathU ) & blockPath ; srcInfo . SID = (
sdiSignalIDU ) & blockSID ; srcInfo . subPath = subPath ; srcInfo . portIndex
= 0 + 1 ; srcInfo . signalName = sigName ; srcInfo . sigSourceUUID = 0 ; rtDW
. fcjl52ojvi . AQHandles = sdiStartAsyncioQueueCreation ( hDT , & srcInfo ,
rt_dataMapInfo . mmi . InstanceMap . fullPath ,
"d093f478-f174-4128-9964-3b45c635e73d" , sigComplexity , & sigDims ,
DIMENSIONS_MODE_FIXED , stCont , "" ) ; sdiCompleteAsyncioQueueCreation (
rtDW . fcjl52ojvi . AQHandles , hDT , & srcInfo ) ; if ( rtDW . fcjl52ojvi .
AQHandles ) { sdiSetSignalSampleTimeString ( rtDW . fcjl52ojvi . AQHandles ,
"Continuous" , 0.0 , ssGetTFinal ( rtS ) ) ; sdiSetSignalRefRate ( rtDW .
fcjl52ojvi . AQHandles , 0.0 ) ; sdiSetRunStartTime ( rtDW . fcjl52ojvi .
AQHandles , ssGetTaskTime ( rtS , 0 ) ) ; sdiAsyncRepoSetSignalExportSettings
( rtDW . fcjl52ojvi . AQHandles , 1 , 0 ) ; sdiAsyncRepoSetSignalExportName (
rtDW . fcjl52ojvi . AQHandles , loggedName , origSigName , propName ) ;
sdiAsyncRepoSetBlockPathDomain ( rtDW . fcjl52ojvi . AQHandles ) ; }
sdiFreeLabel ( sigName ) ; sdiFreeLabel ( loggedName ) ; sdiFreeLabel (
origSigName ) ; sdiFreeLabel ( propName ) ; sdiFreeLabel ( blockPath ) ;
sdiFreeLabel ( blockSID ) ; sdiFreeLabel ( subPath ) ; } } if ( !
isStreamoutAlreadyRegistered ) { { sdiLabelU varName = sdiGetLabelFromChars (
"Energy_In" ) ; sdiRegisterWksVariable ( rtDW . fcjl52ojvi . AQHandles ,
varName , "timeseries" ) ; sdiFreeLabel ( varName ) ; } } } } } { { { bool
isStreamoutAlreadyRegistered = false ; { sdiSignalSourceInfoU srcInfo ;
sdiLabelU loggedName = sdiGetLabelFromChars ( "Gain3" ) ; sdiLabelU
origSigName = sdiGetLabelFromChars ( "" ) ; sdiLabelU propName =
sdiGetLabelFromChars ( "Gain3" ) ; sdiLabelU blockPath = sdiGetLabelFromChars
( "Copy_of_trial_with_intermediate_tanks/To Workspace2" ) ; sdiLabelU
blockSID = sdiGetLabelFromChars ( "" ) ; sdiLabelU subPath =
sdiGetLabelFromChars ( "" ) ; sdiDims sigDims ; sdiLabelU sigName =
sdiGetLabelFromChars ( "Gain3" ) ; sdiAsyncRepoDataTypeHandle hDT =
sdiAsyncRepoGetBuiltInDataTypeHandle ( DATA_TYPE_DOUBLE ) ; { sdiComplexity
sigComplexity = REAL ; sdiSampleTimeContinuity stCont =
SAMPLE_TIME_CONTINUOUS ; int_T sigDimsArray [ 1 ] = { 1 } ; sigDims . nDims =
1 ; sigDims . dimensions = sigDimsArray ; srcInfo . numBlockPathElems = 1 ;
srcInfo . fullBlockPath = ( sdiFullBlkPathU ) & blockPath ; srcInfo . SID = (
sdiSignalIDU ) & blockSID ; srcInfo . subPath = subPath ; srcInfo . portIndex
= 0 + 1 ; srcInfo . signalName = sigName ; srcInfo . sigSourceUUID = 0 ; rtDW
. cpiqtbmrnh . AQHandles = sdiStartAsyncioQueueCreation ( hDT , & srcInfo ,
rt_dataMapInfo . mmi . InstanceMap . fullPath ,
"aa319ba5-14b3-4b41-8a37-bde34ecc0dd5" , sigComplexity , & sigDims ,
DIMENSIONS_MODE_FIXED , stCont , "" ) ; sdiCompleteAsyncioQueueCreation (
rtDW . cpiqtbmrnh . AQHandles , hDT , & srcInfo ) ; if ( rtDW . cpiqtbmrnh .
AQHandles ) { sdiSetSignalSampleTimeString ( rtDW . cpiqtbmrnh . AQHandles ,
"Continuous" , 0.0 , ssGetTFinal ( rtS ) ) ; sdiSetSignalRefRate ( rtDW .
cpiqtbmrnh . AQHandles , 0.0 ) ; sdiSetRunStartTime ( rtDW . cpiqtbmrnh .
AQHandles , ssGetTaskTime ( rtS , 0 ) ) ; sdiAsyncRepoSetSignalExportSettings
( rtDW . cpiqtbmrnh . AQHandles , 1 , 0 ) ; sdiAsyncRepoSetSignalExportName (
rtDW . cpiqtbmrnh . AQHandles , loggedName , origSigName , propName ) ;
sdiAsyncRepoSetBlockPathDomain ( rtDW . cpiqtbmrnh . AQHandles ) ; }
sdiFreeLabel ( sigName ) ; sdiFreeLabel ( loggedName ) ; sdiFreeLabel (
origSigName ) ; sdiFreeLabel ( propName ) ; sdiFreeLabel ( blockPath ) ;
sdiFreeLabel ( blockSID ) ; sdiFreeLabel ( subPath ) ; } } if ( !
isStreamoutAlreadyRegistered ) { { sdiLabelU varName = sdiGetLabelFromChars (
"KE" ) ; sdiRegisterWksVariable ( rtDW . cpiqtbmrnh . AQHandles , varName ,
"timeseries" ) ; sdiFreeLabel ( varName ) ; } } } } } { { { bool
isStreamoutAlreadyRegistered = false ; { sdiSignalSourceInfoU srcInfo ;
sdiLabelU loggedName = sdiGetLabelFromChars ( "Event Time" ) ; sdiLabelU
origSigName = sdiGetLabelFromChars ( "Event Time" ) ; sdiLabelU propName =
sdiGetLabelFromChars ( "Event Time" ) ; sdiLabelU blockPath =
sdiGetLabelFromChars ( "Copy_of_trial_with_intermediate_tanks/To Workspace3"
) ; sdiLabelU blockSID = sdiGetLabelFromChars ( "" ) ; sdiLabelU subPath =
sdiGetLabelFromChars ( "" ) ; sdiDims sigDims ; sdiLabelU sigName =
sdiGetLabelFromChars ( "Event Time" ) ; sdiAsyncRepoDataTypeHandle hDT =
sdiAsyncRepoGetBuiltInDataTypeHandle ( DATA_TYPE_DOUBLE ) ; { sdiComplexity
sigComplexity = REAL ; sdiSampleTimeContinuity stCont = SAMPLE_TIME_DISCRETE
; int_T sigDimsArray [ 1 ] = { 1 } ; sigDims . nDims = 1 ; sigDims .
dimensions = sigDimsArray ; srcInfo . numBlockPathElems = 1 ; srcInfo .
fullBlockPath = ( sdiFullBlkPathU ) & blockPath ; srcInfo . SID = (
sdiSignalIDU ) & blockSID ; srcInfo . subPath = subPath ; srcInfo . portIndex
= 0 + 1 ; srcInfo . signalName = sigName ; srcInfo . sigSourceUUID = 0 ; rtDW
. bya1bx33so . AQHandles = sdiStartAsyncioQueueCreation ( hDT , & srcInfo ,
rt_dataMapInfo . mmi . InstanceMap . fullPath ,
"27e4af32-3aba-4ddc-bf1f-a1ecad6b0d64" , sigComplexity , & sigDims ,
DIMENSIONS_MODE_FIXED , stCont , "" ) ; sdiCompleteAsyncioQueueCreation (
rtDW . bya1bx33so . AQHandles , hDT , & srcInfo ) ; if ( rtDW . bya1bx33so .
AQHandles ) { sdiSetSignalSampleTimeString ( rtDW . bya1bx33so . AQHandles ,
"Continuous" , 0.0 , ssGetTFinal ( rtS ) ) ; sdiSetSignalRefRate ( rtDW .
bya1bx33so . AQHandles , 0.0 ) ; sdiSetRunStartTime ( rtDW . bya1bx33so .
AQHandles , ssGetTaskTime ( rtS , 1 ) ) ; sdiAsyncRepoSetSignalExportSettings
( rtDW . bya1bx33so . AQHandles , 1 , 0 ) ; sdiAsyncRepoSetSignalExportName (
rtDW . bya1bx33so . AQHandles , loggedName , origSigName , propName ) ;
sdiAsyncRepoSetBlockPathDomain ( rtDW . bya1bx33so . AQHandles ) ; }
sdiFreeLabel ( sigName ) ; sdiFreeLabel ( loggedName ) ; sdiFreeLabel (
origSigName ) ; sdiFreeLabel ( propName ) ; sdiFreeLabel ( blockPath ) ;
sdiFreeLabel ( blockSID ) ; sdiFreeLabel ( subPath ) ; } } if ( !
isStreamoutAlreadyRegistered ) { { sdiLabelU varName = sdiGetLabelFromChars (
"Time_settle" ) ; sdiRegisterWksVariable ( rtDW . bya1bx33so . AQHandles ,
varName , "timeseries" ) ; sdiFreeLabel ( varName ) ; } } } } } { { { bool
isStreamoutAlreadyRegistered = false ; { sdiSignalSourceInfoU srcInfo ;
sdiLabelU loggedName = sdiGetLabelFromChars ( "omega" ) ; sdiLabelU
origSigName = sdiGetLabelFromChars ( "omega" ) ; sdiLabelU propName =
sdiGetLabelFromChars ( "omega" ) ; sdiLabelU blockPath = sdiGetLabelFromChars
( "Copy_of_trial_with_intermediate_tanks/To Workspace4" ) ; sdiLabelU
blockSID = sdiGetLabelFromChars ( "" ) ; sdiLabelU subPath =
sdiGetLabelFromChars ( "" ) ; sdiDims sigDims ; sdiLabelU sigName =
sdiGetLabelFromChars ( "omega" ) ; sdiAsyncRepoDataTypeHandle hDT =
sdiAsyncRepoGetBuiltInDataTypeHandle ( DATA_TYPE_DOUBLE ) ; { sdiComplexity
sigComplexity = REAL ; sdiSampleTimeContinuity stCont =
SAMPLE_TIME_CONTINUOUS ; int_T sigDimsArray [ 1 ] = { 1 } ; sigDims . nDims =
1 ; sigDims . dimensions = sigDimsArray ; srcInfo . numBlockPathElems = 1 ;
srcInfo . fullBlockPath = ( sdiFullBlkPathU ) & blockPath ; srcInfo . SID = (
sdiSignalIDU ) & blockSID ; srcInfo . subPath = subPath ; srcInfo . portIndex
= 0 + 1 ; srcInfo . signalName = sigName ; srcInfo . sigSourceUUID = 0 ; rtDW
. h1qcylur1s . AQHandles = sdiStartAsyncioQueueCreation ( hDT , & srcInfo ,
rt_dataMapInfo . mmi . InstanceMap . fullPath ,
"0fc168d7-a242-43a8-b192-ca6c13491d3d" , sigComplexity , & sigDims ,
DIMENSIONS_MODE_FIXED , stCont , "" ) ; sdiCompleteAsyncioQueueCreation (
rtDW . h1qcylur1s . AQHandles , hDT , & srcInfo ) ; if ( rtDW . h1qcylur1s .
AQHandles ) { sdiSetSignalSampleTimeString ( rtDW . h1qcylur1s . AQHandles ,
"Continuous" , 0.0 , ssGetTFinal ( rtS ) ) ; sdiSetSignalRefRate ( rtDW .
h1qcylur1s . AQHandles , 0.0 ) ; sdiSetRunStartTime ( rtDW . h1qcylur1s .
AQHandles , ssGetTaskTime ( rtS , 0 ) ) ; sdiAsyncRepoSetSignalExportSettings
( rtDW . h1qcylur1s . AQHandles , 1 , 0 ) ; sdiAsyncRepoSetSignalExportName (
rtDW . h1qcylur1s . AQHandles , loggedName , origSigName , propName ) ;
sdiAsyncRepoSetBlockPathDomain ( rtDW . h1qcylur1s . AQHandles ) ; }
sdiFreeLabel ( sigName ) ; sdiFreeLabel ( loggedName ) ; sdiFreeLabel (
origSigName ) ; sdiFreeLabel ( propName ) ; sdiFreeLabel ( blockPath ) ;
sdiFreeLabel ( blockSID ) ; sdiFreeLabel ( subPath ) ; } } if ( !
isStreamoutAlreadyRegistered ) { { sdiLabelU varName = sdiGetLabelFromChars (
"omega" ) ; sdiRegisterWksVariable ( rtDW . h1qcylur1s . AQHandles , varName
, "timeseries" ) ; sdiFreeLabel ( varName ) ; } } } } } { { { bool
isStreamoutAlreadyRegistered = false ; { sdiSignalSourceInfoU srcInfo ;
sdiLabelU loggedName = sdiGetLabelFromChars ( "Gain1" ) ; sdiLabelU
origSigName = sdiGetLabelFromChars ( "" ) ; sdiLabelU propName =
sdiGetLabelFromChars ( "Gain1" ) ; sdiLabelU blockPath = sdiGetLabelFromChars
( "Copy_of_trial_with_intermediate_tanks/To Workspace5" ) ; sdiLabelU
blockSID = sdiGetLabelFromChars ( "" ) ; sdiLabelU subPath =
sdiGetLabelFromChars ( "" ) ; sdiDims sigDims ; sdiLabelU sigName =
sdiGetLabelFromChars ( "Gain1" ) ; sdiAsyncRepoDataTypeHandle hDT =
sdiAsyncRepoGetBuiltInDataTypeHandle ( DATA_TYPE_DOUBLE ) ; { sdiComplexity
sigComplexity = REAL ; sdiSampleTimeContinuity stCont =
SAMPLE_TIME_CONTINUOUS ; int_T sigDimsArray [ 1 ] = { 1 } ; sigDims . nDims =
1 ; sigDims . dimensions = sigDimsArray ; srcInfo . numBlockPathElems = 1 ;
srcInfo . fullBlockPath = ( sdiFullBlkPathU ) & blockPath ; srcInfo . SID = (
sdiSignalIDU ) & blockSID ; srcInfo . subPath = subPath ; srcInfo . portIndex
= 0 + 1 ; srcInfo . signalName = sigName ; srcInfo . sigSourceUUID = 0 ; rtDW
. ohqxmkyead . AQHandles = sdiStartAsyncioQueueCreation ( hDT , & srcInfo ,
rt_dataMapInfo . mmi . InstanceMap . fullPath ,
"386511ed-d1b2-4326-95e7-ebf1fe545d99" , sigComplexity , & sigDims ,
DIMENSIONS_MODE_FIXED , stCont , "" ) ; sdiCompleteAsyncioQueueCreation (
rtDW . ohqxmkyead . AQHandles , hDT , & srcInfo ) ; if ( rtDW . ohqxmkyead .
AQHandles ) { sdiSetSignalSampleTimeString ( rtDW . ohqxmkyead . AQHandles ,
"Continuous" , 0.0 , ssGetTFinal ( rtS ) ) ; sdiSetSignalRefRate ( rtDW .
ohqxmkyead . AQHandles , 0.0 ) ; sdiSetRunStartTime ( rtDW . ohqxmkyead .
AQHandles , ssGetTaskTime ( rtS , 0 ) ) ; sdiAsyncRepoSetSignalExportSettings
( rtDW . ohqxmkyead . AQHandles , 1 , 0 ) ; sdiAsyncRepoSetSignalExportName (
rtDW . ohqxmkyead . AQHandles , loggedName , origSigName , propName ) ;
sdiAsyncRepoSetBlockPathDomain ( rtDW . ohqxmkyead . AQHandles ) ; }
sdiFreeLabel ( sigName ) ; sdiFreeLabel ( loggedName ) ; sdiFreeLabel (
origSigName ) ; sdiFreeLabel ( propName ) ; sdiFreeLabel ( blockPath ) ;
sdiFreeLabel ( blockSID ) ; sdiFreeLabel ( subPath ) ; } } if ( !
isStreamoutAlreadyRegistered ) { { sdiLabelU varName = sdiGetLabelFromChars (
"mech_power" ) ; sdiRegisterWksVariable ( rtDW . ohqxmkyead . AQHandles ,
varName , "timeseries" ) ; sdiFreeLabel ( varName ) ; } } } } } { { { bool
isStreamoutAlreadyRegistered = false ; { sdiSignalSourceInfoU srcInfo ;
sdiLabelU loggedName = sdiGetLabelFromChars ( "Product2" ) ; sdiLabelU
origSigName = sdiGetLabelFromChars ( "" ) ; sdiLabelU propName =
sdiGetLabelFromChars ( "Product2" ) ; sdiLabelU blockPath =
sdiGetLabelFromChars ( "Copy_of_trial_with_intermediate_tanks/To Workspace6"
) ; sdiLabelU blockSID = sdiGetLabelFromChars ( "" ) ; sdiLabelU subPath =
sdiGetLabelFromChars ( "" ) ; sdiDims sigDims ; sdiLabelU sigName =
sdiGetLabelFromChars ( "Product2" ) ; sdiAsyncRepoDataTypeHandle hDT =
sdiAsyncRepoGetBuiltInDataTypeHandle ( DATA_TYPE_DOUBLE ) ; { sdiComplexity
sigComplexity = REAL ; sdiSampleTimeContinuity stCont =
SAMPLE_TIME_CONTINUOUS ; int_T sigDimsArray [ 1 ] = { 1 } ; sigDims . nDims =
1 ; sigDims . dimensions = sigDimsArray ; srcInfo . numBlockPathElems = 1 ;
srcInfo . fullBlockPath = ( sdiFullBlkPathU ) & blockPath ; srcInfo . SID = (
sdiSignalIDU ) & blockSID ; srcInfo . subPath = subPath ; srcInfo . portIndex
= 0 + 1 ; srcInfo . signalName = sigName ; srcInfo . sigSourceUUID = 0 ; rtDW
. jnhubbfo5j . AQHandles = sdiStartAsyncioQueueCreation ( hDT , & srcInfo ,
rt_dataMapInfo . mmi . InstanceMap . fullPath ,
"fa7b845f-9b0c-4b12-b9ad-70fd4cea686f" , sigComplexity , & sigDims ,
DIMENSIONS_MODE_FIXED , stCont , "" ) ; sdiCompleteAsyncioQueueCreation (
rtDW . jnhubbfo5j . AQHandles , hDT , & srcInfo ) ; if ( rtDW . jnhubbfo5j .
AQHandles ) { sdiSetSignalSampleTimeString ( rtDW . jnhubbfo5j . AQHandles ,
"Continuous" , 0.0 , ssGetTFinal ( rtS ) ) ; sdiSetSignalRefRate ( rtDW .
jnhubbfo5j . AQHandles , 0.0 ) ; sdiSetRunStartTime ( rtDW . jnhubbfo5j .
AQHandles , ssGetTaskTime ( rtS , 0 ) ) ; sdiAsyncRepoSetSignalExportSettings
( rtDW . jnhubbfo5j . AQHandles , 1 , 0 ) ; sdiAsyncRepoSetSignalExportName (
rtDW . jnhubbfo5j . AQHandles , loggedName , origSigName , propName ) ;
sdiAsyncRepoSetBlockPathDomain ( rtDW . jnhubbfo5j . AQHandles ) ; }
sdiFreeLabel ( sigName ) ; sdiFreeLabel ( loggedName ) ; sdiFreeLabel (
origSigName ) ; sdiFreeLabel ( propName ) ; sdiFreeLabel ( blockPath ) ;
sdiFreeLabel ( blockSID ) ; sdiFreeLabel ( subPath ) ; } } if ( !
isStreamoutAlreadyRegistered ) { { sdiLabelU varName = sdiGetLabelFromChars (
"Regen_power" ) ; sdiRegisterWksVariable ( rtDW . jnhubbfo5j . AQHandles ,
varName , "timeseries" ) ; sdiFreeLabel ( varName ) ; } } } } } { { { bool
isStreamoutAlreadyRegistered = false ; { sdiSignalSourceInfoU srcInfo ;
sdiLabelU loggedName = sdiGetLabelFromChars ( "Add2" ) ; sdiLabelU
origSigName = sdiGetLabelFromChars ( "" ) ; sdiLabelU propName =
sdiGetLabelFromChars ( "Add2" ) ; sdiLabelU blockPath = sdiGetLabelFromChars
( "Copy_of_trial_with_intermediate_tanks/To Workspace7" ) ; sdiLabelU
blockSID = sdiGetLabelFromChars ( "" ) ; sdiLabelU subPath =
sdiGetLabelFromChars ( "" ) ; sdiDims sigDims ; sdiLabelU sigName =
sdiGetLabelFromChars ( "Add2" ) ; sdiAsyncRepoDataTypeHandle hDT =
sdiAsyncRepoGetBuiltInDataTypeHandle ( DATA_TYPE_DOUBLE ) ; { sdiComplexity
sigComplexity = REAL ; sdiSampleTimeContinuity stCont =
SAMPLE_TIME_CONTINUOUS ; int_T sigDimsArray [ 1 ] = { 1 } ; sigDims . nDims =
1 ; sigDims . dimensions = sigDimsArray ; srcInfo . numBlockPathElems = 1 ;
srcInfo . fullBlockPath = ( sdiFullBlkPathU ) & blockPath ; srcInfo . SID = (
sdiSignalIDU ) & blockSID ; srcInfo . subPath = subPath ; srcInfo . portIndex
= 0 + 1 ; srcInfo . signalName = sigName ; srcInfo . sigSourceUUID = 0 ; rtDW
. fokz15rol1 . AQHandles = sdiStartAsyncioQueueCreation ( hDT , & srcInfo ,
rt_dataMapInfo . mmi . InstanceMap . fullPath ,
"6ea11855-cd28-4f98-a6c6-5d04cd5505de" , sigComplexity , & sigDims ,
DIMENSIONS_MODE_FIXED , stCont , "" ) ; sdiCompleteAsyncioQueueCreation (
rtDW . fokz15rol1 . AQHandles , hDT , & srcInfo ) ; if ( rtDW . fokz15rol1 .
AQHandles ) { sdiSetSignalSampleTimeString ( rtDW . fokz15rol1 . AQHandles ,
"Continuous" , 0.0 , ssGetTFinal ( rtS ) ) ; sdiSetSignalRefRate ( rtDW .
fokz15rol1 . AQHandles , 0.0 ) ; sdiSetRunStartTime ( rtDW . fokz15rol1 .
AQHandles , ssGetTaskTime ( rtS , 0 ) ) ; sdiAsyncRepoSetSignalExportSettings
( rtDW . fokz15rol1 . AQHandles , 1 , 0 ) ; sdiAsyncRepoSetSignalExportName (
rtDW . fokz15rol1 . AQHandles , loggedName , origSigName , propName ) ;
sdiAsyncRepoSetBlockPathDomain ( rtDW . fokz15rol1 . AQHandles ) ; }
sdiFreeLabel ( sigName ) ; sdiFreeLabel ( loggedName ) ; sdiFreeLabel (
origSigName ) ; sdiFreeLabel ( propName ) ; sdiFreeLabel ( blockPath ) ;
sdiFreeLabel ( blockSID ) ; sdiFreeLabel ( subPath ) ; } } if ( !
isStreamoutAlreadyRegistered ) { { sdiLabelU varName = sdiGetLabelFromChars (
"Losses" ) ; sdiRegisterWksVariable ( rtDW . fokz15rol1 . AQHandles , varName
, "timeseries" ) ; sdiFreeLabel ( varName ) ; } } } } } { { { bool
isStreamoutAlreadyRegistered = false ; { sdiSignalSourceInfoU srcInfo ;
sdiLabelU loggedName = sdiGetLabelFromChars ( "Manual Switch" ) ; sdiLabelU
origSigName = sdiGetLabelFromChars ( "" ) ; sdiLabelU propName =
sdiGetLabelFromChars ( "Manual Switch" ) ; sdiLabelU blockPath =
sdiGetLabelFromChars ( "Copy_of_trial_with_intermediate_tanks/To Workspace8"
) ; sdiLabelU blockSID = sdiGetLabelFromChars ( "" ) ; sdiLabelU subPath =
sdiGetLabelFromChars ( "" ) ; sdiDims sigDims ; sdiLabelU sigName =
sdiGetLabelFromChars ( "Manual Switch" ) ; sdiAsyncRepoDataTypeHandle hDT =
sdiAsyncRepoGetBuiltInDataTypeHandle ( DATA_TYPE_DOUBLE ) ; { sdiComplexity
sigComplexity = REAL ; sdiSampleTimeContinuity stCont =
SAMPLE_TIME_CONTINUOUS ; int_T sigDimsArray [ 1 ] = { 1 } ; sigDims . nDims =
1 ; sigDims . dimensions = sigDimsArray ; srcInfo . numBlockPathElems = 1 ;
srcInfo . fullBlockPath = ( sdiFullBlkPathU ) & blockPath ; srcInfo . SID = (
sdiSignalIDU ) & blockSID ; srcInfo . subPath = subPath ; srcInfo . portIndex
= 0 + 1 ; srcInfo . signalName = sigName ; srcInfo . sigSourceUUID = 0 ; rtDW
. lj0dkmhahi . AQHandles = sdiStartAsyncioQueueCreation ( hDT , & srcInfo ,
rt_dataMapInfo . mmi . InstanceMap . fullPath ,
"5223a78e-98a0-4f9d-9969-7734d90db92c" , sigComplexity , & sigDims ,
DIMENSIONS_MODE_FIXED , stCont , "" ) ; sdiCompleteAsyncioQueueCreation (
rtDW . lj0dkmhahi . AQHandles , hDT , & srcInfo ) ; if ( rtDW . lj0dkmhahi .
AQHandles ) { sdiSetSignalSampleTimeString ( rtDW . lj0dkmhahi . AQHandles ,
"Continuous" , 0.0 , ssGetTFinal ( rtS ) ) ; sdiSetSignalRefRate ( rtDW .
lj0dkmhahi . AQHandles , 0.0 ) ; sdiSetRunStartTime ( rtDW . lj0dkmhahi .
AQHandles , ssGetTaskTime ( rtS , 0 ) ) ; sdiAsyncRepoSetSignalExportSettings
( rtDW . lj0dkmhahi . AQHandles , 1 , 0 ) ; sdiAsyncRepoSetSignalExportName (
rtDW . lj0dkmhahi . AQHandles , loggedName , origSigName , propName ) ;
sdiAsyncRepoSetBlockPathDomain ( rtDW . lj0dkmhahi . AQHandles ) ; }
sdiFreeLabel ( sigName ) ; sdiFreeLabel ( loggedName ) ; sdiFreeLabel (
origSigName ) ; sdiFreeLabel ( propName ) ; sdiFreeLabel ( blockPath ) ;
sdiFreeLabel ( blockSID ) ; sdiFreeLabel ( subPath ) ; } } if ( !
isStreamoutAlreadyRegistered ) { { sdiLabelU varName = sdiGetLabelFromChars (
"xdot" ) ; sdiRegisterWksVariable ( rtDW . lj0dkmhahi . AQHandles , varName ,
"timeseries" ) ; sdiFreeLabel ( varName ) ; } } } } } { { { bool
isStreamoutAlreadyRegistered = false ; { sdiSignalSourceInfoU srcInfo ;
sdiLabelU loggedName = sdiGetLabelFromChars ( "Displacement" ) ; sdiLabelU
origSigName = sdiGetLabelFromChars ( "Displacement" ) ; sdiLabelU propName =
sdiGetLabelFromChars ( "Displacement" ) ; sdiLabelU blockPath =
sdiGetLabelFromChars ( "Copy_of_trial_with_intermediate_tanks/To Workspace9"
) ; sdiLabelU blockSID = sdiGetLabelFromChars ( "" ) ; sdiLabelU subPath =
sdiGetLabelFromChars ( "" ) ; sdiDims sigDims ; sdiLabelU sigName =
sdiGetLabelFromChars ( "Displacement" ) ; sdiAsyncRepoDataTypeHandle hDT =
sdiAsyncRepoGetBuiltInDataTypeHandle ( DATA_TYPE_DOUBLE ) ; { sdiComplexity
sigComplexity = REAL ; sdiSampleTimeContinuity stCont =
SAMPLE_TIME_CONTINUOUS ; int_T sigDimsArray [ 1 ] = { 1 } ; sigDims . nDims =
1 ; sigDims . dimensions = sigDimsArray ; srcInfo . numBlockPathElems = 1 ;
srcInfo . fullBlockPath = ( sdiFullBlkPathU ) & blockPath ; srcInfo . SID = (
sdiSignalIDU ) & blockSID ; srcInfo . subPath = subPath ; srcInfo . portIndex
= 0 + 1 ; srcInfo . signalName = sigName ; srcInfo . sigSourceUUID = 0 ; rtDW
. oqdb4h0bik . AQHandles = sdiStartAsyncioQueueCreation ( hDT , & srcInfo ,
rt_dataMapInfo . mmi . InstanceMap . fullPath ,
"50d7c59b-c1a0-4eaa-bed3-6fbbe5fec9f1" , sigComplexity , & sigDims ,
DIMENSIONS_MODE_FIXED , stCont , "" ) ; sdiCompleteAsyncioQueueCreation (
rtDW . oqdb4h0bik . AQHandles , hDT , & srcInfo ) ; if ( rtDW . oqdb4h0bik .
AQHandles ) { sdiSetSignalSampleTimeString ( rtDW . oqdb4h0bik . AQHandles ,
"Continuous" , 0.0 , ssGetTFinal ( rtS ) ) ; sdiSetSignalRefRate ( rtDW .
oqdb4h0bik . AQHandles , 0.0 ) ; sdiSetRunStartTime ( rtDW . oqdb4h0bik .
AQHandles , ssGetTaskTime ( rtS , 0 ) ) ; sdiAsyncRepoSetSignalExportSettings
( rtDW . oqdb4h0bik . AQHandles , 1 , 0 ) ; sdiAsyncRepoSetSignalExportName (
rtDW . oqdb4h0bik . AQHandles , loggedName , origSigName , propName ) ;
sdiAsyncRepoSetBlockPathDomain ( rtDW . oqdb4h0bik . AQHandles ) ; }
sdiFreeLabel ( sigName ) ; sdiFreeLabel ( loggedName ) ; sdiFreeLabel (
origSigName ) ; sdiFreeLabel ( propName ) ; sdiFreeLabel ( blockPath ) ;
sdiFreeLabel ( blockSID ) ; sdiFreeLabel ( subPath ) ; } } if ( !
isStreamoutAlreadyRegistered ) { { sdiLabelU varName = sdiGetLabelFromChars (
"x" ) ; sdiRegisterWksVariable ( rtDW . oqdb4h0bik . AQHandles , varName ,
"timeseries" ) ; sdiFreeLabel ( varName ) ; } } } } } rtB . edzhheuutq = rtP
. InitialSpeed_Value ; rtB . fgm5y4v5wu = rtP . param . P_H ; rtB .
ciwvymzmbd = rtP . param . P_M ; rtB . bb2zjb32zc = rtP . param . P_M ;
MdlInitialize ( ) ; } void MdlOutputs ( int_T tid ) { real_T P0 ; real_T Pg ;
ZCEventType zcEvent ; srClearBC ( rtDW . iln55qrgxk ) ; rtB . mzfn3lhia4 =
rtX . eer5ulgpq1 ; rtB . jswgtkntlr = rtP . Gain_Gain * rtB . mzfn3lhia4 ; if
( ssIsModeUpdateTimeStep ( rtS ) ) { if ( rtDW . aku4xs4fdc != 0 ) { rtX .
gvqywnxquo = rtB . bb2zjb32zc ; } rtB . pyky0iaydb = rtX . gvqywnxquo ; }
else { rtB . pyky0iaydb = rtX . gvqywnxquo ; } rtB . cwl31aims5 = rtB .
pyky0iaydb - rtP . param . P_H ; if ( ssIsModeUpdateTimeStep ( rtS ) ) { rtDW
. awwuosv0rp = ( rtB . cwl31aims5 >= 0.0 ) ; } if ( rtDW . awwuosv0rp > 0 ) {
rtB . dk50fnzsks = rtB . cwl31aims5 ; } else { rtB . dk50fnzsks = - rtB .
cwl31aims5 ; } if ( ssIsModeUpdateTimeStep ( rtS ) ) { rtDW . j5zzepflri = (
( rtB . dk50fnzsks >= rtP . Relay_OnVal ) || ( ( ! ( rtB . dk50fnzsks <= rtP
. Relay_OffVal ) ) && rtDW . j5zzepflri ) ) ; } if ( rtDW . j5zzepflri ) {
rtB . lr3evndp05 = rtP . Relay_YOn ; } else { rtB . lr3evndp05 = rtP .
Relay_YOff ; } rtB . gaymsnblvn = rtX . mui3hk2kdn ; if (
ssIsModeUpdateTimeStep ( rtS ) ) { if ( rtDW . aigzt14viy != 0 ) { rtX .
kllv1nez4n = rtB . fgm5y4v5wu ; } rtB . hfumqnwosi = rtX . kllv1nez4n ; if (
rtDW . jwuiy2atz3 != 0 ) { rtX . orxgnzk15p = rtB . ciwvymzmbd ; } rtB .
jcnaw1aobb = rtX . orxgnzk15p ; } else { rtB . hfumqnwosi = rtX . kllv1nez4n
; rtB . jcnaw1aobb = rtX . orxgnzk15p ; } rtDW . loaiy2cacx = bjcsxykqxt ; Pg
= rtB . hfumqnwosi - rtP . param . P_L ; P0 = rtP . param . P_H - rtP . param
. P_L ; rtB . m0rkotwez0 = ( ( muDoubleScalarExp ( Pg / rtP . param . beta )
- ( Pg / rtP . param . beta + 1.0 ) ) * rtP . param . beta - (
muDoubleScalarExp ( P0 / rtP . param . beta ) - ( P0 / rtP . param . beta +
1.0 ) ) * rtP . param . beta ) * rtP . param . V3_0 ; rtDW . kvbq5jwbvw =
bjcsxykqxt ; Pg = rtB . jcnaw1aobb - rtP . param . P_L ; P0 = rtP . param .
P_M - rtP . param . P_L ; rtB . gnjhdq1qvh = ( ( muDoubleScalarExp ( Pg / rtP
. param . beta ) - ( Pg / rtP . param . beta + 1.0 ) ) * rtP . param . beta -
( muDoubleScalarExp ( P0 / rtP . param . beta ) - ( P0 / rtP . param . beta +
1.0 ) ) * rtP . param . beta ) * rtP . param . V4_0 ; rtB . c5jlsr0y2m = rtX
. apttuwn22n ; rtB . gjajnpeavu = rtX . avebadcay4 ; rtB . onszijeutj = rtX .
o2kfetkxay ; rtB . o5jgaxymjd = rtX . fvfvpkek2n ; rtB . aqvwe35511 = rtX .
aqqzk12fdv ; rtB . ga4spwnzmb = rtX . e30fqjkckb ; rtB . iqcx4ji2gy = rtX .
mizgi5ztop ; Pg = ( ( rtB . o5jgaxymjd + rtB . aqvwe35511 ) + rtB .
ga4spwnzmb ) + rtB . iqcx4ji2gy ; if ( ssIsModeUpdateTimeStep ( rtS ) ) { if
( rtDW . isc2qgvxqr != 0 ) { rtX . iikwtjyqxf = rtB . edzhheuutq ; } rtB .
b5jzfqhswx = rtX . iikwtjyqxf ; } else { rtB . b5jzfqhswx = rtX . iikwtjyqxf
; } rtB . o5psjsgjbl = ( rtP . param . J_elec + rtP . param . J_hyd ) * ( rtB
. b5jzfqhswx * rtB . b5jzfqhswx ) * rtP . Gain3_Gain ; rtB . plcvs3tnfq = ( (
( ( rtB . m0rkotwez0 + rtB . gnjhdq1qvh ) + rtB . gjajnpeavu ) + rtB .
onszijeutj ) + Pg ) + rtB . o5psjsgjbl ; rtB . j0krl2hmjm = rtX . c5bbd2cwaf
; rtB . l0b3kthfn5 = rtB . ci3hl3yrta + rtB . j0krl2hmjm ; rtB . g5clq5nhx1 =
rtX . loqzkyynb5 ; rtB . mfk254enll = rtP . TransferFcn_C [ 0 ] * rtX .
klye04tpmo [ 0 ] ; rtB . mfk254enll += rtP . TransferFcn_C [ 1 ] * rtX .
klye04tpmo [ 1 ] ; if ( ssIsModeUpdateTimeStep ( rtS ) ) { if ( rtB .
mfk254enll >= rtP . param . max_Avt ) { rtDW . g40nougyde = 1 ; } else if (
rtB . mfk254enll > rtP . Saturation_LowerSat ) { rtDW . g40nougyde = 0 ; }
else { rtDW . g40nougyde = - 1 ; } } if ( rtDW . g40nougyde == 1 ) { rtB .
id3dy3zvhz = rtP . param . max_Avt ; } else if ( rtDW . g40nougyde == - 1 ) {
rtB . id3dy3zvhz = rtP . Saturation_LowerSat ; } else { rtB . id3dy3zvhz =
rtB . mfk254enll ; } rtDW . apkag20puf = bjcsxykqxt ; P0 = rtP . param . P_H
- rtB . hfumqnwosi ; rtB . aufytnaglh = rtP . param . Cd * rtB . id3dy3zvhz *
muDoubleScalarSqrt ( muDoubleScalarAbs ( P0 ) ) * muDoubleScalarSign ( P0 ) ;
rtB . exiuxpqhzf = rtP . TransferFcn_C_nnurr522qa [ 0 ] * rtX . ht5behjiqe [
0 ] ; rtB . exiuxpqhzf += rtP . TransferFcn_C_nnurr522qa [ 1 ] * rtX .
ht5behjiqe [ 1 ] ; if ( ssIsModeUpdateTimeStep ( rtS ) ) { if ( rtB .
exiuxpqhzf >= rtP . param . max_Avt ) { rtDW . cangmlit2r = 1 ; } else if (
rtB . exiuxpqhzf > rtP . Saturation_LowerSat_l0ksmwxpmq ) { rtDW . cangmlit2r
= 0 ; } else { rtDW . cangmlit2r = - 1 ; } } if ( rtDW . cangmlit2r == 1 ) {
rtB . jjqotddbo4 = rtP . param . max_Avt ; } else if ( rtDW . cangmlit2r == -
1 ) { rtB . jjqotddbo4 = rtP . Saturation_LowerSat_l0ksmwxpmq ; } else { rtB
. jjqotddbo4 = rtB . exiuxpqhzf ; } rtDW . glwrm0q1go = bjcsxykqxt ; P0 = rtB
. hfumqnwosi - rtB . pyky0iaydb ; rtB . miudztndok = rtP . param . Cd * rtB .
jjqotddbo4 * muDoubleScalarSqrt ( muDoubleScalarAbs ( P0 ) ) *
muDoubleScalarSign ( P0 ) ; rtB . hwx2zv05iw = rtP . TransferFcn_C_j01bq5y01n
[ 0 ] * rtX . aa2dkshjtd [ 0 ] ; rtB . hwx2zv05iw += rtP .
TransferFcn_C_j01bq5y01n [ 1 ] * rtX . aa2dkshjtd [ 1 ] ; if (
ssIsModeUpdateTimeStep ( rtS ) ) { if ( rtB . hwx2zv05iw >= rtP . param .
max_Avt ) { rtDW . od3s0wtt1r = 1 ; } else if ( rtB . hwx2zv05iw > rtP .
Saturation_LowerSat_m0gx2ayxie ) { rtDW . od3s0wtt1r = 0 ; } else { rtDW .
od3s0wtt1r = - 1 ; } } if ( rtDW . od3s0wtt1r == 1 ) { rtB . n3lttcqcta = rtP
. param . max_Avt ; } else if ( rtDW . od3s0wtt1r == - 1 ) { rtB . n3lttcqcta
= rtP . Saturation_LowerSat_m0gx2ayxie ; } else { rtB . n3lttcqcta = rtB .
hwx2zv05iw ; } rtDW . inlqdh4kh1 = bjcsxykqxt ; P0 = rtB . jcnaw1aobb - rtB .
pyky0iaydb ; rtB . eks3p0husl = rtP . param . Cd * rtB . n3lttcqcta *
muDoubleScalarSqrt ( muDoubleScalarAbs ( P0 ) ) * muDoubleScalarSign ( P0 ) ;
rtB . klsduu5tgr = rtP . TransferFcn_C_kavhf31tp1 [ 0 ] * rtX . kco1vq4cxg [
0 ] ; rtB . klsduu5tgr += rtP . TransferFcn_C_kavhf31tp1 [ 1 ] * rtX .
kco1vq4cxg [ 1 ] ; if ( ssIsModeUpdateTimeStep ( rtS ) ) { if ( rtB .
klsduu5tgr >= rtP . param . max_Avt ) { rtDW . ppnlaclhj1 = 1 ; } else if (
rtB . klsduu5tgr > rtP . Saturation_LowerSat_l5sm1vdo10 ) { rtDW . ppnlaclhj1
= 0 ; } else { rtDW . ppnlaclhj1 = - 1 ; } } if ( rtDW . ppnlaclhj1 == 1 ) {
rtB . huebprkurf = rtP . param . max_Avt ; } else if ( rtDW . ppnlaclhj1 == -
1 ) { rtB . huebprkurf = rtP . Saturation_LowerSat_l5sm1vdo10 ; } else { rtB
. huebprkurf = rtB . klsduu5tgr ; } rtDW . atxpvhc5n0 = bjcsxykqxt ; P0 = rtP
. param . P_H - rtB . jcnaw1aobb ; rtB . b0rncimetq = rtP . param . Cd * rtB
. huebprkurf * muDoubleScalarSqrt ( muDoubleScalarAbs ( P0 ) ) *
muDoubleScalarSign ( P0 ) ; rtDW . je5l4cyavv = bjcsxykqxt ; rtB . eqrybinqmz
= rtP . param . D / 6.2831853071795862 * rtB . b5jzfqhswx ; rtB . limtet2woe
= ssGetT ( rtS ) ; if ( ssIsSampleHit ( rtS , 1 , 0 ) &&
ssIsModeUpdateTimeStep ( rtS ) ) { zcEvent = rt_ZCFcn ( RISING_ZERO_CROSSING
, & rtPrevZCX . ah03fydvt5 , ( rtB . lr3evndp05 ) ) ; if ( zcEvent !=
NO_ZCEVENT ) { rtB . oguodtbfqn = rtB . limtet2woe ; rtDW . iln55qrgxk = 4 ;
} } rtB . ap3x1johma = rtX . fklbo31ygf ; { if ( rtDW . homjuevopz .
AQHandles && ssGetLogOutput ( rtS ) ) { sdiWriteSignal ( rtDW . homjuevopz .
AQHandles , ssGetTaskTime ( rtS , 0 ) , ( char * ) & rtB . pyky0iaydb + 0 ) ;
} } { if ( rtDW . h1nrd1cc1r . AQHandles && ssGetLogOutput ( rtS ) ) {
sdiWriteSignal ( rtDW . h1nrd1cc1r . AQHandles , ssGetTaskTime ( rtS , 0 ) ,
( char * ) & rtB . gjajnpeavu + 0 ) ; } } if ( ssIsSampleHit ( rtS , 2 , 0 )
) { rtB . iiwlxh5nit = rtDW . j0cswhuhs0 ; } rtDW . ebwtwofxpd = bjcsxykqxt ;
if ( rtB . lr3evndp05 == 1.0 ) { if ( rtB . b5jzfqhswx < 200.0 ) { rtB .
hx5qshmezi = rtB . iiwlxh5nit ; } else { rtB . hx5qshmezi = 5.0 ; } } else if
( rtB . b5jzfqhswx >= 300.0 ) { rtB . hx5qshmezi = rtB . iiwlxh5nit ; } else
{ rtB . hx5qshmezi = - 5.0 ; } { if ( rtDW . l444m1qlra . AQHandles &&
ssGetLogOutput ( rtS ) ) { sdiWriteSignal ( rtDW . l444m1qlra . AQHandles ,
ssGetTaskTime ( rtS , 0 ) , ( char * ) & rtB . hx5qshmezi + 0 ) ; } } { if (
rtDW . fcjl52ojvi . AQHandles && ssGetLogOutput ( rtS ) ) { sdiWriteSignal (
rtDW . fcjl52ojvi . AQHandles , ssGetTaskTime ( rtS , 0 ) , ( char * ) & rtB
. j0krl2hmjm + 0 ) ; } } { if ( rtDW . cpiqtbmrnh . AQHandles &&
ssGetLogOutput ( rtS ) ) { sdiWriteSignal ( rtDW . cpiqtbmrnh . AQHandles ,
ssGetTaskTime ( rtS , 0 ) , ( char * ) & rtB . o5psjsgjbl + 0 ) ; } } if (
ssIsSampleHit ( rtS , 1 , 0 ) ) { { if ( rtDW . bya1bx33so . AQHandles &&
ssGetLogOutput ( rtS ) ) { sdiWriteSignal ( rtDW . bya1bx33so . AQHandles ,
ssGetTaskTime ( rtS , 1 ) , ( char * ) & rtB . oguodtbfqn + 0 ) ; } } } { if
( rtDW . h1qcylur1s . AQHandles && ssGetLogOutput ( rtS ) ) { sdiWriteSignal
( rtDW . h1qcylur1s . AQHandles , ssGetTaskTime ( rtS , 0 ) , ( char * ) &
rtB . b5jzfqhswx + 0 ) ; } } rtDW . hamy5qs5qe = bjcsxykqxt ; P0 = rtP .
param . D / 6.2831853071795862 * ( rtB . hfumqnwosi - rtB . jcnaw1aobb ) ;
rtB . jkrmpymttx = 1.0 / ( rtP . param . J_elec + rtP . param . J_hyd ) * (
P0 - rtB . hx5qshmezi ) ; rtB . mrxnnzezb3 = P0 ; rtB . eiqmkhbtyt = rtB .
jkrmpymttx * rtB . b5jzfqhswx ; rtB . btfjf0eyh0 = ( rtP . param . J_elec +
rtP . param . J_hyd ) * rtB . eiqmkhbtyt ; { if ( rtDW . ohqxmkyead .
AQHandles && ssGetLogOutput ( rtS ) ) { sdiWriteSignal ( rtDW . ohqxmkyead .
AQHandles , ssGetTaskTime ( rtS , 0 ) , ( char * ) & rtB . btfjf0eyh0 + 0 ) ;
} } rtB . cbizxpnxe4 = rtB . hx5qshmezi * rtB . b5jzfqhswx ; { if ( rtDW .
jnhubbfo5j . AQHandles && ssGetLogOutput ( rtS ) ) { sdiWriteSignal ( rtDW .
jnhubbfo5j . AQHandles , ssGetTaskTime ( rtS , 0 ) , ( char * ) & rtB .
cbizxpnxe4 + 0 ) ; } } { if ( rtDW . fokz15rol1 . AQHandles && ssGetLogOutput
( rtS ) ) { sdiWriteSignal ( rtDW . fokz15rol1 . AQHandles , ssGetTaskTime (
rtS , 0 ) , ( char * ) & Pg + 0 ) ; } } if ( rtP .
ManualSwitch_CurrentSetting == 1 ) { rtB . gzd0bvtkjr = muDoubleScalarSin (
rtP . Velocity_Freq * ssGetTaskTime ( rtS , 0 ) + rtP . Velocity_Phase ) *
rtP . Velocity_Amp + rtP . Velocity_Bias ; } else { rtB . gzd0bvtkjr = rtP .
Constant_Value ; } { if ( rtDW . lj0dkmhahi . AQHandles && ssGetLogOutput (
rtS ) ) { sdiWriteSignal ( rtDW . lj0dkmhahi . AQHandles , ssGetTaskTime (
rtS , 0 ) , ( char * ) & rtB . gzd0bvtkjr + 0 ) ; } } rtB . m3tgvpb0je = rtX
. kwoozlzy2w ; { if ( rtDW . oqdb4h0bik . AQHandles && ssGetLogOutput ( rtS )
) { sdiWriteSignal ( rtDW . oqdb4h0bik . AQHandles , ssGetTaskTime ( rtS , 0
) , ( char * ) & rtB . m3tgvpb0je + 0 ) ; } } rtB . kfmldcnjnq = rtB .
miudztndok + rtB . eks3p0husl ; rtB . hdgz1v45gg = rtB . aufytnaglh + rtB .
b0rncimetq ; rtDW . e0qqpzfruh = bjcsxykqxt ; rtB . fsfzk25tf0 = rtP . param
. beta / ( rtP . param . Acap * rtB . m3tgvpb0je + rtP . param . V1_0 ) * ( (
rtB . miudztndok + rtB . eks3p0husl ) - rtP . param . Acap * rtB . gzd0bvtkjr
) ; if ( ssIsSampleHit ( rtS , 1 , 0 ) ) { rtDW . audz3xgmq3 = bjcsxykqxt ;
if ( rtB . lr3evndp05 == 1.0 ) { rtB . hmcqeqif2o = rtP . param . max_Avt ;
rtB . owtjgmweb5 = rtP . param . max_Avt ; rtB . jpwjwhvvnr = 0.0 ; rtB .
pju2ibjy53 = rtP . param . max_Avt ; } else { rtB . hmcqeqif2o = rtP . param
. max_Avt ; rtB . owtjgmweb5 = 0.0 ; rtB . jpwjwhvvnr = rtP . param . max_Avt
; rtB . pju2ibjy53 = 0.0 ; } } rtDW . epcaebihmu = bjcsxykqxt ; Pg = rtB .
pyky0iaydb - rtP . param . P_L ; rtB . k35lmpwnot = ( ( muDoubleScalarExp (
Pg / rtP . param . beta ) - ( Pg / rtP . param . beta + 1.0 ) ) * rtP . param
. beta + Pg ) * ( rtB . miudztndok + rtB . eks3p0husl ) ; rtDW . d34s0zbduf =
bjcsxykqxt ; rtB . folkmnz0t4 = ( muDoubleScalarExp ( ( rtB . hfumqnwosi -
rtP . param . P_L ) / rtP . param . beta ) - 1.0 ) * rtP . param . beta * ( (
rtB . aufytnaglh - rtB . miudztndok ) - rtB . eqrybinqmz ) ; rtDW .
p1ctkshy11 = bjcsxykqxt ; rtB . n2srwnob2z = ( muDoubleScalarExp ( ( rtB .
jcnaw1aobb - rtP . param . P_L ) / rtP . param . beta ) - 1.0 ) * rtP . param
. beta * ( ( rtB . eqrybinqmz + rtB . b0rncimetq ) - rtB . eks3p0husl ) ; rtB
. lzxak0a220 = rtB . pyky0iaydb - rtP . param . P_L ; rtB . faofsw34yz = rtB
. kfmldcnjnq * rtB . lzxak0a220 ; rtB . oz3rry1khh = rtB . hdgz1v45gg * rtP .
param . P_H ; rtB . jhxf0vx510 = rtB . hfumqnwosi - rtB . pyky0iaydb ; rtB .
ahjuj0uqhz = rtB . jhxf0vx510 * rtB . miudztndok ; rtB . eidzbo2axi = rtB .
jcnaw1aobb - rtB . pyky0iaydb ; rtB . dzr21vacs5 = rtB . eidzbo2axi * rtB .
eks3p0husl ; rtB . g5dn334cgf = rtP . param . P_H - rtB . hfumqnwosi ; rtB .
hllpacqnya = rtB . g5dn334cgf * rtB . aufytnaglh ; rtB . fad5yhgoij = rtP .
param . P_H - rtB . jcnaw1aobb ; rtB . nzchyezcdl = rtB . fad5yhgoij * rtB .
b0rncimetq ; rtDW . a4pqll03x1 = bjcsxykqxt ; rtB . afrqwefb1y = ( ( rtB .
aufytnaglh - rtB . miudztndok ) - rtB . eqrybinqmz ) * ( rtP . param . beta /
rtP . param . V3_0 ) ; rtDW . mpaqmqlwie = bjcsxykqxt ; rtB . djfhnx1214 = (
( rtB . b0rncimetq + rtB . eqrybinqmz ) - rtB . eks3p0husl ) * ( rtP . param
. beta / rtP . param . V4_0 ) ; rtB . ouvddzwv14 = rtP .
TransferFcn_C_lu1wzht1on [ 0 ] * rtX . mvzjw32hje [ 0 ] ; rtB . ouvddzwv14 +=
rtP . TransferFcn_C_lu1wzht1on [ 1 ] * rtX . mvzjw32hje [ 1 ] ; if (
ssIsModeUpdateTimeStep ( rtS ) ) { if ( rtB . ouvddzwv14 >= rtP . param .
max_Avt ) { rtDW . exujqdty5l = 1 ; } else if ( rtB . ouvddzwv14 > rtP .
Saturation_LowerSat_mclkng5niv ) { rtDW . exujqdty5l = 0 ; } else { rtDW .
exujqdty5l = - 1 ; } } if ( rtDW . exujqdty5l == 1 ) { rtB . ixdckelvxz = rtP
. param . max_Avt ; } else if ( rtDW . exujqdty5l == - 1 ) { rtB . ixdckelvxz
= rtP . Saturation_LowerSat_mclkng5niv ; } else { rtB . ixdckelvxz = rtB .
ouvddzwv14 ; } UNUSED_PARAMETER ( tid ) ; } void MdlOutputsTID3 ( int_T tid )
{ rtB . edzhheuutq = rtP . InitialSpeed_Value ; rtB . ci3hl3yrta = ( rtP .
param . J_elec + rtP . param . J_hyd ) * ( rtB . edzhheuutq * rtB .
edzhheuutq ) * rtP . Gain5_Gain ; rtB . fgm5y4v5wu = rtP . param . P_H ; rtB
. ciwvymzmbd = rtP . param . P_M ; rtB . bb2zjb32zc = rtP . param . P_M ;
UNUSED_PARAMETER ( tid ) ; } void MdlUpdate ( int_T tid ) { rtDW . aku4xs4fdc
= 0 ; rtDW . aigzt14viy = 0 ; rtDW . jwuiy2atz3 = 0 ; rtDW . isc2qgvxqr = 0 ;
if ( ssIsSampleHit ( rtS , 2 , 0 ) ) { rtDW . j0cswhuhs0 = rtB . mrxnnzezb3 ;
} UNUSED_PARAMETER ( tid ) ; } void MdlUpdateTID3 ( int_T tid ) {
UNUSED_PARAMETER ( tid ) ; } void MdlDerivatives ( void ) { XDot * _rtXdot ;
_rtXdot = ( ( XDot * ) ssGetdX ( rtS ) ) ; _rtXdot -> eer5ulgpq1 = rtB .
b5jzfqhswx ; _rtXdot -> gvqywnxquo = rtB . fsfzk25tf0 ; _rtXdot -> mui3hk2kdn
= rtB . n2srwnob2z ; _rtXdot -> kllv1nez4n = rtB . afrqwefb1y ; _rtXdot ->
orxgnzk15p = rtB . djfhnx1214 ; _rtXdot -> apttuwn22n = rtB . folkmnz0t4 ;
_rtXdot -> avebadcay4 = rtB . cbizxpnxe4 ; _rtXdot -> o2kfetkxay = rtB .
k35lmpwnot ; _rtXdot -> fvfvpkek2n = rtB . hllpacqnya ; _rtXdot -> aqqzk12fdv
= rtB . ahjuj0uqhz ; _rtXdot -> e30fqjkckb = rtB . dzr21vacs5 ; _rtXdot ->
mizgi5ztop = rtB . nzchyezcdl ; _rtXdot -> iikwtjyqxf = rtB . jkrmpymttx ;
_rtXdot -> c5bbd2cwaf = rtB . oz3rry1khh ; _rtXdot -> loqzkyynb5 = rtB .
btfjf0eyh0 ; _rtXdot -> klye04tpmo [ 0 ] = rtP . TransferFcn_A [ 0 ] * rtX .
klye04tpmo [ 0 ] ; _rtXdot -> klye04tpmo [ 0 ] += rtP . TransferFcn_A [ 1 ] *
rtX . klye04tpmo [ 1 ] ; _rtXdot -> klye04tpmo [ 1 ] = rtX . klye04tpmo [ 0 ]
; _rtXdot -> klye04tpmo [ 0 ] += rtB . hmcqeqif2o ; _rtXdot -> ht5behjiqe [ 0
] = rtP . TransferFcn_A_msxihde20w [ 0 ] * rtX . ht5behjiqe [ 0 ] ; _rtXdot
-> ht5behjiqe [ 0 ] += rtP . TransferFcn_A_msxihde20w [ 1 ] * rtX .
ht5behjiqe [ 1 ] ; _rtXdot -> ht5behjiqe [ 1 ] = rtX . ht5behjiqe [ 0 ] ;
_rtXdot -> ht5behjiqe [ 0 ] += rtB . owtjgmweb5 ; _rtXdot -> aa2dkshjtd [ 0 ]
= rtP . TransferFcn_A_bwuzoqq3iu [ 0 ] * rtX . aa2dkshjtd [ 0 ] ; _rtXdot ->
aa2dkshjtd [ 0 ] += rtP . TransferFcn_A_bwuzoqq3iu [ 1 ] * rtX . aa2dkshjtd [
1 ] ; _rtXdot -> aa2dkshjtd [ 1 ] = rtX . aa2dkshjtd [ 0 ] ; _rtXdot ->
aa2dkshjtd [ 0 ] += rtB . jpwjwhvvnr ; _rtXdot -> kco1vq4cxg [ 0 ] = rtP .
TransferFcn_A_cqt2tjfue3 [ 0 ] * rtX . kco1vq4cxg [ 0 ] ; _rtXdot ->
kco1vq4cxg [ 0 ] += rtP . TransferFcn_A_cqt2tjfue3 [ 1 ] * rtX . kco1vq4cxg [
1 ] ; _rtXdot -> kco1vq4cxg [ 1 ] = rtX . kco1vq4cxg [ 0 ] ; _rtXdot ->
kco1vq4cxg [ 0 ] += rtB . pju2ibjy53 ; _rtXdot -> fklbo31ygf = rtB .
faofsw34yz ; _rtXdot -> kwoozlzy2w = rtB . gzd0bvtkjr ; _rtXdot -> mvzjw32hje
[ 0 ] = rtP . TransferFcn_A_m2qs34lzon [ 0 ] * rtX . mvzjw32hje [ 0 ] ;
_rtXdot -> mvzjw32hje [ 0 ] += rtP . TransferFcn_A_m2qs34lzon [ 1 ] * rtX .
mvzjw32hje [ 1 ] ; _rtXdot -> mvzjw32hje [ 1 ] = rtX . mvzjw32hje [ 0 ] ;
_rtXdot -> mvzjw32hje [ 0 ] += rtP . param . max_Avt ; } void MdlProjection (
void ) { } void MdlZeroCrossings ( void ) { ZCV * _rtZCSV ; _rtZCSV = ( ( ZCV
* ) ssGetSolverZcSignalVector ( rtS ) ) ; _rtZCSV -> krptktw1n2 = rtB .
cwl31aims5 ; if ( rtDW . j5zzepflri ) { _rtZCSV -> py1rahqile = rtB .
dk50fnzsks - rtP . Relay_OffVal ; } else { _rtZCSV -> py1rahqile = rtB .
dk50fnzsks - rtP . Relay_OnVal ; } _rtZCSV -> a0e1ektqfm = rtB . mfk254enll -
rtP . param . max_Avt ; _rtZCSV -> lonlr0pekt = rtB . mfk254enll - rtP .
Saturation_LowerSat ; _rtZCSV -> ghimdihoai = rtB . exiuxpqhzf - rtP . param
. max_Avt ; _rtZCSV -> hen3ks1v4p = rtB . exiuxpqhzf - rtP .
Saturation_LowerSat_l0ksmwxpmq ; _rtZCSV -> kij5fz4qdb = rtB . hwx2zv05iw -
rtP . param . max_Avt ; _rtZCSV -> l4cui55vk5 = rtB . hwx2zv05iw - rtP .
Saturation_LowerSat_m0gx2ayxie ; _rtZCSV -> h2qdorbmuq = rtB . klsduu5tgr -
rtP . param . max_Avt ; _rtZCSV -> pseezcu2gs = rtB . klsduu5tgr - rtP .
Saturation_LowerSat_l5sm1vdo10 ; _rtZCSV -> i3mlogp5q3 = rtB . ouvddzwv14 -
rtP . param . max_Avt ; _rtZCSV -> azuoncxh5k = rtB . ouvddzwv14 - rtP .
Saturation_LowerSat_mclkng5niv ; } void MdlTerminate ( void ) { { if ( rtDW .
homjuevopz . AQHandles ) { sdiTerminateStreaming ( & rtDW . homjuevopz .
AQHandles ) ; } } { if ( rtDW . h1nrd1cc1r . AQHandles ) {
sdiTerminateStreaming ( & rtDW . h1nrd1cc1r . AQHandles ) ; } } { if ( rtDW .
l444m1qlra . AQHandles ) { sdiTerminateStreaming ( & rtDW . l444m1qlra .
AQHandles ) ; } } { if ( rtDW . fcjl52ojvi . AQHandles ) {
sdiTerminateStreaming ( & rtDW . fcjl52ojvi . AQHandles ) ; } } { if ( rtDW .
cpiqtbmrnh . AQHandles ) { sdiTerminateStreaming ( & rtDW . cpiqtbmrnh .
AQHandles ) ; } } { if ( rtDW . bya1bx33so . AQHandles ) {
sdiTerminateStreaming ( & rtDW . bya1bx33so . AQHandles ) ; } } { if ( rtDW .
h1qcylur1s . AQHandles ) { sdiTerminateStreaming ( & rtDW . h1qcylur1s .
AQHandles ) ; } } { if ( rtDW . ohqxmkyead . AQHandles ) {
sdiTerminateStreaming ( & rtDW . ohqxmkyead . AQHandles ) ; } } { if ( rtDW .
jnhubbfo5j . AQHandles ) { sdiTerminateStreaming ( & rtDW . jnhubbfo5j .
AQHandles ) ; } } { if ( rtDW . fokz15rol1 . AQHandles ) {
sdiTerminateStreaming ( & rtDW . fokz15rol1 . AQHandles ) ; } } { if ( rtDW .
lj0dkmhahi . AQHandles ) { sdiTerminateStreaming ( & rtDW . lj0dkmhahi .
AQHandles ) ; } } { if ( rtDW . oqdb4h0bik . AQHandles ) {
sdiTerminateStreaming ( & rtDW . oqdb4h0bik . AQHandles ) ; } } } static void
mr_Copy_of_trial_with_intermediate_tanks_cacheDataAsMxArray ( mxArray *
destArray , mwIndex i , int j , const void * srcData , size_t numBytes ) ;
static void mr_Copy_of_trial_with_intermediate_tanks_cacheDataAsMxArray (
mxArray * destArray , mwIndex i , int j , const void * srcData , size_t
numBytes ) { mxArray * newArray = mxCreateUninitNumericMatrix ( ( size_t ) 1
, numBytes , mxUINT8_CLASS , mxREAL ) ; memcpy ( ( uint8_T * ) mxGetData (
newArray ) , ( const uint8_T * ) srcData , numBytes ) ; mxSetFieldByNumber (
destArray , i , j , newArray ) ; } static void
mr_Copy_of_trial_with_intermediate_tanks_restoreDataFromMxArray ( void *
destData , const mxArray * srcArray , mwIndex i , int j , size_t numBytes ) ;
static void mr_Copy_of_trial_with_intermediate_tanks_restoreDataFromMxArray (
void * destData , const mxArray * srcArray , mwIndex i , int j , size_t
numBytes ) { memcpy ( ( uint8_T * ) destData , ( const uint8_T * ) mxGetData
( mxGetFieldByNumber ( srcArray , i , j ) ) , numBytes ) ; } static void
mr_Copy_of_trial_with_intermediate_tanks_cacheBitFieldToMxArray ( mxArray *
destArray , mwIndex i , int j , uint_T bitVal ) ; static void
mr_Copy_of_trial_with_intermediate_tanks_cacheBitFieldToMxArray ( mxArray *
destArray , mwIndex i , int j , uint_T bitVal ) { mxSetFieldByNumber (
destArray , i , j , mxCreateDoubleScalar ( ( real_T ) bitVal ) ) ; } static
uint_T mr_Copy_of_trial_with_intermediate_tanks_extractBitFieldFromMxArray (
const mxArray * srcArray , mwIndex i , int j , uint_T numBits ) ; static
uint_T mr_Copy_of_trial_with_intermediate_tanks_extractBitFieldFromMxArray (
const mxArray * srcArray , mwIndex i , int j , uint_T numBits ) { const
uint_T varVal = ( uint_T ) mxGetScalar ( mxGetFieldByNumber ( srcArray , i ,
j ) ) ; return varVal & ( ( 1u << numBits ) - 1u ) ; } static void
mr_Copy_of_trial_with_intermediate_tanks_cacheDataToMxArrayWithOffset (
mxArray * destArray , mwIndex i , int j , mwIndex offset , const void *
srcData , size_t numBytes ) ; static void
mr_Copy_of_trial_with_intermediate_tanks_cacheDataToMxArrayWithOffset (
mxArray * destArray , mwIndex i , int j , mwIndex offset , const void *
srcData , size_t numBytes ) { uint8_T * varData = ( uint8_T * ) mxGetData (
mxGetFieldByNumber ( destArray , i , j ) ) ; memcpy ( ( uint8_T * ) & varData
[ offset * numBytes ] , ( const uint8_T * ) srcData , numBytes ) ; } static
void
mr_Copy_of_trial_with_intermediate_tanks_restoreDataFromMxArrayWithOffset (
void * destData , const mxArray * srcArray , mwIndex i , int j , mwIndex
offset , size_t numBytes ) ; static void
mr_Copy_of_trial_with_intermediate_tanks_restoreDataFromMxArrayWithOffset (
void * destData , const mxArray * srcArray , mwIndex i , int j , mwIndex
offset , size_t numBytes ) { const uint8_T * varData = ( const uint8_T * )
mxGetData ( mxGetFieldByNumber ( srcArray , i , j ) ) ; memcpy ( ( uint8_T *
) destData , ( const uint8_T * ) & varData [ offset * numBytes ] , numBytes )
; } static void
mr_Copy_of_trial_with_intermediate_tanks_cacheBitFieldToCellArrayWithOffset (
mxArray * destArray , mwIndex i , int j , mwIndex offset , uint_T fieldVal )
; static void
mr_Copy_of_trial_with_intermediate_tanks_cacheBitFieldToCellArrayWithOffset (
mxArray * destArray , mwIndex i , int j , mwIndex offset , uint_T fieldVal )
{ mxSetCell ( mxGetFieldByNumber ( destArray , i , j ) , offset ,
mxCreateDoubleScalar ( ( real_T ) fieldVal ) ) ; } static uint_T
mr_Copy_of_trial_with_intermediate_tanks_extractBitFieldFromCellArrayWithOffset
( const mxArray * srcArray , mwIndex i , int j , mwIndex offset , uint_T
numBits ) ; static uint_T
mr_Copy_of_trial_with_intermediate_tanks_extractBitFieldFromCellArrayWithOffset
( const mxArray * srcArray , mwIndex i , int j , mwIndex offset , uint_T
numBits ) { const uint_T fieldVal = ( uint_T ) mxGetScalar ( mxGetCell (
mxGetFieldByNumber ( srcArray , i , j ) , offset ) ) ; return fieldVal & ( (
1u << numBits ) - 1u ) ; } mxArray *
mr_Copy_of_trial_with_intermediate_tanks_GetDWork ( ) { static const char_T *
ssDWFieldNames [ 3 ] = { "rtB" , "rtDW" , "rtPrevZCX" , } ; mxArray * ssDW =
mxCreateStructMatrix ( 1 , 1 , 3 , ssDWFieldNames ) ;
mr_Copy_of_trial_with_intermediate_tanks_cacheDataAsMxArray ( ssDW , 0 , 0 ,
( const void * ) & ( rtB ) , sizeof ( rtB ) ) ; { static const char_T *
rtdwDataFieldNames [ 61 ] = { "rtDW.j0cswhuhs0" , "rtDW.atxpvhc5n0" ,
"rtDW.inlqdh4kh1" , "rtDW.glwrm0q1go" , "rtDW.apkag20puf" , "rtDW.mpaqmqlwie"
, "rtDW.a4pqll03x1" , "rtDW.p1ctkshy11" , "rtDW.d34s0zbduf" ,
"rtDW.epcaebihmu" , "rtDW.kvbq5jwbvw" , "rtDW.audz3xgmq3" , "rtDW.loaiy2cacx"
, "rtDW.e0qqpzfruh" , "rtDW.hamy5qs5qe" , "rtDW.je5l4cyavv" ,
"rtDW.ebwtwofxpd" , "rtDW.aku4xs4fdc" , "rtDW.aigzt14viy" , "rtDW.jwuiy2atz3"
, "rtDW.isc2qgvxqr" , "rtDW.awwuosv0rp" , "rtDW.g40nougyde" ,
"rtDW.cangmlit2r" , "rtDW.od3s0wtt1r" , "rtDW.ppnlaclhj1" , "rtDW.exujqdty5l"
, "rtDW.iln55qrgxk" , "rtDW.i213auut4h" , "rtDW.e13oplclke" ,
"rtDW.laquyfl1kp" , "rtDW.f45mx2s03m" , "rtDW.gstdmwepbc" , "rtDW.ikvzyivi0j"
, "rtDW.libi2dsvky" , "rtDW.ga2xty0o4o" , "rtDW.e1b1tx1d4l" ,
"rtDW.p0eq4zkt3w" , "rtDW.oahmxad2ft" , "rtDW.ksvy5jh2vc" , "rtDW.fyaulsukmr"
, "rtDW.hdus1cnzep" , "rtDW.amiin4kdmj" , "rtDW.nwxy1ani3r" ,
"rtDW.j5zzepflri" , "rtDW.hk25kkxqwm" , "rtDW.ffd0olavlx" , "rtDW.eiohhcwuqa"
, "rtDW.eyyp3zs3le" , "rtDW.nakeaw2rna" , "rtDW.n4ne2vzdw3" ,
"rtDW.ljfwvqlys1" , "rtDW.ct4cduudvv" , "rtDW.jkvzly4a1o" , "rtDW.fxawe3htlm"
, "rtDW.irq32iym0t" , "rtDW.jz3v2t25p2" , "rtDW.m0i52xetkk" ,
"rtDW.dfvcjzwy4h" , "rtDW.nas3vtbw5t" , "rtDW.awqg33it32" , } ; mxArray *
rtdwData = mxCreateStructMatrix ( 1 , 1 , 61 , rtdwDataFieldNames ) ;
mr_Copy_of_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 ,
0 , ( const void * ) & ( rtDW . j0cswhuhs0 ) , sizeof ( rtDW . j0cswhuhs0 ) )
; mr_Copy_of_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0
, 1 , ( const void * ) & ( rtDW . atxpvhc5n0 ) , sizeof ( rtDW . atxpvhc5n0 )
) ; mr_Copy_of_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData ,
0 , 2 , ( const void * ) & ( rtDW . inlqdh4kh1 ) , sizeof ( rtDW . inlqdh4kh1
) ) ; mr_Copy_of_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData
, 0 , 3 , ( const void * ) & ( rtDW . glwrm0q1go ) , sizeof ( rtDW .
glwrm0q1go ) ) ; mr_Copy_of_trial_with_intermediate_tanks_cacheDataAsMxArray
( rtdwData , 0 , 4 , ( const void * ) & ( rtDW . apkag20puf ) , sizeof ( rtDW
. apkag20puf ) ) ;
mr_Copy_of_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 ,
5 , ( const void * ) & ( rtDW . mpaqmqlwie ) , sizeof ( rtDW . mpaqmqlwie ) )
; mr_Copy_of_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0
, 6 , ( const void * ) & ( rtDW . a4pqll03x1 ) , sizeof ( rtDW . a4pqll03x1 )
) ; mr_Copy_of_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData ,
0 , 7 , ( const void * ) & ( rtDW . p1ctkshy11 ) , sizeof ( rtDW . p1ctkshy11
) ) ; mr_Copy_of_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData
, 0 , 8 , ( const void * ) & ( rtDW . d34s0zbduf ) , sizeof ( rtDW .
d34s0zbduf ) ) ; mr_Copy_of_trial_with_intermediate_tanks_cacheDataAsMxArray
( rtdwData , 0 , 9 , ( const void * ) & ( rtDW . epcaebihmu ) , sizeof ( rtDW
. epcaebihmu ) ) ;
mr_Copy_of_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 ,
10 , ( const void * ) & ( rtDW . kvbq5jwbvw ) , sizeof ( rtDW . kvbq5jwbvw )
) ; mr_Copy_of_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData ,
0 , 11 , ( const void * ) & ( rtDW . audz3xgmq3 ) , sizeof ( rtDW .
audz3xgmq3 ) ) ; mr_Copy_of_trial_with_intermediate_tanks_cacheDataAsMxArray
( rtdwData , 0 , 12 , ( const void * ) & ( rtDW . loaiy2cacx ) , sizeof (
rtDW . loaiy2cacx ) ) ;
mr_Copy_of_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 ,
13 , ( const void * ) & ( rtDW . e0qqpzfruh ) , sizeof ( rtDW . e0qqpzfruh )
) ; mr_Copy_of_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData ,
0 , 14 , ( const void * ) & ( rtDW . hamy5qs5qe ) , sizeof ( rtDW .
hamy5qs5qe ) ) ; mr_Copy_of_trial_with_intermediate_tanks_cacheDataAsMxArray
( rtdwData , 0 , 15 , ( const void * ) & ( rtDW . je5l4cyavv ) , sizeof (
rtDW . je5l4cyavv ) ) ;
mr_Copy_of_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 ,
16 , ( const void * ) & ( rtDW . ebwtwofxpd ) , sizeof ( rtDW . ebwtwofxpd )
) ; mr_Copy_of_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData ,
0 , 17 , ( const void * ) & ( rtDW . aku4xs4fdc ) , sizeof ( rtDW .
aku4xs4fdc ) ) ; mr_Copy_of_trial_with_intermediate_tanks_cacheDataAsMxArray
( rtdwData , 0 , 18 , ( const void * ) & ( rtDW . aigzt14viy ) , sizeof (
rtDW . aigzt14viy ) ) ;
mr_Copy_of_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 ,
19 , ( const void * ) & ( rtDW . jwuiy2atz3 ) , sizeof ( rtDW . jwuiy2atz3 )
) ; mr_Copy_of_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData ,
0 , 20 , ( const void * ) & ( rtDW . isc2qgvxqr ) , sizeof ( rtDW .
isc2qgvxqr ) ) ; mr_Copy_of_trial_with_intermediate_tanks_cacheDataAsMxArray
( rtdwData , 0 , 21 , ( const void * ) & ( rtDW . awwuosv0rp ) , sizeof (
rtDW . awwuosv0rp ) ) ;
mr_Copy_of_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 ,
22 , ( const void * ) & ( rtDW . g40nougyde ) , sizeof ( rtDW . g40nougyde )
) ; mr_Copy_of_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData ,
0 , 23 , ( const void * ) & ( rtDW . cangmlit2r ) , sizeof ( rtDW .
cangmlit2r ) ) ; mr_Copy_of_trial_with_intermediate_tanks_cacheDataAsMxArray
( rtdwData , 0 , 24 , ( const void * ) & ( rtDW . od3s0wtt1r ) , sizeof (
rtDW . od3s0wtt1r ) ) ;
mr_Copy_of_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 ,
25 , ( const void * ) & ( rtDW . ppnlaclhj1 ) , sizeof ( rtDW . ppnlaclhj1 )
) ; mr_Copy_of_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData ,
0 , 26 , ( const void * ) & ( rtDW . exujqdty5l ) , sizeof ( rtDW .
exujqdty5l ) ) ; mr_Copy_of_trial_with_intermediate_tanks_cacheDataAsMxArray
( rtdwData , 0 , 27 , ( const void * ) & ( rtDW . iln55qrgxk ) , sizeof (
rtDW . iln55qrgxk ) ) ;
mr_Copy_of_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 ,
28 , ( const void * ) & ( rtDW . i213auut4h ) , sizeof ( rtDW . i213auut4h )
) ; mr_Copy_of_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData ,
0 , 29 , ( const void * ) & ( rtDW . e13oplclke ) , sizeof ( rtDW .
e13oplclke ) ) ; mr_Copy_of_trial_with_intermediate_tanks_cacheDataAsMxArray
( rtdwData , 0 , 30 , ( const void * ) & ( rtDW . laquyfl1kp ) , sizeof (
rtDW . laquyfl1kp ) ) ;
mr_Copy_of_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 ,
31 , ( const void * ) & ( rtDW . f45mx2s03m ) , sizeof ( rtDW . f45mx2s03m )
) ; mr_Copy_of_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData ,
0 , 32 , ( const void * ) & ( rtDW . gstdmwepbc ) , sizeof ( rtDW .
gstdmwepbc ) ) ; mr_Copy_of_trial_with_intermediate_tanks_cacheDataAsMxArray
( rtdwData , 0 , 33 , ( const void * ) & ( rtDW . ikvzyivi0j ) , sizeof (
rtDW . ikvzyivi0j ) ) ;
mr_Copy_of_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 ,
34 , ( const void * ) & ( rtDW . libi2dsvky ) , sizeof ( rtDW . libi2dsvky )
) ; mr_Copy_of_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData ,
0 , 35 , ( const void * ) & ( rtDW . ga2xty0o4o ) , sizeof ( rtDW .
ga2xty0o4o ) ) ; mr_Copy_of_trial_with_intermediate_tanks_cacheDataAsMxArray
( rtdwData , 0 , 36 , ( const void * ) & ( rtDW . e1b1tx1d4l ) , sizeof (
rtDW . e1b1tx1d4l ) ) ;
mr_Copy_of_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 ,
37 , ( const void * ) & ( rtDW . p0eq4zkt3w ) , sizeof ( rtDW . p0eq4zkt3w )
) ; mr_Copy_of_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData ,
0 , 38 , ( const void * ) & ( rtDW . oahmxad2ft ) , sizeof ( rtDW .
oahmxad2ft ) ) ; mr_Copy_of_trial_with_intermediate_tanks_cacheDataAsMxArray
( rtdwData , 0 , 39 , ( const void * ) & ( rtDW . ksvy5jh2vc ) , sizeof (
rtDW . ksvy5jh2vc ) ) ;
mr_Copy_of_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 ,
40 , ( const void * ) & ( rtDW . fyaulsukmr ) , sizeof ( rtDW . fyaulsukmr )
) ; mr_Copy_of_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData ,
0 , 41 , ( const void * ) & ( rtDW . hdus1cnzep ) , sizeof ( rtDW .
hdus1cnzep ) ) ; mr_Copy_of_trial_with_intermediate_tanks_cacheDataAsMxArray
( rtdwData , 0 , 42 , ( const void * ) & ( rtDW . amiin4kdmj ) , sizeof (
rtDW . amiin4kdmj ) ) ;
mr_Copy_of_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 ,
43 , ( const void * ) & ( rtDW . nwxy1ani3r ) , sizeof ( rtDW . nwxy1ani3r )
) ; mr_Copy_of_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData ,
0 , 44 , ( const void * ) & ( rtDW . j5zzepflri ) , sizeof ( rtDW .
j5zzepflri ) ) ; mr_Copy_of_trial_with_intermediate_tanks_cacheDataAsMxArray
( rtdwData , 0 , 45 , ( const void * ) & ( rtDW . hk25kkxqwm ) , sizeof (
rtDW . hk25kkxqwm ) ) ;
mr_Copy_of_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 ,
46 , ( const void * ) & ( rtDW . ffd0olavlx ) , sizeof ( rtDW . ffd0olavlx )
) ; mr_Copy_of_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData ,
0 , 47 , ( const void * ) & ( rtDW . eiohhcwuqa ) , sizeof ( rtDW .
eiohhcwuqa ) ) ; mr_Copy_of_trial_with_intermediate_tanks_cacheDataAsMxArray
( rtdwData , 0 , 48 , ( const void * ) & ( rtDW . eyyp3zs3le ) , sizeof (
rtDW . eyyp3zs3le ) ) ;
mr_Copy_of_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 ,
49 , ( const void * ) & ( rtDW . nakeaw2rna ) , sizeof ( rtDW . nakeaw2rna )
) ; mr_Copy_of_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData ,
0 , 50 , ( const void * ) & ( rtDW . n4ne2vzdw3 ) , sizeof ( rtDW .
n4ne2vzdw3 ) ) ; mr_Copy_of_trial_with_intermediate_tanks_cacheDataAsMxArray
( rtdwData , 0 , 51 , ( const void * ) & ( rtDW . ljfwvqlys1 ) , sizeof (
rtDW . ljfwvqlys1 ) ) ;
mr_Copy_of_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 ,
52 , ( const void * ) & ( rtDW . ct4cduudvv ) , sizeof ( rtDW . ct4cduudvv )
) ; mr_Copy_of_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData ,
0 , 53 , ( const void * ) & ( rtDW . jkvzly4a1o ) , sizeof ( rtDW .
jkvzly4a1o ) ) ; mr_Copy_of_trial_with_intermediate_tanks_cacheDataAsMxArray
( rtdwData , 0 , 54 , ( const void * ) & ( rtDW . fxawe3htlm ) , sizeof (
rtDW . fxawe3htlm ) ) ;
mr_Copy_of_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 ,
55 , ( const void * ) & ( rtDW . irq32iym0t ) , sizeof ( rtDW . irq32iym0t )
) ; mr_Copy_of_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData ,
0 , 56 , ( const void * ) & ( rtDW . jz3v2t25p2 ) , sizeof ( rtDW .
jz3v2t25p2 ) ) ; mr_Copy_of_trial_with_intermediate_tanks_cacheDataAsMxArray
( rtdwData , 0 , 57 , ( const void * ) & ( rtDW . m0i52xetkk ) , sizeof (
rtDW . m0i52xetkk ) ) ;
mr_Copy_of_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 ,
58 , ( const void * ) & ( rtDW . dfvcjzwy4h ) , sizeof ( rtDW . dfvcjzwy4h )
) ; mr_Copy_of_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData ,
0 , 59 , ( const void * ) & ( rtDW . nas3vtbw5t ) , sizeof ( rtDW .
nas3vtbw5t ) ) ; mr_Copy_of_trial_with_intermediate_tanks_cacheDataAsMxArray
( rtdwData , 0 , 60 , ( const void * ) & ( rtDW . awqg33it32 ) , sizeof (
rtDW . awqg33it32 ) ) ; mxSetFieldByNumber ( ssDW , 0 , 1 , rtdwData ) ; }
mr_Copy_of_trial_with_intermediate_tanks_cacheDataAsMxArray ( ssDW , 0 , 2 ,
( const void * ) & ( rtPrevZCX ) , sizeof ( rtPrevZCX ) ) ; return ssDW ; }
void mr_Copy_of_trial_with_intermediate_tanks_SetDWork ( const mxArray * ssDW
) { ( void ) ssDW ;
mr_Copy_of_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * )
& ( rtB ) , ssDW , 0 , 0 , sizeof ( rtB ) ) ; { const mxArray * rtdwData =
mxGetFieldByNumber ( ssDW , 0 , 1 ) ;
mr_Copy_of_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * )
& ( rtDW . j0cswhuhs0 ) , rtdwData , 0 , 0 , sizeof ( rtDW . j0cswhuhs0 ) ) ;
mr_Copy_of_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * )
& ( rtDW . atxpvhc5n0 ) , rtdwData , 0 , 1 , sizeof ( rtDW . atxpvhc5n0 ) ) ;
mr_Copy_of_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * )
& ( rtDW . inlqdh4kh1 ) , rtdwData , 0 , 2 , sizeof ( rtDW . inlqdh4kh1 ) ) ;
mr_Copy_of_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * )
& ( rtDW . glwrm0q1go ) , rtdwData , 0 , 3 , sizeof ( rtDW . glwrm0q1go ) ) ;
mr_Copy_of_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * )
& ( rtDW . apkag20puf ) , rtdwData , 0 , 4 , sizeof ( rtDW . apkag20puf ) ) ;
mr_Copy_of_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * )
& ( rtDW . mpaqmqlwie ) , rtdwData , 0 , 5 , sizeof ( rtDW . mpaqmqlwie ) ) ;
mr_Copy_of_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * )
& ( rtDW . a4pqll03x1 ) , rtdwData , 0 , 6 , sizeof ( rtDW . a4pqll03x1 ) ) ;
mr_Copy_of_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * )
& ( rtDW . p1ctkshy11 ) , rtdwData , 0 , 7 , sizeof ( rtDW . p1ctkshy11 ) ) ;
mr_Copy_of_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * )
& ( rtDW . d34s0zbduf ) , rtdwData , 0 , 8 , sizeof ( rtDW . d34s0zbduf ) ) ;
mr_Copy_of_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * )
& ( rtDW . epcaebihmu ) , rtdwData , 0 , 9 , sizeof ( rtDW . epcaebihmu ) ) ;
mr_Copy_of_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * )
& ( rtDW . kvbq5jwbvw ) , rtdwData , 0 , 10 , sizeof ( rtDW . kvbq5jwbvw ) )
; mr_Copy_of_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void *
) & ( rtDW . audz3xgmq3 ) , rtdwData , 0 , 11 , sizeof ( rtDW . audz3xgmq3 )
) ; mr_Copy_of_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void
* ) & ( rtDW . loaiy2cacx ) , rtdwData , 0 , 12 , sizeof ( rtDW . loaiy2cacx
) ) ; mr_Copy_of_trial_with_intermediate_tanks_restoreDataFromMxArray ( (
void * ) & ( rtDW . e0qqpzfruh ) , rtdwData , 0 , 13 , sizeof ( rtDW .
e0qqpzfruh ) ) ;
mr_Copy_of_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * )
& ( rtDW . hamy5qs5qe ) , rtdwData , 0 , 14 , sizeof ( rtDW . hamy5qs5qe ) )
; mr_Copy_of_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void *
) & ( rtDW . je5l4cyavv ) , rtdwData , 0 , 15 , sizeof ( rtDW . je5l4cyavv )
) ; mr_Copy_of_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void
* ) & ( rtDW . ebwtwofxpd ) , rtdwData , 0 , 16 , sizeof ( rtDW . ebwtwofxpd
) ) ; mr_Copy_of_trial_with_intermediate_tanks_restoreDataFromMxArray ( (
void * ) & ( rtDW . aku4xs4fdc ) , rtdwData , 0 , 17 , sizeof ( rtDW .
aku4xs4fdc ) ) ;
mr_Copy_of_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * )
& ( rtDW . aigzt14viy ) , rtdwData , 0 , 18 , sizeof ( rtDW . aigzt14viy ) )
; mr_Copy_of_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void *
) & ( rtDW . jwuiy2atz3 ) , rtdwData , 0 , 19 , sizeof ( rtDW . jwuiy2atz3 )
) ; mr_Copy_of_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void
* ) & ( rtDW . isc2qgvxqr ) , rtdwData , 0 , 20 , sizeof ( rtDW . isc2qgvxqr
) ) ; mr_Copy_of_trial_with_intermediate_tanks_restoreDataFromMxArray ( (
void * ) & ( rtDW . awwuosv0rp ) , rtdwData , 0 , 21 , sizeof ( rtDW .
awwuosv0rp ) ) ;
mr_Copy_of_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * )
& ( rtDW . g40nougyde ) , rtdwData , 0 , 22 , sizeof ( rtDW . g40nougyde ) )
; mr_Copy_of_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void *
) & ( rtDW . cangmlit2r ) , rtdwData , 0 , 23 , sizeof ( rtDW . cangmlit2r )
) ; mr_Copy_of_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void
* ) & ( rtDW . od3s0wtt1r ) , rtdwData , 0 , 24 , sizeof ( rtDW . od3s0wtt1r
) ) ; mr_Copy_of_trial_with_intermediate_tanks_restoreDataFromMxArray ( (
void * ) & ( rtDW . ppnlaclhj1 ) , rtdwData , 0 , 25 , sizeof ( rtDW .
ppnlaclhj1 ) ) ;
mr_Copy_of_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * )
& ( rtDW . exujqdty5l ) , rtdwData , 0 , 26 , sizeof ( rtDW . exujqdty5l ) )
; mr_Copy_of_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void *
) & ( rtDW . iln55qrgxk ) , rtdwData , 0 , 27 , sizeof ( rtDW . iln55qrgxk )
) ; mr_Copy_of_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void
* ) & ( rtDW . i213auut4h ) , rtdwData , 0 , 28 , sizeof ( rtDW . i213auut4h
) ) ; mr_Copy_of_trial_with_intermediate_tanks_restoreDataFromMxArray ( (
void * ) & ( rtDW . e13oplclke ) , rtdwData , 0 , 29 , sizeof ( rtDW .
e13oplclke ) ) ;
mr_Copy_of_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * )
& ( rtDW . laquyfl1kp ) , rtdwData , 0 , 30 , sizeof ( rtDW . laquyfl1kp ) )
; mr_Copy_of_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void *
) & ( rtDW . f45mx2s03m ) , rtdwData , 0 , 31 , sizeof ( rtDW . f45mx2s03m )
) ; mr_Copy_of_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void
* ) & ( rtDW . gstdmwepbc ) , rtdwData , 0 , 32 , sizeof ( rtDW . gstdmwepbc
) ) ; mr_Copy_of_trial_with_intermediate_tanks_restoreDataFromMxArray ( (
void * ) & ( rtDW . ikvzyivi0j ) , rtdwData , 0 , 33 , sizeof ( rtDW .
ikvzyivi0j ) ) ;
mr_Copy_of_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * )
& ( rtDW . libi2dsvky ) , rtdwData , 0 , 34 , sizeof ( rtDW . libi2dsvky ) )
; mr_Copy_of_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void *
) & ( rtDW . ga2xty0o4o ) , rtdwData , 0 , 35 , sizeof ( rtDW . ga2xty0o4o )
) ; mr_Copy_of_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void
* ) & ( rtDW . e1b1tx1d4l ) , rtdwData , 0 , 36 , sizeof ( rtDW . e1b1tx1d4l
) ) ; mr_Copy_of_trial_with_intermediate_tanks_restoreDataFromMxArray ( (
void * ) & ( rtDW . p0eq4zkt3w ) , rtdwData , 0 , 37 , sizeof ( rtDW .
p0eq4zkt3w ) ) ;
mr_Copy_of_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * )
& ( rtDW . oahmxad2ft ) , rtdwData , 0 , 38 , sizeof ( rtDW . oahmxad2ft ) )
; mr_Copy_of_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void *
) & ( rtDW . ksvy5jh2vc ) , rtdwData , 0 , 39 , sizeof ( rtDW . ksvy5jh2vc )
) ; mr_Copy_of_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void
* ) & ( rtDW . fyaulsukmr ) , rtdwData , 0 , 40 , sizeof ( rtDW . fyaulsukmr
) ) ; mr_Copy_of_trial_with_intermediate_tanks_restoreDataFromMxArray ( (
void * ) & ( rtDW . hdus1cnzep ) , rtdwData , 0 , 41 , sizeof ( rtDW .
hdus1cnzep ) ) ;
mr_Copy_of_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * )
& ( rtDW . amiin4kdmj ) , rtdwData , 0 , 42 , sizeof ( rtDW . amiin4kdmj ) )
; mr_Copy_of_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void *
) & ( rtDW . nwxy1ani3r ) , rtdwData , 0 , 43 , sizeof ( rtDW . nwxy1ani3r )
) ; mr_Copy_of_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void
* ) & ( rtDW . j5zzepflri ) , rtdwData , 0 , 44 , sizeof ( rtDW . j5zzepflri
) ) ; mr_Copy_of_trial_with_intermediate_tanks_restoreDataFromMxArray ( (
void * ) & ( rtDW . hk25kkxqwm ) , rtdwData , 0 , 45 , sizeof ( rtDW .
hk25kkxqwm ) ) ;
mr_Copy_of_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * )
& ( rtDW . ffd0olavlx ) , rtdwData , 0 , 46 , sizeof ( rtDW . ffd0olavlx ) )
; mr_Copy_of_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void *
) & ( rtDW . eiohhcwuqa ) , rtdwData , 0 , 47 , sizeof ( rtDW . eiohhcwuqa )
) ; mr_Copy_of_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void
* ) & ( rtDW . eyyp3zs3le ) , rtdwData , 0 , 48 , sizeof ( rtDW . eyyp3zs3le
) ) ; mr_Copy_of_trial_with_intermediate_tanks_restoreDataFromMxArray ( (
void * ) & ( rtDW . nakeaw2rna ) , rtdwData , 0 , 49 , sizeof ( rtDW .
nakeaw2rna ) ) ;
mr_Copy_of_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * )
& ( rtDW . n4ne2vzdw3 ) , rtdwData , 0 , 50 , sizeof ( rtDW . n4ne2vzdw3 ) )
; mr_Copy_of_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void *
) & ( rtDW . ljfwvqlys1 ) , rtdwData , 0 , 51 , sizeof ( rtDW . ljfwvqlys1 )
) ; mr_Copy_of_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void
* ) & ( rtDW . ct4cduudvv ) , rtdwData , 0 , 52 , sizeof ( rtDW . ct4cduudvv
) ) ; mr_Copy_of_trial_with_intermediate_tanks_restoreDataFromMxArray ( (
void * ) & ( rtDW . jkvzly4a1o ) , rtdwData , 0 , 53 , sizeof ( rtDW .
jkvzly4a1o ) ) ;
mr_Copy_of_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * )
& ( rtDW . fxawe3htlm ) , rtdwData , 0 , 54 , sizeof ( rtDW . fxawe3htlm ) )
; mr_Copy_of_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void *
) & ( rtDW . irq32iym0t ) , rtdwData , 0 , 55 , sizeof ( rtDW . irq32iym0t )
) ; mr_Copy_of_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void
* ) & ( rtDW . jz3v2t25p2 ) , rtdwData , 0 , 56 , sizeof ( rtDW . jz3v2t25p2
) ) ; mr_Copy_of_trial_with_intermediate_tanks_restoreDataFromMxArray ( (
void * ) & ( rtDW . m0i52xetkk ) , rtdwData , 0 , 57 , sizeof ( rtDW .
m0i52xetkk ) ) ;
mr_Copy_of_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * )
& ( rtDW . dfvcjzwy4h ) , rtdwData , 0 , 58 , sizeof ( rtDW . dfvcjzwy4h ) )
; mr_Copy_of_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void *
) & ( rtDW . nas3vtbw5t ) , rtdwData , 0 , 59 , sizeof ( rtDW . nas3vtbw5t )
) ; mr_Copy_of_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void
* ) & ( rtDW . awqg33it32 ) , rtdwData , 0 , 60 , sizeof ( rtDW . awqg33it32
) ) ; } mr_Copy_of_trial_with_intermediate_tanks_restoreDataFromMxArray ( (
void * ) & ( rtPrevZCX ) , ssDW , 0 , 2 , sizeof ( rtPrevZCX ) ) ; } mxArray
* mr_Copy_of_trial_with_intermediate_tanks_GetSimStateDisallowedBlocks ( ) {
mxArray * data = mxCreateCellMatrix ( 8 , 3 ) ; mwIndex subs [ 2 ] , offset ;
{ static const char_T * blockType [ 8 ] = { "Scope" , "Scope" , "Scope" ,
"Scope" , "Scope" , "Scope" , "Scope" , "Scope" , } ; static const char_T *
blockPath [ 8 ] = { "Copy_of_trial_with_intermediate_tanks/Energy Comparison"
, "Copy_of_trial_with_intermediate_tanks/Scope" ,
"Copy_of_trial_with_intermediate_tanks/Scope1" ,
"Copy_of_trial_with_intermediate_tanks/Scope2" ,
"Copy_of_trial_with_intermediate_tanks/Scope3" ,
"Copy_of_trial_with_intermediate_tanks/Scope4" ,
"Copy_of_trial_with_intermediate_tanks/Scope5" ,
"Copy_of_trial_with_intermediate_tanks/Scope6" , } ; static const int reason
[ 8 ] = { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , } ; for ( subs [ 0 ] = 0 ; subs [ 0
] < 8 ; ++ ( subs [ 0 ] ) ) { subs [ 1 ] = 0 ; offset = mxCalcSingleSubscript
( data , 2 , subs ) ; mxSetCell ( data , offset , mxCreateString ( blockType
[ subs [ 0 ] ] ) ) ; subs [ 1 ] = 1 ; offset = mxCalcSingleSubscript ( data ,
2 , subs ) ; mxSetCell ( data , offset , mxCreateString ( blockPath [ subs [
0 ] ] ) ) ; subs [ 1 ] = 2 ; offset = mxCalcSingleSubscript ( data , 2 , subs
) ; mxSetCell ( data , offset , mxCreateDoubleScalar ( ( real_T ) reason [
subs [ 0 ] ] ) ) ; } } return data ; } void MdlInitializeSizes ( void ) {
ssSetNumContStates ( rtS , 27 ) ; ssSetNumPeriodicContStates ( rtS , 0 ) ;
ssSetNumY ( rtS , 0 ) ; ssSetNumU ( rtS , 0 ) ; ssSetDirectFeedThrough ( rtS
, 0 ) ; ssSetNumSampleTimes ( rtS , 3 ) ; ssSetNumBlocks ( rtS , 146 ) ;
ssSetNumBlockIO ( rtS , 79 ) ; ssSetNumBlockParams ( rtS , 56 ) ; } void
MdlInitializeSampleTimes ( void ) { ssSetSampleTime ( rtS , 0 , 0.0 ) ;
ssSetSampleTime ( rtS , 1 , 0.0 ) ; ssSetSampleTime ( rtS , 2 , 1.0E-6 ) ;
ssSetOffsetTime ( rtS , 0 , 0.0 ) ; ssSetOffsetTime ( rtS , 1 , 1.0 ) ;
ssSetOffsetTime ( rtS , 2 , 0.0 ) ; } void raccel_set_checksum ( ) {
ssSetChecksumVal ( rtS , 0 , 2453326653U ) ; ssSetChecksumVal ( rtS , 1 ,
3278769589U ) ; ssSetChecksumVal ( rtS , 2 , 1795305681U ) ; ssSetChecksumVal
( rtS , 3 , 4235931226U ) ; }
#if defined(_MSC_VER)
#pragma optimize( "", off )
#endif
SimStruct * raccel_register_model ( ssExecutionInfo * executionInfo ) {
static struct _ssMdlInfo mdlInfo ; static struct _ssBlkInfo2 blkInfo2 ;
static struct _ssBlkInfoSLSize blkInfoSLSize ; rt_modelMapInfoPtr = & (
rt_dataMapInfo . mmi ) ; executionInfo -> gblObjects_ . numToFiles = 0 ;
executionInfo -> gblObjects_ . numFrFiles = 0 ; executionInfo -> gblObjects_
. numFrWksBlocks = 0 ; executionInfo -> gblObjects_ . numModelInputs = 0 ;
executionInfo -> gblObjects_ . numRootInportBlks = 0 ; executionInfo ->
gblObjects_ . inportDataTypeIdx = NULL ; executionInfo -> gblObjects_ .
inportDims = NULL ; executionInfo -> gblObjects_ . inportComplex = NULL ;
executionInfo -> gblObjects_ . inportInterpoFlag = NULL ; executionInfo ->
gblObjects_ . inportContinuous = NULL ; ( void ) memset ( ( char_T * ) rtS ,
0 , sizeof ( SimStruct ) ) ; ( void ) memset ( ( char_T * ) & mdlInfo , 0 ,
sizeof ( struct _ssMdlInfo ) ) ; ( void ) memset ( ( char_T * ) & blkInfo2 ,
0 , sizeof ( struct _ssBlkInfo2 ) ) ; ( void ) memset ( ( char_T * ) &
blkInfoSLSize , 0 , sizeof ( struct _ssBlkInfoSLSize ) ) ; ssSetBlkInfo2Ptr (
rtS , & blkInfo2 ) ; ssSetBlkInfoSLSizePtr ( rtS , & blkInfoSLSize ) ;
ssSetMdlInfoPtr ( rtS , & mdlInfo ) ; ssSetExecutionInfo ( rtS ,
executionInfo ) ; slsaAllocOPModelData ( rtS ) ; { static time_T mdlPeriod [
NSAMPLE_TIMES ] ; static time_T mdlOffset [ NSAMPLE_TIMES ] ; static time_T
mdlTaskTimes [ NSAMPLE_TIMES ] ; static int_T mdlTsMap [ NSAMPLE_TIMES ] ;
static int_T mdlSampleHits [ NSAMPLE_TIMES ] ; static boolean_T
mdlTNextWasAdjustedPtr [ NSAMPLE_TIMES ] ; static int_T mdlPerTaskSampleHits
[ NSAMPLE_TIMES * NSAMPLE_TIMES ] ; static time_T mdlTimeOfNextSampleHit [
NSAMPLE_TIMES ] ; { int_T i ; for ( i = 0 ; i < NSAMPLE_TIMES ; i ++ ) {
mdlPeriod [ i ] = 0.0 ; mdlOffset [ i ] = 0.0 ; mdlTaskTimes [ i ] = 0.0 ;
mdlTsMap [ i ] = i ; mdlSampleHits [ i ] = 1 ; } } ssSetSampleTimePtr ( rtS ,
& mdlPeriod [ 0 ] ) ; ssSetOffsetTimePtr ( rtS , & mdlOffset [ 0 ] ) ;
ssSetSampleTimeTaskIDPtr ( rtS , & mdlTsMap [ 0 ] ) ; ssSetTPtr ( rtS , &
mdlTaskTimes [ 0 ] ) ; ssSetSampleHitPtr ( rtS , & mdlSampleHits [ 0 ] ) ;
ssSetTNextWasAdjustedPtr ( rtS , & mdlTNextWasAdjustedPtr [ 0 ] ) ;
ssSetPerTaskSampleHitsPtr ( rtS , & mdlPerTaskSampleHits [ 0 ] ) ;
ssSetTimeOfNextSampleHitPtr ( rtS , & mdlTimeOfNextSampleHit [ 0 ] ) ; }
ssSetSolverMode ( rtS , SOLVER_MODE_SINGLETASKING ) ; { ssSetBlockIO ( rtS ,
( ( void * ) & rtB ) ) ; ( void ) memset ( ( ( void * ) & rtB ) , 0 , sizeof
( B ) ) ; } { real_T * x = ( real_T * ) & rtX ; ssSetContStates ( rtS , x ) ;
( void ) memset ( ( void * ) x , 0 , sizeof ( X ) ) ; } { void * dwork = (
void * ) & rtDW ; ssSetRootDWork ( rtS , dwork ) ; ( void ) memset ( dwork ,
0 , sizeof ( DW ) ) ; } { static DataTypeTransInfo dtInfo ; ( void ) memset (
( char_T * ) & dtInfo , 0 , sizeof ( dtInfo ) ) ; ssSetModelMappingInfo ( rtS
, & dtInfo ) ; dtInfo . numDataTypes = 24 ; dtInfo . dataTypeSizes = &
rtDataTypeSizes [ 0 ] ; dtInfo . dataTypeNames = & rtDataTypeNames [ 0 ] ;
dtInfo . BTransTable = & rtBTransTable ; dtInfo . PTransTable = &
rtPTransTable ; dtInfo . dataTypeInfoTable = rtDataTypeInfoTable ; }
Copy_of_trial_with_intermediate_tanks_InitializeDataMapInfo ( ) ;
ssSetIsRapidAcceleratorActive ( rtS , true ) ; ssSetRootSS ( rtS , rtS ) ;
ssSetVersion ( rtS , SIMSTRUCT_VERSION_LEVEL2 ) ; ssSetModelName ( rtS ,
"Copy_of_trial_with_intermediate_tanks" ) ; ssSetPath ( rtS ,
"Copy_of_trial_with_intermediate_tanks" ) ; ssSetTStart ( rtS , 0.0 ) ;
ssSetTFinal ( rtS , 0.2 ) ; { static RTWLogInfo rt_DataLoggingInfo ;
rt_DataLoggingInfo . loggingInterval = ( NULL ) ; ssSetRTWLogInfo ( rtS , &
rt_DataLoggingInfo ) ; } { { static int_T rt_LoggedStateWidths [ ] = { 1 , 1
, 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 2 , 2 , 2 , 2 , 1 , 1 ,
2 , 1 } ; static int_T rt_LoggedStateNumDimensions [ ] = { 1 , 1 , 1 , 1 , 1
, 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 } ;
static int_T rt_LoggedStateDimensions [ ] = { 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 ,
1 , 1 , 1 , 1 , 1 , 1 , 1 , 2 , 2 , 2 , 2 , 1 , 1 , 2 , 1 } ; static
boolean_T rt_LoggedStateIsVarDims [ ] = { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ,
0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 } ; static
BuiltInDTypeId rt_LoggedStateDataTypeIds [ ] = { SS_DOUBLE , SS_DOUBLE ,
SS_DOUBLE , SS_DOUBLE , SS_DOUBLE , SS_DOUBLE , SS_DOUBLE , SS_DOUBLE ,
SS_DOUBLE , SS_DOUBLE , SS_DOUBLE , SS_DOUBLE , SS_DOUBLE , SS_DOUBLE ,
SS_DOUBLE , SS_DOUBLE , SS_DOUBLE , SS_DOUBLE , SS_DOUBLE , SS_DOUBLE ,
SS_DOUBLE , SS_DOUBLE , SS_DOUBLE } ; static int_T
rt_LoggedStateComplexSignals [ ] = { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ,
0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 } ; static
RTWPreprocessingFcnPtr rt_LoggingStatePreprocessingFcnPtrs [ ] = { ( NULL ) ,
( NULL ) , ( NULL ) , ( NULL ) , ( NULL ) , ( NULL ) , ( NULL ) , ( NULL ) ,
( NULL ) , ( NULL ) , ( NULL ) , ( NULL ) , ( NULL ) , ( NULL ) , ( NULL ) ,
( NULL ) , ( NULL ) , ( NULL ) , ( NULL ) , ( NULL ) , ( NULL ) , ( NULL ) ,
( NULL ) } ; static const char_T * rt_LoggedStateLabels [ ] = { "CSTATE" ,
"CSTATE" , "CSTATE" , "CSTATE" , "CSTATE" , "CSTATE" , "CSTATE" , "CSTATE" ,
"CSTATE" , "CSTATE" , "CSTATE" , "CSTATE" , "CSTATE" , "CSTATE" , "CSTATE" ,
"CSTATE" , "CSTATE" , "CSTATE" , "CSTATE" , "CSTATE" , "CSTATE" , "CSTATE" ,
"DSTATE" } ; static const char_T * rt_LoggedStateBlockNames [ ] = {
"Copy_of_trial_with_intermediate_tanks/Integrator4" ,
"Copy_of_trial_with_intermediate_tanks/Integrator2" ,
"Copy_of_trial_with_intermediate_tanks/Integrator13" ,
"Copy_of_trial_with_intermediate_tanks/Integrator" ,
"Copy_of_trial_with_intermediate_tanks/Integrator1" ,
"Copy_of_trial_with_intermediate_tanks/Integrator14" ,
"Copy_of_trial_with_intermediate_tanks/Integrator7" ,
"Copy_of_trial_with_intermediate_tanks/Integrator12" ,
"Copy_of_trial_with_intermediate_tanks/Integrator10" ,
"Copy_of_trial_with_intermediate_tanks/Integrator8" ,
"Copy_of_trial_with_intermediate_tanks/Integrator9" ,
"Copy_of_trial_with_intermediate_tanks/Integrator11" ,
"Copy_of_trial_with_intermediate_tanks/Integrator3" ,
"Copy_of_trial_with_intermediate_tanks/Integrator6" ,
"Copy_of_trial_with_intermediate_tanks/Integrator16" ,
"Copy_of_trial_with_intermediate_tanks/Dynamics of Valve1/Transfer Fcn" ,
"Copy_of_trial_with_intermediate_tanks/Dynamics of Valve2/Transfer Fcn" ,
"Copy_of_trial_with_intermediate_tanks/Dynamics of Valve3/Transfer Fcn" ,
"Copy_of_trial_with_intermediate_tanks/Dynamics of Valve4/Transfer Fcn" ,
"Copy_of_trial_with_intermediate_tanks/Integrator5" ,
"Copy_of_trial_with_intermediate_tanks/Integrator15" ,
"Copy_of_trial_with_intermediate_tanks/Dynamics of Valve/Transfer Fcn" ,
"Copy_of_trial_with_intermediate_tanks/Delay" } ; static const char_T *
rt_LoggedStateNames [ ] = { "" , "" , "" , "" , "" , "" , "" , "" , "" , "" ,
"" , "" , "" , "" , "" , "" , "" , "" , "" , "" , "" , "" , "DSTATE" } ;
static boolean_T rt_LoggedStateCrossMdlRef [ ] = { 0 , 0 , 0 , 0 , 0 , 0 , 0
, 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 } ; static
RTWLogDataTypeConvert rt_RTWLogDataTypeConvert [ ] = { { 0 , SS_DOUBLE ,
SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 , 0.0 } , { 0 , SS_DOUBLE , SS_DOUBLE , 0 , 0
, 0 , 1.0 , 0 , 0.0 } , { 0 , SS_DOUBLE , SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 ,
0.0 } , { 0 , SS_DOUBLE , SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 , 0.0 } , { 0 ,
SS_DOUBLE , SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 , 0.0 } , { 0 , SS_DOUBLE ,
SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 , 0.0 } , { 0 , SS_DOUBLE , SS_DOUBLE , 0 , 0
, 0 , 1.0 , 0 , 0.0 } , { 0 , SS_DOUBLE , SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 ,
0.0 } , { 0 , SS_DOUBLE , SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 , 0.0 } , { 0 ,
SS_DOUBLE , SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 , 0.0 } , { 0 , SS_DOUBLE ,
SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 , 0.0 } , { 0 , SS_DOUBLE , SS_DOUBLE , 0 , 0
, 0 , 1.0 , 0 , 0.0 } , { 0 , SS_DOUBLE , SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 ,
0.0 } , { 0 , SS_DOUBLE , SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 , 0.0 } , { 0 ,
SS_DOUBLE , SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 , 0.0 } , { 0 , SS_DOUBLE ,
SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 , 0.0 } , { 0 , SS_DOUBLE , SS_DOUBLE , 0 , 0
, 0 , 1.0 , 0 , 0.0 } , { 0 , SS_DOUBLE , SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 ,
0.0 } , { 0 , SS_DOUBLE , SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 , 0.0 } , { 0 ,
SS_DOUBLE , SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 , 0.0 } , { 0 , SS_DOUBLE ,
SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 , 0.0 } , { 0 , SS_DOUBLE , SS_DOUBLE , 0 , 0
, 0 , 1.0 , 0 , 0.0 } , { 0 , SS_DOUBLE , SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 ,
0.0 } } ; static int_T rt_LoggedStateIdxList [ ] = { 0 , 1 , 2 , 3 , 4 , 5 ,
6 , 7 , 8 , 9 , 10 , 11 , 12 , 13 , 14 , 15 , 16 , 17 , 18 , 19 , 20 , 21 , 0
} ; static RTWLogSignalInfo rt_LoggedStateSignalInfo = { 23 ,
rt_LoggedStateWidths , rt_LoggedStateNumDimensions , rt_LoggedStateDimensions
, rt_LoggedStateIsVarDims , ( NULL ) , ( NULL ) , rt_LoggedStateDataTypeIds ,
rt_LoggedStateComplexSignals , ( NULL ) , rt_LoggingStatePreprocessingFcnPtrs
, { rt_LoggedStateLabels } , ( NULL ) , ( NULL ) , ( NULL ) , {
rt_LoggedStateBlockNames } , { rt_LoggedStateNames } ,
rt_LoggedStateCrossMdlRef , rt_RTWLogDataTypeConvert , rt_LoggedStateIdxList
} ; static void * rt_LoggedStateSignalPtrs [ 23 ] ; rtliSetLogXSignalPtrs (
ssGetRTWLogInfo ( rtS ) , ( LogSignalPtrsType ) rt_LoggedStateSignalPtrs ) ;
rtliSetLogXSignalInfo ( ssGetRTWLogInfo ( rtS ) , & rt_LoggedStateSignalInfo
) ; rt_LoggedStateSignalPtrs [ 0 ] = ( void * ) & rtX . eer5ulgpq1 ;
rt_LoggedStateSignalPtrs [ 1 ] = ( void * ) & rtX . gvqywnxquo ;
rt_LoggedStateSignalPtrs [ 2 ] = ( void * ) & rtX . mui3hk2kdn ;
rt_LoggedStateSignalPtrs [ 3 ] = ( void * ) & rtX . kllv1nez4n ;
rt_LoggedStateSignalPtrs [ 4 ] = ( void * ) & rtX . orxgnzk15p ;
rt_LoggedStateSignalPtrs [ 5 ] = ( void * ) & rtX . apttuwn22n ;
rt_LoggedStateSignalPtrs [ 6 ] = ( void * ) & rtX . avebadcay4 ;
rt_LoggedStateSignalPtrs [ 7 ] = ( void * ) & rtX . o2kfetkxay ;
rt_LoggedStateSignalPtrs [ 8 ] = ( void * ) & rtX . fvfvpkek2n ;
rt_LoggedStateSignalPtrs [ 9 ] = ( void * ) & rtX . aqqzk12fdv ;
rt_LoggedStateSignalPtrs [ 10 ] = ( void * ) & rtX . e30fqjkckb ;
rt_LoggedStateSignalPtrs [ 11 ] = ( void * ) & rtX . mizgi5ztop ;
rt_LoggedStateSignalPtrs [ 12 ] = ( void * ) & rtX . iikwtjyqxf ;
rt_LoggedStateSignalPtrs [ 13 ] = ( void * ) & rtX . c5bbd2cwaf ;
rt_LoggedStateSignalPtrs [ 14 ] = ( void * ) & rtX . loqzkyynb5 ;
rt_LoggedStateSignalPtrs [ 15 ] = ( void * ) & rtX . klye04tpmo [ 0 ] ;
rt_LoggedStateSignalPtrs [ 16 ] = ( void * ) & rtX . ht5behjiqe [ 0 ] ;
rt_LoggedStateSignalPtrs [ 17 ] = ( void * ) & rtX . aa2dkshjtd [ 0 ] ;
rt_LoggedStateSignalPtrs [ 18 ] = ( void * ) & rtX . kco1vq4cxg [ 0 ] ;
rt_LoggedStateSignalPtrs [ 19 ] = ( void * ) & rtX . fklbo31ygf ;
rt_LoggedStateSignalPtrs [ 20 ] = ( void * ) & rtX . kwoozlzy2w ;
rt_LoggedStateSignalPtrs [ 21 ] = ( void * ) & rtX . mvzjw32hje [ 0 ] ;
rt_LoggedStateSignalPtrs [ 22 ] = ( void * ) & rtDW . j0cswhuhs0 ; }
rtliSetLogT ( ssGetRTWLogInfo ( rtS ) , "tout" ) ; rtliSetLogX (
ssGetRTWLogInfo ( rtS ) , "" ) ; rtliSetLogXFinal ( ssGetRTWLogInfo ( rtS ) ,
"xFinal" ) ; rtliSetLogVarNameModifier ( ssGetRTWLogInfo ( rtS ) , "none" ) ;
rtliSetLogFormat ( ssGetRTWLogInfo ( rtS ) , 4 ) ; rtliSetLogMaxRows (
ssGetRTWLogInfo ( rtS ) , 0 ) ; rtliSetLogDecimation ( ssGetRTWLogInfo ( rtS
) , 1 ) ; rtliSetLogY ( ssGetRTWLogInfo ( rtS ) , "" ) ;
rtliSetLogYSignalInfo ( ssGetRTWLogInfo ( rtS ) , ( NULL ) ) ;
rtliSetLogYSignalPtrs ( ssGetRTWLogInfo ( rtS ) , ( NULL ) ) ; } { static
struct _ssStatesInfo2 statesInfo2 ; ssSetStatesInfo2 ( rtS , & statesInfo2 )
; } { static ssPeriodicStatesInfo periodicStatesInfo ;
ssSetPeriodicStatesInfo ( rtS , & periodicStatesInfo ) ; } { static
ssJacobianPerturbationBounds jacobianPerturbationBounds ;
ssSetJacobianPerturbationBounds ( rtS , & jacobianPerturbationBounds ) ; } {
static ssSolverInfo slvrInfo ; static boolean_T contStatesDisabled [ 27 ] ;
static real_T absTol [ 27 ] = { 1.0000000000000001E-11 ,
1.0000000000000001E-11 , 1.0000000000000001E-11 , 1.0000000000000001E-11 ,
1.0000000000000001E-11 , 1.0000000000000001E-11 , 1.0000000000000001E-11 ,
1.0000000000000001E-11 , 1.0000000000000001E-11 , 1.0000000000000001E-11 ,
1.0000000000000001E-11 , 1.0000000000000001E-11 , 1.0000000000000001E-11 ,
1.0000000000000001E-11 , 1.0000000000000001E-11 , 1.0000000000000001E-11 ,
1.0000000000000001E-11 , 1.0000000000000001E-11 , 1.0000000000000001E-11 ,
1.0000000000000001E-11 , 1.0000000000000001E-11 , 1.0000000000000001E-11 ,
1.0000000000000001E-11 , 1.0000000000000001E-11 , 1.0000000000000001E-11 ,
1.0000000000000001E-11 , 1.0000000000000001E-11 } ; static uint8_T
absTolControl [ 27 ] = { 0U , 0U , 0U , 0U , 0U , 0U , 0U , 0U , 0U , 0U , 0U
, 0U , 0U , 0U , 0U , 0U , 0U , 0U , 0U , 0U , 0U , 0U , 0U , 0U , 0U , 0U ,
0U } ; static real_T contStateJacPerturbBoundMinVec [ 27 ] ; static real_T
contStateJacPerturbBoundMaxVec [ 27 ] ; static uint8_T zcAttributes [ 13 ] =
{ ( ZC_EVENT_ALL ) , ( ZC_EVENT_ALL ) , ( ZC_EVENT_ALL ) , ( ZC_EVENT_ALL ) ,
( ZC_EVENT_ALL ) , ( ZC_EVENT_ALL ) , ( ZC_EVENT_ALL ) , ( ZC_EVENT_ALL ) , (
ZC_EVENT_ALL ) , ( ZC_EVENT_ALL ) , ( ZC_EVENT_ALL ) , ( ZC_EVENT_ALL ) , (
0xc0 | ZC_EVENT_ALL_UP ) } ; static ssNonContDerivSigInfo nonContDerivSigInfo
[ 6 ] = { { 1 * sizeof ( real_T ) , ( char * ) ( & rtB . iiwlxh5nit ) , (
NULL ) } , { 1 * sizeof ( real_T ) , ( char * ) ( & rtB . lr3evndp05 ) , (
NULL ) } , { 1 * sizeof ( real_T ) , ( char * ) ( & rtB . pju2ibjy53 ) , (
NULL ) } , { 1 * sizeof ( real_T ) , ( char * ) ( & rtB . jpwjwhvvnr ) , (
NULL ) } , { 1 * sizeof ( real_T ) , ( char * ) ( & rtB . owtjgmweb5 ) , (
NULL ) } , { 1 * sizeof ( real_T ) , ( char * ) ( & rtB . hmcqeqif2o ) , (
NULL ) } } ; { int i ; for ( i = 0 ; i < 27 ; ++ i ) {
contStateJacPerturbBoundMinVec [ i ] = 0 ; contStateJacPerturbBoundMaxVec [ i
] = rtGetInf ( ) ; } } ssSetSolverRelTol ( rtS , 1.0E-8 ) ; ssSetStepSize (
rtS , 0.0 ) ; ssSetMinStepSize ( rtS , 0.0 ) ; ssSetMaxNumMinSteps ( rtS , -
1 ) ; ssSetMinStepViolatedError ( rtS , 0 ) ; ssSetMaxStepSize ( rtS , 1.0E-6
) ; ssSetSolverMaxOrder ( rtS , - 1 ) ; ssSetSolverRefineFactor ( rtS , 1 ) ;
ssSetOutputTimes ( rtS , ( NULL ) ) ; ssSetNumOutputTimes ( rtS , 0 ) ;
ssSetOutputTimesOnly ( rtS , 0 ) ; ssSetOutputTimesIndex ( rtS , 0 ) ;
ssSetZCCacheNeedsReset ( rtS , 0 ) ; ssSetDerivCacheNeedsReset ( rtS , 0 ) ;
ssSetNumNonContDerivSigInfos ( rtS , 6 ) ; ssSetNonContDerivSigInfos ( rtS ,
nonContDerivSigInfo ) ; ssSetSolverInfo ( rtS , & slvrInfo ) ;
ssSetSolverName ( rtS , "ode23s" ) ; ssSetVariableStepSolver ( rtS , 1 ) ;
ssSetSolverConsistencyChecking ( rtS , 0 ) ; ssSetSolverAdaptiveZcDetection (
rtS , 0 ) ; ssSetSolverRobustResetMethod ( rtS , 0 ) ; ssSetAbsTolVector (
rtS , absTol ) ; ssSetAbsTolControlVector ( rtS , absTolControl ) ;
ssSetSolverAbsTol_Obsolete ( rtS , absTol ) ;
ssSetSolverAbsTolControl_Obsolete ( rtS , absTolControl ) ;
ssSetJacobianPerturbationBoundsMinVec ( rtS , contStateJacPerturbBoundMinVec
) ; ssSetJacobianPerturbationBoundsMaxVec ( rtS ,
contStateJacPerturbBoundMaxVec ) ; ssSetSolverStateProjection ( rtS , 0 ) ;
ssSetSolverMassMatrixType ( rtS , ( ssMatrixType ) 0 ) ;
ssSetSolverMassMatrixNzMax ( rtS , 0 ) ; ssSetModelOutputs ( rtS , MdlOutputs
) ; ssSetModelUpdate ( rtS , MdlUpdate ) ; ssSetModelDerivatives ( rtS ,
MdlDerivatives ) ; ssSetSolverZcSignalAttrib ( rtS , zcAttributes ) ;
ssSetSolverNumZcSignals ( rtS , 13 ) ; ssSetModelZeroCrossings ( rtS ,
MdlZeroCrossings ) ; ssSetSolverConsecutiveZCsStepRelTol ( rtS ,
2.8421709430404007E-13 ) ; ssSetSolverMaxConsecutiveZCs ( rtS , 10000 ) ;
ssSetSolverConsecutiveZCsError ( rtS , 2 ) ; ssSetSolverMaskedZcDiagnostic (
rtS , 1 ) ; ssSetSolverIgnoredZcDiagnostic ( rtS , 1 ) ;
ssSetSolverMaxConsecutiveMinStep ( rtS , 1 ) ;
ssSetSolverShapePreserveControl ( rtS , 2 ) ; ssSetTNextTid ( rtS , INT_MIN )
; ssSetTNext ( rtS , rtMinusInf ) ; ssSetSolverNeedsReset ( rtS ) ;
ssSetNumNonsampledZCs ( rtS , 12 ) ; ssSetContStateDisabled ( rtS ,
contStatesDisabled ) ; ssSetSolverMaxConsecutiveMinStep ( rtS , 1 ) ; } {
ZCSigState * zc = ( ZCSigState * ) & rtPrevZCX ; ssSetPrevZCSigState ( rtS ,
zc ) ; } { rtPrevZCX . ah03fydvt5 = UNINITIALIZED_ZCSIG ; } ssSetChecksumVal
( rtS , 0 , 2453326653U ) ; ssSetChecksumVal ( rtS , 1 , 3278769589U ) ;
ssSetChecksumVal ( rtS , 2 , 1795305681U ) ; ssSetChecksumVal ( rtS , 3 ,
4235931226U ) ; { static const sysRanDType rtAlwaysEnabled =
SUBSYS_RAN_BC_ENABLE ; static RTWExtModeInfo rt_ExtModeInfo ; static const
sysRanDType * systemRan [ 18 ] ; gblRTWExtModeInfo = & rt_ExtModeInfo ;
ssSetRTWExtModeInfo ( rtS , & rt_ExtModeInfo ) ;
rteiSetSubSystemActiveVectorAddresses ( & rt_ExtModeInfo , systemRan ) ;
systemRan [ 0 ] = & rtAlwaysEnabled ; systemRan [ 1 ] = & rtAlwaysEnabled ;
systemRan [ 2 ] = & rtAlwaysEnabled ; systemRan [ 3 ] = & rtAlwaysEnabled ;
systemRan [ 4 ] = & rtAlwaysEnabled ; systemRan [ 5 ] = & rtAlwaysEnabled ;
systemRan [ 6 ] = & rtAlwaysEnabled ; systemRan [ 7 ] = & rtAlwaysEnabled ;
systemRan [ 8 ] = & rtAlwaysEnabled ; systemRan [ 9 ] = & rtAlwaysEnabled ;
systemRan [ 10 ] = & rtAlwaysEnabled ; systemRan [ 11 ] = & rtAlwaysEnabled ;
systemRan [ 12 ] = & rtAlwaysEnabled ; systemRan [ 13 ] = ( sysRanDType * ) &
rtDW . iln55qrgxk ; systemRan [ 14 ] = & rtAlwaysEnabled ; systemRan [ 15 ] =
& rtAlwaysEnabled ; systemRan [ 16 ] = & rtAlwaysEnabled ; systemRan [ 17 ] =
& rtAlwaysEnabled ; rteiSetModelMappingInfoPtr ( ssGetRTWExtModeInfo ( rtS )
, & ssGetModelMappingInfo ( rtS ) ) ; rteiSetChecksumsPtr (
ssGetRTWExtModeInfo ( rtS ) , ssGetChecksums ( rtS ) ) ; rteiSetTPtr (
ssGetRTWExtModeInfo ( rtS ) , ssGetTPtr ( rtS ) ) ; }
slsaDisallowedBlocksForSimTargetOP ( rtS ,
mr_Copy_of_trial_with_intermediate_tanks_GetSimStateDisallowedBlocks ) ;
slsaGetWorkFcnForSimTargetOP ( rtS ,
mr_Copy_of_trial_with_intermediate_tanks_GetDWork ) ;
slsaSetWorkFcnForSimTargetOP ( rtS ,
mr_Copy_of_trial_with_intermediate_tanks_SetDWork ) ;
rt_RapidReadMatFileAndUpdateParams ( rtS ) ; if ( ssGetErrorStatus ( rtS ) )
{ return rtS ; } return rtS ; }
#if defined(_MSC_VER)
#pragma optimize( "", on )
#endif
void MdlOutputsParameterSampleTime ( int_T tid ) { MdlOutputsTID3 ( tid ) ; }
