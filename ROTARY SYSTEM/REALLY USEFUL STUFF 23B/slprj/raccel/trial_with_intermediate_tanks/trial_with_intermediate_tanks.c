#include "trial_with_intermediate_tanks.h"
#include "rtwtypes.h"
#include "mwmathutil.h"
#include "trial_with_intermediate_tanks_private.h"
#include "rt_logging_mmi.h"
#include "trial_with_intermediate_tanks_capi.h"
#include "zero_crossing_types.h"
#include "trial_with_intermediate_tanks_dt.h"
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
#define oswwtj4xei (-1)
B rtB ; X rtX ; DW rtDW ; PrevZCX rtPrevZCX ; static SimStruct model_S ;
SimStruct * const rtS = & model_S ; void MdlInitialize ( void ) { boolean_T
tmp ; rtX . njnhyztmul = rtP . Integrator4_IC ; rtDW . di0xkhax5v = 1 ; if (
ssIsFirstInitCond ( rtS ) ) { rtX . cxtidlrvgx = 1.0E+7 ; tmp =
slIsRapidAcceleratorSimulating ( ) ; if ( tmp ) { tmp =
ssGetGlobalInitialStatesAvailable ( rtS ) ; rtDW . di0xkhax5v = ! tmp ; }
else { rtDW . di0xkhax5v = 1 ; } rtX . h1nzd5kpbz = 2.0E+7 ; } rtX .
anhq5q4adb = rtP . Integrator13_IC ; rtDW . dol3avfpna = 1 ; if (
ssIsFirstInitCond ( rtS ) ) { tmp = slIsRapidAcceleratorSimulating ( ) ; if (
tmp ) { tmp = ssGetGlobalInitialStatesAvailable ( rtS ) ; rtDW . dol3avfpna =
! tmp ; } else { rtDW . dol3avfpna = 1 ; } rtX . pmj0r0j4v2 = 1.0E+7 ; } rtDW
. mfw5n1ea1y = 1 ; if ( ssIsFirstInitCond ( rtS ) ) { tmp =
slIsRapidAcceleratorSimulating ( ) ; if ( tmp ) { tmp =
ssGetGlobalInitialStatesAvailable ( rtS ) ; rtDW . mfw5n1ea1y = ! tmp ; }
else { rtDW . mfw5n1ea1y = 1 ; } rtX . ith00xbmqs = 200.0 ; } rtX .
dusu5owudw = rtP . Integrator14_IC ; rtX . decqh12jsd = rtP . Integrator7_IC
; rtX . f3sxyupuw4 = rtP . Integrator12_IC ; rtX . javki5lnyj = rtP .
Integrator10_IC ; rtX . d1tt3gwruv = rtP . Integrator8_IC ; rtX . erk0rjtwtw
= rtP . Integrator9_IC ; rtX . m2od42nnln = rtP . Integrator11_IC ; rtDW .
oj5q0jsluk = 1 ; if ( ssIsFirstInitCond ( rtS ) ) { tmp =
slIsRapidAcceleratorSimulating ( ) ; if ( tmp ) { tmp =
ssGetGlobalInitialStatesAvailable ( rtS ) ; rtDW . oj5q0jsluk = ! tmp ; }
else { rtDW . oj5q0jsluk = 1 ; } } rtX . djreaptdva = rtP . Integrator6_IC ;
rtX . pnlrwk52wo = rtP . Integrator16_IC ; rtX . mh2ntyj1vq = rtP .
Integrator5_IC ; rtDW . mog1ypd0na = rtP . Delay_InitialCondition ; rtX .
kh4uwqwwas = rtP . Integrator15_IC ; rtX . lxsclehbdn [ 0 ] = 0.0 ; rtX .
herd5wji3e [ 0 ] = 0.0 ; rtX . ekhpkuwiae [ 0 ] = 0.0 ; rtX . k2eyyio13c [ 0
] = 0.0 ; rtX . l1m0lcksjp [ 0 ] = 0.0 ; rtX . lxsclehbdn [ 1 ] = 0.0 ; rtX .
herd5wji3e [ 1 ] = 0.0 ; rtX . ekhpkuwiae [ 1 ] = 0.0 ; rtX . k2eyyio13c [ 1
] = 0.0 ; rtX . l1m0lcksjp [ 1 ] = 0.0 ; rtDW . f40yghkee4 = false ; rtDW .
brdpyuefrq = oswwtj4xei ; rtDW . cbdxrybiau = false ; rtDW . bx2j5ck3fy =
oswwtj4xei ; rtDW . nhmw0vweas = false ; rtDW . nfkviszx2j = oswwtj4xei ;
rtDW . mx1mzmmudm = false ; rtDW . ojcciswcd3 = oswwtj4xei ; rtDW .
cifi3jtean = false ; rtDW . brka54cpul = oswwtj4xei ; rtDW . g3swjv3kh4 =
false ; rtDW . d1ijg3op2k = oswwtj4xei ; rtDW . jvf3pz3ask = false ; rtDW .
kjn3oyfrac = oswwtj4xei ; rtB . a1drd1cazf = rtP . Out1_Y0 ; rtDW .
cy4ormxj51 = false ; rtDW . d3q0cyta13 = oswwtj4xei ; rtDW . ob1pklob1o =
false ; rtDW . ncj5ml0ihc = oswwtj4xei ; rtDW . pke0ukp1xy = false ; rtDW .
esun0jx3xw = oswwtj4xei ; rtDW . dia2paiwqh = false ; rtDW . fzispyfeqr =
oswwtj4xei ; rtDW . jelszbraal = false ; rtDW . puivctwyq5 = oswwtj4xei ;
rtDW . dwzumrr4bn = false ; rtDW . egiikq3opb = oswwtj4xei ; rtDW .
cwu3pvv3pv = false ; rtDW . mevs4smziy = oswwtj4xei ; rtDW . au1n1jxucn =
false ; rtDW . c2jtsctdjh = oswwtj4xei ; rtDW . enb1yfnvht = false ; rtDW .
ndyrr3tae0 = oswwtj4xei ; } void MdlStart ( void ) { { bool
externalInputIsInDatasetFormat = false ; void * pISigstreamManager =
rt_GetISigstreamManager ( rtS ) ;
rtwISigstreamManagerGetInputIsInDatasetFormat ( pISigstreamManager , &
externalInputIsInDatasetFormat ) ; if ( externalInputIsInDatasetFormat ) { }
} { { { bool isStreamoutAlreadyRegistered = false ; { sdiSignalSourceInfoU
srcInfo ; sdiLabelU loggedName = sdiGetLabelFromChars ( "P3" ) ; sdiLabelU
origSigName = sdiGetLabelFromChars ( "P3" ) ; sdiLabelU propName =
sdiGetLabelFromChars ( "P3" ) ; sdiLabelU blockPath = sdiGetLabelFromChars (
"trial_with_intermediate_tanks/To Workspace" ) ; sdiLabelU blockSID =
sdiGetLabelFromChars ( "" ) ; sdiLabelU subPath = sdiGetLabelFromChars ( "" )
; sdiDims sigDims ; sdiLabelU sigName = sdiGetLabelFromChars ( "P3" ) ;
sdiAsyncRepoDataTypeHandle hDT = sdiAsyncRepoGetBuiltInDataTypeHandle (
DATA_TYPE_DOUBLE ) ; { sdiComplexity sigComplexity = REAL ;
sdiSampleTimeContinuity stCont = SAMPLE_TIME_CONTINUOUS ; int_T sigDimsArray
[ 1 ] = { 1 } ; sigDims . nDims = 1 ; sigDims . dimensions = sigDimsArray ;
srcInfo . numBlockPathElems = 1 ; srcInfo . fullBlockPath = ( sdiFullBlkPathU
) & blockPath ; srcInfo . SID = ( sdiSignalIDU ) & blockSID ; srcInfo .
subPath = subPath ; srcInfo . portIndex = 0 + 1 ; srcInfo . signalName =
sigName ; srcInfo . sigSourceUUID = 0 ; rtDW . cpa2ih4ymh . AQHandles =
sdiStartAsyncioQueueCreation ( hDT , & srcInfo , rt_dataMapInfo . mmi .
InstanceMap . fullPath , "d6ef521d-7619-47c3-ad8b-0a874888afc6" ,
sigComplexity , & sigDims , DIMENSIONS_MODE_FIXED , stCont , "" ) ;
sdiCompleteAsyncioQueueCreation ( rtDW . cpa2ih4ymh . AQHandles , hDT , &
srcInfo ) ; if ( rtDW . cpa2ih4ymh . AQHandles ) {
sdiSetSignalSampleTimeString ( rtDW . cpa2ih4ymh . AQHandles , "Continuous" ,
0.0 , ssGetTFinal ( rtS ) ) ; sdiSetSignalRefRate ( rtDW . cpa2ih4ymh .
AQHandles , 0.0 ) ; sdiSetRunStartTime ( rtDW . cpa2ih4ymh . AQHandles ,
ssGetTaskTime ( rtS , 0 ) ) ; sdiAsyncRepoSetSignalExportSettings ( rtDW .
cpa2ih4ymh . AQHandles , 1 , 0 ) ; sdiAsyncRepoSetSignalExportName ( rtDW .
cpa2ih4ymh . AQHandles , loggedName , origSigName , propName ) ;
sdiAsyncRepoSetBlockPathDomain ( rtDW . cpa2ih4ymh . AQHandles ) ; }
sdiFreeLabel ( sigName ) ; sdiFreeLabel ( loggedName ) ; sdiFreeLabel (
origSigName ) ; sdiFreeLabel ( propName ) ; sdiFreeLabel ( blockPath ) ;
sdiFreeLabel ( blockSID ) ; sdiFreeLabel ( subPath ) ; } } if ( !
isStreamoutAlreadyRegistered ) { { sdiLabelU varName = sdiGetLabelFromChars (
"P1" ) ; sdiRegisterWksVariable ( rtDW . cpa2ih4ymh . AQHandles , varName ,
"timeseries" ) ; sdiFreeLabel ( varName ) ; } } } } } { { { bool
isStreamoutAlreadyRegistered = false ; { sdiSignalSourceInfoU srcInfo ;
sdiLabelU loggedName = sdiGetLabelFromChars ( "Integrator7" ) ; sdiLabelU
origSigName = sdiGetLabelFromChars ( "" ) ; sdiLabelU propName =
sdiGetLabelFromChars ( "Integrator7" ) ; sdiLabelU blockPath =
sdiGetLabelFromChars ( "trial_with_intermediate_tanks/To Workspace1" ) ;
sdiLabelU blockSID = sdiGetLabelFromChars ( "" ) ; sdiLabelU subPath =
sdiGetLabelFromChars ( "" ) ; sdiDims sigDims ; sdiLabelU sigName =
sdiGetLabelFromChars ( "Integrator7" ) ; sdiAsyncRepoDataTypeHandle hDT =
sdiAsyncRepoGetBuiltInDataTypeHandle ( DATA_TYPE_DOUBLE ) ; { sdiComplexity
sigComplexity = REAL ; sdiSampleTimeContinuity stCont =
SAMPLE_TIME_CONTINUOUS ; int_T sigDimsArray [ 1 ] = { 1 } ; sigDims . nDims =
1 ; sigDims . dimensions = sigDimsArray ; srcInfo . numBlockPathElems = 1 ;
srcInfo . fullBlockPath = ( sdiFullBlkPathU ) & blockPath ; srcInfo . SID = (
sdiSignalIDU ) & blockSID ; srcInfo . subPath = subPath ; srcInfo . portIndex
= 0 + 1 ; srcInfo . signalName = sigName ; srcInfo . sigSourceUUID = 0 ; rtDW
. igylchmu0t . AQHandles = sdiStartAsyncioQueueCreation ( hDT , & srcInfo ,
rt_dataMapInfo . mmi . InstanceMap . fullPath ,
"f7e8adef-296f-466c-a210-e16474237178" , sigComplexity , & sigDims ,
DIMENSIONS_MODE_FIXED , stCont , "" ) ; sdiCompleteAsyncioQueueCreation (
rtDW . igylchmu0t . AQHandles , hDT , & srcInfo ) ; if ( rtDW . igylchmu0t .
AQHandles ) { sdiSetSignalSampleTimeString ( rtDW . igylchmu0t . AQHandles ,
"Continuous" , 0.0 , ssGetTFinal ( rtS ) ) ; sdiSetSignalRefRate ( rtDW .
igylchmu0t . AQHandles , 0.0 ) ; sdiSetRunStartTime ( rtDW . igylchmu0t .
AQHandles , ssGetTaskTime ( rtS , 0 ) ) ; sdiAsyncRepoSetSignalExportSettings
( rtDW . igylchmu0t . AQHandles , 1 , 0 ) ; sdiAsyncRepoSetSignalExportName (
rtDW . igylchmu0t . AQHandles , loggedName , origSigName , propName ) ;
sdiAsyncRepoSetBlockPathDomain ( rtDW . igylchmu0t . AQHandles ) ; }
sdiFreeLabel ( sigName ) ; sdiFreeLabel ( loggedName ) ; sdiFreeLabel (
origSigName ) ; sdiFreeLabel ( propName ) ; sdiFreeLabel ( blockPath ) ;
sdiFreeLabel ( blockSID ) ; sdiFreeLabel ( subPath ) ; } } if ( !
isStreamoutAlreadyRegistered ) { { sdiLabelU varName = sdiGetLabelFromChars (
"Regen" ) ; sdiRegisterWksVariable ( rtDW . igylchmu0t . AQHandles , varName
, "timeseries" ) ; sdiFreeLabel ( varName ) ; } } } } } { { { bool
isStreamoutAlreadyRegistered = false ; { sdiSignalSourceInfoU srcInfo ;
sdiLabelU loggedName = sdiGetLabelFromChars ( "Electric Torque Controller" )
; sdiLabelU origSigName = sdiGetLabelFromChars ( "" ) ; sdiLabelU propName =
sdiGetLabelFromChars ( "Electric Torque Controller" ) ; sdiLabelU blockPath =
sdiGetLabelFromChars ( "trial_with_intermediate_tanks/To Workspace10" ) ;
sdiLabelU blockSID = sdiGetLabelFromChars ( "" ) ; sdiLabelU subPath =
sdiGetLabelFromChars ( "" ) ; sdiDims sigDims ; sdiLabelU sigName =
sdiGetLabelFromChars ( "Electric Torque Controller" ) ;
sdiAsyncRepoDataTypeHandle hDT = sdiAsyncRepoGetBuiltInDataTypeHandle (
DATA_TYPE_DOUBLE ) ; { sdiComplexity sigComplexity = REAL ;
sdiSampleTimeContinuity stCont = SAMPLE_TIME_CONTINUOUS ; int_T sigDimsArray
[ 1 ] = { 1 } ; sigDims . nDims = 1 ; sigDims . dimensions = sigDimsArray ;
srcInfo . numBlockPathElems = 1 ; srcInfo . fullBlockPath = ( sdiFullBlkPathU
) & blockPath ; srcInfo . SID = ( sdiSignalIDU ) & blockSID ; srcInfo .
subPath = subPath ; srcInfo . portIndex = 0 + 1 ; srcInfo . signalName =
sigName ; srcInfo . sigSourceUUID = 0 ; rtDW . aqjm4jt1o5 . AQHandles =
sdiStartAsyncioQueueCreation ( hDT , & srcInfo , rt_dataMapInfo . mmi .
InstanceMap . fullPath , "525b40c8-4150-4512-8b99-66ac61306d67" ,
sigComplexity , & sigDims , DIMENSIONS_MODE_FIXED , stCont , "" ) ;
sdiCompleteAsyncioQueueCreation ( rtDW . aqjm4jt1o5 . AQHandles , hDT , &
srcInfo ) ; if ( rtDW . aqjm4jt1o5 . AQHandles ) {
sdiSetSignalSampleTimeString ( rtDW . aqjm4jt1o5 . AQHandles , "Continuous" ,
0.0 , ssGetTFinal ( rtS ) ) ; sdiSetSignalRefRate ( rtDW . aqjm4jt1o5 .
AQHandles , 0.0 ) ; sdiSetRunStartTime ( rtDW . aqjm4jt1o5 . AQHandles ,
ssGetTaskTime ( rtS , 0 ) ) ; sdiAsyncRepoSetSignalExportSettings ( rtDW .
aqjm4jt1o5 . AQHandles , 1 , 0 ) ; sdiAsyncRepoSetSignalExportName ( rtDW .
aqjm4jt1o5 . AQHandles , loggedName , origSigName , propName ) ;
sdiAsyncRepoSetBlockPathDomain ( rtDW . aqjm4jt1o5 . AQHandles ) ; }
sdiFreeLabel ( sigName ) ; sdiFreeLabel ( loggedName ) ; sdiFreeLabel (
origSigName ) ; sdiFreeLabel ( propName ) ; sdiFreeLabel ( blockPath ) ;
sdiFreeLabel ( blockSID ) ; sdiFreeLabel ( subPath ) ; } } if ( !
isStreamoutAlreadyRegistered ) { { sdiLabelU varName = sdiGetLabelFromChars (
"T_elec" ) ; sdiRegisterWksVariable ( rtDW . aqjm4jt1o5 . AQHandles , varName
, "timeseries" ) ; sdiFreeLabel ( varName ) ; } } } } } { { { bool
isStreamoutAlreadyRegistered = false ; { sdiSignalSourceInfoU srcInfo ;
sdiLabelU loggedName = sdiGetLabelFromChars ( "Integrator6" ) ; sdiLabelU
origSigName = sdiGetLabelFromChars ( "" ) ; sdiLabelU propName =
sdiGetLabelFromChars ( "Integrator6" ) ; sdiLabelU blockPath =
sdiGetLabelFromChars ( "trial_with_intermediate_tanks/To Workspace11" ) ;
sdiLabelU blockSID = sdiGetLabelFromChars ( "" ) ; sdiLabelU subPath =
sdiGetLabelFromChars ( "" ) ; sdiDims sigDims ; sdiLabelU sigName =
sdiGetLabelFromChars ( "Integrator6" ) ; sdiAsyncRepoDataTypeHandle hDT =
sdiAsyncRepoGetBuiltInDataTypeHandle ( DATA_TYPE_DOUBLE ) ; { sdiComplexity
sigComplexity = REAL ; sdiSampleTimeContinuity stCont =
SAMPLE_TIME_CONTINUOUS ; int_T sigDimsArray [ 1 ] = { 1 } ; sigDims . nDims =
1 ; sigDims . dimensions = sigDimsArray ; srcInfo . numBlockPathElems = 1 ;
srcInfo . fullBlockPath = ( sdiFullBlkPathU ) & blockPath ; srcInfo . SID = (
sdiSignalIDU ) & blockSID ; srcInfo . subPath = subPath ; srcInfo . portIndex
= 0 + 1 ; srcInfo . signalName = sigName ; srcInfo . sigSourceUUID = 0 ; rtDW
. gzhbrbv4c4 . AQHandles = sdiStartAsyncioQueueCreation ( hDT , & srcInfo ,
rt_dataMapInfo . mmi . InstanceMap . fullPath ,
"2fc4b6e7-bfef-4707-bf79-e0dd9f32e63f" , sigComplexity , & sigDims ,
DIMENSIONS_MODE_FIXED , stCont , "" ) ; sdiCompleteAsyncioQueueCreation (
rtDW . gzhbrbv4c4 . AQHandles , hDT , & srcInfo ) ; if ( rtDW . gzhbrbv4c4 .
AQHandles ) { sdiSetSignalSampleTimeString ( rtDW . gzhbrbv4c4 . AQHandles ,
"Continuous" , 0.0 , ssGetTFinal ( rtS ) ) ; sdiSetSignalRefRate ( rtDW .
gzhbrbv4c4 . AQHandles , 0.0 ) ; sdiSetRunStartTime ( rtDW . gzhbrbv4c4 .
AQHandles , ssGetTaskTime ( rtS , 0 ) ) ; sdiAsyncRepoSetSignalExportSettings
( rtDW . gzhbrbv4c4 . AQHandles , 1 , 0 ) ; sdiAsyncRepoSetSignalExportName (
rtDW . gzhbrbv4c4 . AQHandles , loggedName , origSigName , propName ) ;
sdiAsyncRepoSetBlockPathDomain ( rtDW . gzhbrbv4c4 . AQHandles ) ; }
sdiFreeLabel ( sigName ) ; sdiFreeLabel ( loggedName ) ; sdiFreeLabel (
origSigName ) ; sdiFreeLabel ( propName ) ; sdiFreeLabel ( blockPath ) ;
sdiFreeLabel ( blockSID ) ; sdiFreeLabel ( subPath ) ; } } if ( !
isStreamoutAlreadyRegistered ) { { sdiLabelU varName = sdiGetLabelFromChars (
"Energy_In" ) ; sdiRegisterWksVariable ( rtDW . gzhbrbv4c4 . AQHandles ,
varName , "timeseries" ) ; sdiFreeLabel ( varName ) ; } } } } } { { { bool
isStreamoutAlreadyRegistered = false ; { sdiSignalSourceInfoU srcInfo ;
sdiLabelU loggedName = sdiGetLabelFromChars ( "Gain3" ) ; sdiLabelU
origSigName = sdiGetLabelFromChars ( "" ) ; sdiLabelU propName =
sdiGetLabelFromChars ( "Gain3" ) ; sdiLabelU blockPath = sdiGetLabelFromChars
( "trial_with_intermediate_tanks/To Workspace2" ) ; sdiLabelU blockSID =
sdiGetLabelFromChars ( "" ) ; sdiLabelU subPath = sdiGetLabelFromChars ( "" )
; sdiDims sigDims ; sdiLabelU sigName = sdiGetLabelFromChars ( "Gain3" ) ;
sdiAsyncRepoDataTypeHandle hDT = sdiAsyncRepoGetBuiltInDataTypeHandle (
DATA_TYPE_DOUBLE ) ; { sdiComplexity sigComplexity = REAL ;
sdiSampleTimeContinuity stCont = SAMPLE_TIME_CONTINUOUS ; int_T sigDimsArray
[ 1 ] = { 1 } ; sigDims . nDims = 1 ; sigDims . dimensions = sigDimsArray ;
srcInfo . numBlockPathElems = 1 ; srcInfo . fullBlockPath = ( sdiFullBlkPathU
) & blockPath ; srcInfo . SID = ( sdiSignalIDU ) & blockSID ; srcInfo .
subPath = subPath ; srcInfo . portIndex = 0 + 1 ; srcInfo . signalName =
sigName ; srcInfo . sigSourceUUID = 0 ; rtDW . htile1lzbu . AQHandles =
sdiStartAsyncioQueueCreation ( hDT , & srcInfo , rt_dataMapInfo . mmi .
InstanceMap . fullPath , "d473f28d-961c-4d52-886e-1e91ef496d5b" ,
sigComplexity , & sigDims , DIMENSIONS_MODE_FIXED , stCont , "" ) ;
sdiCompleteAsyncioQueueCreation ( rtDW . htile1lzbu . AQHandles , hDT , &
srcInfo ) ; if ( rtDW . htile1lzbu . AQHandles ) {
sdiSetSignalSampleTimeString ( rtDW . htile1lzbu . AQHandles , "Continuous" ,
0.0 , ssGetTFinal ( rtS ) ) ; sdiSetSignalRefRate ( rtDW . htile1lzbu .
AQHandles , 0.0 ) ; sdiSetRunStartTime ( rtDW . htile1lzbu . AQHandles ,
ssGetTaskTime ( rtS , 0 ) ) ; sdiAsyncRepoSetSignalExportSettings ( rtDW .
htile1lzbu . AQHandles , 1 , 0 ) ; sdiAsyncRepoSetSignalExportName ( rtDW .
htile1lzbu . AQHandles , loggedName , origSigName , propName ) ;
sdiAsyncRepoSetBlockPathDomain ( rtDW . htile1lzbu . AQHandles ) ; }
sdiFreeLabel ( sigName ) ; sdiFreeLabel ( loggedName ) ; sdiFreeLabel (
origSigName ) ; sdiFreeLabel ( propName ) ; sdiFreeLabel ( blockPath ) ;
sdiFreeLabel ( blockSID ) ; sdiFreeLabel ( subPath ) ; } } if ( !
isStreamoutAlreadyRegistered ) { { sdiLabelU varName = sdiGetLabelFromChars (
"KE" ) ; sdiRegisterWksVariable ( rtDW . htile1lzbu . AQHandles , varName ,
"timeseries" ) ; sdiFreeLabel ( varName ) ; } } } } } { { { bool
isStreamoutAlreadyRegistered = false ; { sdiSignalSourceInfoU srcInfo ;
sdiLabelU loggedName = sdiGetLabelFromChars ( "Event Time" ) ; sdiLabelU
origSigName = sdiGetLabelFromChars ( "Event Time" ) ; sdiLabelU propName =
sdiGetLabelFromChars ( "Event Time" ) ; sdiLabelU blockPath =
sdiGetLabelFromChars ( "trial_with_intermediate_tanks/To Workspace3" ) ;
sdiLabelU blockSID = sdiGetLabelFromChars ( "" ) ; sdiLabelU subPath =
sdiGetLabelFromChars ( "" ) ; sdiDims sigDims ; sdiLabelU sigName =
sdiGetLabelFromChars ( "Event Time" ) ; sdiAsyncRepoDataTypeHandle hDT =
sdiAsyncRepoGetBuiltInDataTypeHandle ( DATA_TYPE_DOUBLE ) ; { sdiComplexity
sigComplexity = REAL ; sdiSampleTimeContinuity stCont = SAMPLE_TIME_DISCRETE
; int_T sigDimsArray [ 1 ] = { 1 } ; sigDims . nDims = 1 ; sigDims .
dimensions = sigDimsArray ; srcInfo . numBlockPathElems = 1 ; srcInfo .
fullBlockPath = ( sdiFullBlkPathU ) & blockPath ; srcInfo . SID = (
sdiSignalIDU ) & blockSID ; srcInfo . subPath = subPath ; srcInfo . portIndex
= 0 + 1 ; srcInfo . signalName = sigName ; srcInfo . sigSourceUUID = 0 ; rtDW
. p1ayluai4b . AQHandles = sdiStartAsyncioQueueCreation ( hDT , & srcInfo ,
rt_dataMapInfo . mmi . InstanceMap . fullPath ,
"4a933b4d-db2a-4d20-bc88-f86bd946adbb" , sigComplexity , & sigDims ,
DIMENSIONS_MODE_FIXED , stCont , "" ) ; sdiCompleteAsyncioQueueCreation (
rtDW . p1ayluai4b . AQHandles , hDT , & srcInfo ) ; if ( rtDW . p1ayluai4b .
AQHandles ) { sdiSetSignalSampleTimeString ( rtDW . p1ayluai4b . AQHandles ,
"Continuous" , 0.0 , ssGetTFinal ( rtS ) ) ; sdiSetSignalRefRate ( rtDW .
p1ayluai4b . AQHandles , 0.0 ) ; sdiSetRunStartTime ( rtDW . p1ayluai4b .
AQHandles , ssGetTaskTime ( rtS , 1 ) ) ; sdiAsyncRepoSetSignalExportSettings
( rtDW . p1ayluai4b . AQHandles , 1 , 0 ) ; sdiAsyncRepoSetSignalExportName (
rtDW . p1ayluai4b . AQHandles , loggedName , origSigName , propName ) ;
sdiAsyncRepoSetBlockPathDomain ( rtDW . p1ayluai4b . AQHandles ) ; }
sdiFreeLabel ( sigName ) ; sdiFreeLabel ( loggedName ) ; sdiFreeLabel (
origSigName ) ; sdiFreeLabel ( propName ) ; sdiFreeLabel ( blockPath ) ;
sdiFreeLabel ( blockSID ) ; sdiFreeLabel ( subPath ) ; } } if ( !
isStreamoutAlreadyRegistered ) { { sdiLabelU varName = sdiGetLabelFromChars (
"Time_settle" ) ; sdiRegisterWksVariable ( rtDW . p1ayluai4b . AQHandles ,
varName , "timeseries" ) ; sdiFreeLabel ( varName ) ; } } } } } { { { bool
isStreamoutAlreadyRegistered = false ; { sdiSignalSourceInfoU srcInfo ;
sdiLabelU loggedName = sdiGetLabelFromChars ( "omega" ) ; sdiLabelU
origSigName = sdiGetLabelFromChars ( "omega" ) ; sdiLabelU propName =
sdiGetLabelFromChars ( "omega" ) ; sdiLabelU blockPath = sdiGetLabelFromChars
( "trial_with_intermediate_tanks/To Workspace4" ) ; sdiLabelU blockSID =
sdiGetLabelFromChars ( "" ) ; sdiLabelU subPath = sdiGetLabelFromChars ( "" )
; sdiDims sigDims ; sdiLabelU sigName = sdiGetLabelFromChars ( "omega" ) ;
sdiAsyncRepoDataTypeHandle hDT = sdiAsyncRepoGetBuiltInDataTypeHandle (
DATA_TYPE_DOUBLE ) ; { sdiComplexity sigComplexity = REAL ;
sdiSampleTimeContinuity stCont = SAMPLE_TIME_CONTINUOUS ; int_T sigDimsArray
[ 1 ] = { 1 } ; sigDims . nDims = 1 ; sigDims . dimensions = sigDimsArray ;
srcInfo . numBlockPathElems = 1 ; srcInfo . fullBlockPath = ( sdiFullBlkPathU
) & blockPath ; srcInfo . SID = ( sdiSignalIDU ) & blockSID ; srcInfo .
subPath = subPath ; srcInfo . portIndex = 0 + 1 ; srcInfo . signalName =
sigName ; srcInfo . sigSourceUUID = 0 ; rtDW . gjtuyi3wgw . AQHandles =
sdiStartAsyncioQueueCreation ( hDT , & srcInfo , rt_dataMapInfo . mmi .
InstanceMap . fullPath , "3f5bbb3d-2321-459d-bdf3-8b84e1b4dc5b" ,
sigComplexity , & sigDims , DIMENSIONS_MODE_FIXED , stCont , "" ) ;
sdiCompleteAsyncioQueueCreation ( rtDW . gjtuyi3wgw . AQHandles , hDT , &
srcInfo ) ; if ( rtDW . gjtuyi3wgw . AQHandles ) {
sdiSetSignalSampleTimeString ( rtDW . gjtuyi3wgw . AQHandles , "Continuous" ,
0.0 , ssGetTFinal ( rtS ) ) ; sdiSetSignalRefRate ( rtDW . gjtuyi3wgw .
AQHandles , 0.0 ) ; sdiSetRunStartTime ( rtDW . gjtuyi3wgw . AQHandles ,
ssGetTaskTime ( rtS , 0 ) ) ; sdiAsyncRepoSetSignalExportSettings ( rtDW .
gjtuyi3wgw . AQHandles , 1 , 0 ) ; sdiAsyncRepoSetSignalExportName ( rtDW .
gjtuyi3wgw . AQHandles , loggedName , origSigName , propName ) ;
sdiAsyncRepoSetBlockPathDomain ( rtDW . gjtuyi3wgw . AQHandles ) ; }
sdiFreeLabel ( sigName ) ; sdiFreeLabel ( loggedName ) ; sdiFreeLabel (
origSigName ) ; sdiFreeLabel ( propName ) ; sdiFreeLabel ( blockPath ) ;
sdiFreeLabel ( blockSID ) ; sdiFreeLabel ( subPath ) ; } } if ( !
isStreamoutAlreadyRegistered ) { { sdiLabelU varName = sdiGetLabelFromChars (
"omega" ) ; sdiRegisterWksVariable ( rtDW . gjtuyi3wgw . AQHandles , varName
, "timeseries" ) ; sdiFreeLabel ( varName ) ; } } } } } { { { bool
isStreamoutAlreadyRegistered = false ; { sdiSignalSourceInfoU srcInfo ;
sdiLabelU loggedName = sdiGetLabelFromChars ( "Gain1" ) ; sdiLabelU
origSigName = sdiGetLabelFromChars ( "" ) ; sdiLabelU propName =
sdiGetLabelFromChars ( "Gain1" ) ; sdiLabelU blockPath = sdiGetLabelFromChars
( "trial_with_intermediate_tanks/To Workspace5" ) ; sdiLabelU blockSID =
sdiGetLabelFromChars ( "" ) ; sdiLabelU subPath = sdiGetLabelFromChars ( "" )
; sdiDims sigDims ; sdiLabelU sigName = sdiGetLabelFromChars ( "Gain1" ) ;
sdiAsyncRepoDataTypeHandle hDT = sdiAsyncRepoGetBuiltInDataTypeHandle (
DATA_TYPE_DOUBLE ) ; { sdiComplexity sigComplexity = REAL ;
sdiSampleTimeContinuity stCont = SAMPLE_TIME_CONTINUOUS ; int_T sigDimsArray
[ 1 ] = { 1 } ; sigDims . nDims = 1 ; sigDims . dimensions = sigDimsArray ;
srcInfo . numBlockPathElems = 1 ; srcInfo . fullBlockPath = ( sdiFullBlkPathU
) & blockPath ; srcInfo . SID = ( sdiSignalIDU ) & blockSID ; srcInfo .
subPath = subPath ; srcInfo . portIndex = 0 + 1 ; srcInfo . signalName =
sigName ; srcInfo . sigSourceUUID = 0 ; rtDW . jd4fubbir0 . AQHandles =
sdiStartAsyncioQueueCreation ( hDT , & srcInfo , rt_dataMapInfo . mmi .
InstanceMap . fullPath , "0285a53f-fbfb-43ee-a2f3-421b4ab86fc4" ,
sigComplexity , & sigDims , DIMENSIONS_MODE_FIXED , stCont , "" ) ;
sdiCompleteAsyncioQueueCreation ( rtDW . jd4fubbir0 . AQHandles , hDT , &
srcInfo ) ; if ( rtDW . jd4fubbir0 . AQHandles ) {
sdiSetSignalSampleTimeString ( rtDW . jd4fubbir0 . AQHandles , "Continuous" ,
0.0 , ssGetTFinal ( rtS ) ) ; sdiSetSignalRefRate ( rtDW . jd4fubbir0 .
AQHandles , 0.0 ) ; sdiSetRunStartTime ( rtDW . jd4fubbir0 . AQHandles ,
ssGetTaskTime ( rtS , 0 ) ) ; sdiAsyncRepoSetSignalExportSettings ( rtDW .
jd4fubbir0 . AQHandles , 1 , 0 ) ; sdiAsyncRepoSetSignalExportName ( rtDW .
jd4fubbir0 . AQHandles , loggedName , origSigName , propName ) ;
sdiAsyncRepoSetBlockPathDomain ( rtDW . jd4fubbir0 . AQHandles ) ; }
sdiFreeLabel ( sigName ) ; sdiFreeLabel ( loggedName ) ; sdiFreeLabel (
origSigName ) ; sdiFreeLabel ( propName ) ; sdiFreeLabel ( blockPath ) ;
sdiFreeLabel ( blockSID ) ; sdiFreeLabel ( subPath ) ; } } if ( !
isStreamoutAlreadyRegistered ) { { sdiLabelU varName = sdiGetLabelFromChars (
"mech_power" ) ; sdiRegisterWksVariable ( rtDW . jd4fubbir0 . AQHandles ,
varName , "timeseries" ) ; sdiFreeLabel ( varName ) ; } } } } } { { { bool
isStreamoutAlreadyRegistered = false ; { sdiSignalSourceInfoU srcInfo ;
sdiLabelU loggedName = sdiGetLabelFromChars ( "Product2" ) ; sdiLabelU
origSigName = sdiGetLabelFromChars ( "" ) ; sdiLabelU propName =
sdiGetLabelFromChars ( "Product2" ) ; sdiLabelU blockPath =
sdiGetLabelFromChars ( "trial_with_intermediate_tanks/To Workspace6" ) ;
sdiLabelU blockSID = sdiGetLabelFromChars ( "" ) ; sdiLabelU subPath =
sdiGetLabelFromChars ( "" ) ; sdiDims sigDims ; sdiLabelU sigName =
sdiGetLabelFromChars ( "Product2" ) ; sdiAsyncRepoDataTypeHandle hDT =
sdiAsyncRepoGetBuiltInDataTypeHandle ( DATA_TYPE_DOUBLE ) ; { sdiComplexity
sigComplexity = REAL ; sdiSampleTimeContinuity stCont =
SAMPLE_TIME_CONTINUOUS ; int_T sigDimsArray [ 1 ] = { 1 } ; sigDims . nDims =
1 ; sigDims . dimensions = sigDimsArray ; srcInfo . numBlockPathElems = 1 ;
srcInfo . fullBlockPath = ( sdiFullBlkPathU ) & blockPath ; srcInfo . SID = (
sdiSignalIDU ) & blockSID ; srcInfo . subPath = subPath ; srcInfo . portIndex
= 0 + 1 ; srcInfo . signalName = sigName ; srcInfo . sigSourceUUID = 0 ; rtDW
. kbuynfluan . AQHandles = sdiStartAsyncioQueueCreation ( hDT , & srcInfo ,
rt_dataMapInfo . mmi . InstanceMap . fullPath ,
"d2bddb2c-813e-4052-9a7a-d879ca4fe9ba" , sigComplexity , & sigDims ,
DIMENSIONS_MODE_FIXED , stCont , "" ) ; sdiCompleteAsyncioQueueCreation (
rtDW . kbuynfluan . AQHandles , hDT , & srcInfo ) ; if ( rtDW . kbuynfluan .
AQHandles ) { sdiSetSignalSampleTimeString ( rtDW . kbuynfluan . AQHandles ,
"Continuous" , 0.0 , ssGetTFinal ( rtS ) ) ; sdiSetSignalRefRate ( rtDW .
kbuynfluan . AQHandles , 0.0 ) ; sdiSetRunStartTime ( rtDW . kbuynfluan .
AQHandles , ssGetTaskTime ( rtS , 0 ) ) ; sdiAsyncRepoSetSignalExportSettings
( rtDW . kbuynfluan . AQHandles , 1 , 0 ) ; sdiAsyncRepoSetSignalExportName (
rtDW . kbuynfluan . AQHandles , loggedName , origSigName , propName ) ;
sdiAsyncRepoSetBlockPathDomain ( rtDW . kbuynfluan . AQHandles ) ; }
sdiFreeLabel ( sigName ) ; sdiFreeLabel ( loggedName ) ; sdiFreeLabel (
origSigName ) ; sdiFreeLabel ( propName ) ; sdiFreeLabel ( blockPath ) ;
sdiFreeLabel ( blockSID ) ; sdiFreeLabel ( subPath ) ; } } if ( !
isStreamoutAlreadyRegistered ) { { sdiLabelU varName = sdiGetLabelFromChars (
"Regen_power" ) ; sdiRegisterWksVariable ( rtDW . kbuynfluan . AQHandles ,
varName , "timeseries" ) ; sdiFreeLabel ( varName ) ; } } } } } { { { bool
isStreamoutAlreadyRegistered = false ; { sdiSignalSourceInfoU srcInfo ;
sdiLabelU loggedName = sdiGetLabelFromChars ( "Add2" ) ; sdiLabelU
origSigName = sdiGetLabelFromChars ( "" ) ; sdiLabelU propName =
sdiGetLabelFromChars ( "Add2" ) ; sdiLabelU blockPath = sdiGetLabelFromChars
( "trial_with_intermediate_tanks/To Workspace7" ) ; sdiLabelU blockSID =
sdiGetLabelFromChars ( "" ) ; sdiLabelU subPath = sdiGetLabelFromChars ( "" )
; sdiDims sigDims ; sdiLabelU sigName = sdiGetLabelFromChars ( "Add2" ) ;
sdiAsyncRepoDataTypeHandle hDT = sdiAsyncRepoGetBuiltInDataTypeHandle (
DATA_TYPE_DOUBLE ) ; { sdiComplexity sigComplexity = REAL ;
sdiSampleTimeContinuity stCont = SAMPLE_TIME_CONTINUOUS ; int_T sigDimsArray
[ 1 ] = { 1 } ; sigDims . nDims = 1 ; sigDims . dimensions = sigDimsArray ;
srcInfo . numBlockPathElems = 1 ; srcInfo . fullBlockPath = ( sdiFullBlkPathU
) & blockPath ; srcInfo . SID = ( sdiSignalIDU ) & blockSID ; srcInfo .
subPath = subPath ; srcInfo . portIndex = 0 + 1 ; srcInfo . signalName =
sigName ; srcInfo . sigSourceUUID = 0 ; rtDW . j5pyh2o05c . AQHandles =
sdiStartAsyncioQueueCreation ( hDT , & srcInfo , rt_dataMapInfo . mmi .
InstanceMap . fullPath , "57216774-7067-4811-bd25-0bbd2b0c8c64" ,
sigComplexity , & sigDims , DIMENSIONS_MODE_FIXED , stCont , "" ) ;
sdiCompleteAsyncioQueueCreation ( rtDW . j5pyh2o05c . AQHandles , hDT , &
srcInfo ) ; if ( rtDW . j5pyh2o05c . AQHandles ) {
sdiSetSignalSampleTimeString ( rtDW . j5pyh2o05c . AQHandles , "Continuous" ,
0.0 , ssGetTFinal ( rtS ) ) ; sdiSetSignalRefRate ( rtDW . j5pyh2o05c .
AQHandles , 0.0 ) ; sdiSetRunStartTime ( rtDW . j5pyh2o05c . AQHandles ,
ssGetTaskTime ( rtS , 0 ) ) ; sdiAsyncRepoSetSignalExportSettings ( rtDW .
j5pyh2o05c . AQHandles , 1 , 0 ) ; sdiAsyncRepoSetSignalExportName ( rtDW .
j5pyh2o05c . AQHandles , loggedName , origSigName , propName ) ;
sdiAsyncRepoSetBlockPathDomain ( rtDW . j5pyh2o05c . AQHandles ) ; }
sdiFreeLabel ( sigName ) ; sdiFreeLabel ( loggedName ) ; sdiFreeLabel (
origSigName ) ; sdiFreeLabel ( propName ) ; sdiFreeLabel ( blockPath ) ;
sdiFreeLabel ( blockSID ) ; sdiFreeLabel ( subPath ) ; } } if ( !
isStreamoutAlreadyRegistered ) { { sdiLabelU varName = sdiGetLabelFromChars (
"Losses" ) ; sdiRegisterWksVariable ( rtDW . j5pyh2o05c . AQHandles , varName
, "timeseries" ) ; sdiFreeLabel ( varName ) ; } } } } } { { { bool
isStreamoutAlreadyRegistered = false ; { sdiSignalSourceInfoU srcInfo ;
sdiLabelU loggedName = sdiGetLabelFromChars ( "Manual Switch" ) ; sdiLabelU
origSigName = sdiGetLabelFromChars ( "" ) ; sdiLabelU propName =
sdiGetLabelFromChars ( "Manual Switch" ) ; sdiLabelU blockPath =
sdiGetLabelFromChars ( "trial_with_intermediate_tanks/To Workspace8" ) ;
sdiLabelU blockSID = sdiGetLabelFromChars ( "" ) ; sdiLabelU subPath =
sdiGetLabelFromChars ( "" ) ; sdiDims sigDims ; sdiLabelU sigName =
sdiGetLabelFromChars ( "Manual Switch" ) ; sdiAsyncRepoDataTypeHandle hDT =
sdiAsyncRepoGetBuiltInDataTypeHandle ( DATA_TYPE_DOUBLE ) ; { sdiComplexity
sigComplexity = REAL ; sdiSampleTimeContinuity stCont =
SAMPLE_TIME_CONTINUOUS ; int_T sigDimsArray [ 1 ] = { 1 } ; sigDims . nDims =
1 ; sigDims . dimensions = sigDimsArray ; srcInfo . numBlockPathElems = 1 ;
srcInfo . fullBlockPath = ( sdiFullBlkPathU ) & blockPath ; srcInfo . SID = (
sdiSignalIDU ) & blockSID ; srcInfo . subPath = subPath ; srcInfo . portIndex
= 0 + 1 ; srcInfo . signalName = sigName ; srcInfo . sigSourceUUID = 0 ; rtDW
. p0fxbp34np . AQHandles = sdiStartAsyncioQueueCreation ( hDT , & srcInfo ,
rt_dataMapInfo . mmi . InstanceMap . fullPath ,
"0bfc729f-df12-42b2-9506-8076cd6086bd" , sigComplexity , & sigDims ,
DIMENSIONS_MODE_FIXED , stCont , "" ) ; sdiCompleteAsyncioQueueCreation (
rtDW . p0fxbp34np . AQHandles , hDT , & srcInfo ) ; if ( rtDW . p0fxbp34np .
AQHandles ) { sdiSetSignalSampleTimeString ( rtDW . p0fxbp34np . AQHandles ,
"Continuous" , 0.0 , ssGetTFinal ( rtS ) ) ; sdiSetSignalRefRate ( rtDW .
p0fxbp34np . AQHandles , 0.0 ) ; sdiSetRunStartTime ( rtDW . p0fxbp34np .
AQHandles , ssGetTaskTime ( rtS , 0 ) ) ; sdiAsyncRepoSetSignalExportSettings
( rtDW . p0fxbp34np . AQHandles , 1 , 0 ) ; sdiAsyncRepoSetSignalExportName (
rtDW . p0fxbp34np . AQHandles , loggedName , origSigName , propName ) ;
sdiAsyncRepoSetBlockPathDomain ( rtDW . p0fxbp34np . AQHandles ) ; }
sdiFreeLabel ( sigName ) ; sdiFreeLabel ( loggedName ) ; sdiFreeLabel (
origSigName ) ; sdiFreeLabel ( propName ) ; sdiFreeLabel ( blockPath ) ;
sdiFreeLabel ( blockSID ) ; sdiFreeLabel ( subPath ) ; } } if ( !
isStreamoutAlreadyRegistered ) { { sdiLabelU varName = sdiGetLabelFromChars (
"xdot" ) ; sdiRegisterWksVariable ( rtDW . p0fxbp34np . AQHandles , varName ,
"timeseries" ) ; sdiFreeLabel ( varName ) ; } } } } } { { { bool
isStreamoutAlreadyRegistered = false ; { sdiSignalSourceInfoU srcInfo ;
sdiLabelU loggedName = sdiGetLabelFromChars ( "Displacement" ) ; sdiLabelU
origSigName = sdiGetLabelFromChars ( "Displacement" ) ; sdiLabelU propName =
sdiGetLabelFromChars ( "Displacement" ) ; sdiLabelU blockPath =
sdiGetLabelFromChars ( "trial_with_intermediate_tanks/To Workspace9" ) ;
sdiLabelU blockSID = sdiGetLabelFromChars ( "" ) ; sdiLabelU subPath =
sdiGetLabelFromChars ( "" ) ; sdiDims sigDims ; sdiLabelU sigName =
sdiGetLabelFromChars ( "Displacement" ) ; sdiAsyncRepoDataTypeHandle hDT =
sdiAsyncRepoGetBuiltInDataTypeHandle ( DATA_TYPE_DOUBLE ) ; { sdiComplexity
sigComplexity = REAL ; sdiSampleTimeContinuity stCont =
SAMPLE_TIME_CONTINUOUS ; int_T sigDimsArray [ 1 ] = { 1 } ; sigDims . nDims =
1 ; sigDims . dimensions = sigDimsArray ; srcInfo . numBlockPathElems = 1 ;
srcInfo . fullBlockPath = ( sdiFullBlkPathU ) & blockPath ; srcInfo . SID = (
sdiSignalIDU ) & blockSID ; srcInfo . subPath = subPath ; srcInfo . portIndex
= 0 + 1 ; srcInfo . signalName = sigName ; srcInfo . sigSourceUUID = 0 ; rtDW
. n4qxc3kzkl . AQHandles = sdiStartAsyncioQueueCreation ( hDT , & srcInfo ,
rt_dataMapInfo . mmi . InstanceMap . fullPath ,
"0c0e8ea8-4303-4bcd-b8fd-2358765afa61" , sigComplexity , & sigDims ,
DIMENSIONS_MODE_FIXED , stCont , "" ) ; sdiCompleteAsyncioQueueCreation (
rtDW . n4qxc3kzkl . AQHandles , hDT , & srcInfo ) ; if ( rtDW . n4qxc3kzkl .
AQHandles ) { sdiSetSignalSampleTimeString ( rtDW . n4qxc3kzkl . AQHandles ,
"Continuous" , 0.0 , ssGetTFinal ( rtS ) ) ; sdiSetSignalRefRate ( rtDW .
n4qxc3kzkl . AQHandles , 0.0 ) ; sdiSetRunStartTime ( rtDW . n4qxc3kzkl .
AQHandles , ssGetTaskTime ( rtS , 0 ) ) ; sdiAsyncRepoSetSignalExportSettings
( rtDW . n4qxc3kzkl . AQHandles , 1 , 0 ) ; sdiAsyncRepoSetSignalExportName (
rtDW . n4qxc3kzkl . AQHandles , loggedName , origSigName , propName ) ;
sdiAsyncRepoSetBlockPathDomain ( rtDW . n4qxc3kzkl . AQHandles ) ; }
sdiFreeLabel ( sigName ) ; sdiFreeLabel ( loggedName ) ; sdiFreeLabel (
origSigName ) ; sdiFreeLabel ( propName ) ; sdiFreeLabel ( blockPath ) ;
sdiFreeLabel ( blockSID ) ; sdiFreeLabel ( subPath ) ; } } if ( !
isStreamoutAlreadyRegistered ) { { sdiLabelU varName = sdiGetLabelFromChars (
"x" ) ; sdiRegisterWksVariable ( rtDW . n4qxc3kzkl . AQHandles , varName ,
"timeseries" ) ; sdiFreeLabel ( varName ) ; } } } } } rtB . okpqgwd3h2 = rtP
. InitialSpeed_Value ; rtB . lv3onhl4x1 = rtP . param . P_H ; rtB .
fmn4xw0byt = rtP . param . P_M ; rtB . aeifxkbw2r = rtP . param . P_M ;
MdlInitialize ( ) ; } void MdlOutputs ( int_T tid ) { real_T P0 ; real_T Pg ;
ZCEventType zcEvent ; srClearBC ( rtDW . astl0qyquk ) ; rtB . du3hz2fyzs =
rtX . njnhyztmul ; rtB . frh4wmn3vf = rtP . Gain_Gain * rtB . du3hz2fyzs ; if
( ssIsModeUpdateTimeStep ( rtS ) ) { if ( rtDW . di0xkhax5v != 0 ) { rtX .
cxtidlrvgx = rtB . aeifxkbw2r ; } rtB . mcwb4pv1pd = rtX . cxtidlrvgx ; }
else { rtB . mcwb4pv1pd = rtX . cxtidlrvgx ; } rtB . cqrrgwzxyv = rtB .
mcwb4pv1pd - rtP . param . P_H ; if ( ssIsModeUpdateTimeStep ( rtS ) ) { rtDW
. dfi4io44l1 = ( rtB . cqrrgwzxyv >= 0.0 ) ; } if ( rtDW . dfi4io44l1 > 0 ) {
rtB . dxeimajz1j = rtB . cqrrgwzxyv ; } else { rtB . dxeimajz1j = - rtB .
cqrrgwzxyv ; } if ( ssIsModeUpdateTimeStep ( rtS ) ) { rtDW . e3b1amt0gc = (
( rtB . dxeimajz1j >= rtP . Relay_OnVal ) || ( ( ! ( rtB . dxeimajz1j <= rtP
. Relay_OffVal ) ) && rtDW . e3b1amt0gc ) ) ; } if ( rtDW . e3b1amt0gc ) {
rtB . a5ddnr2k0h = rtP . Relay_YOn ; } else { rtB . a5ddnr2k0h = rtP .
Relay_YOff ; } rtB . blk35wfutl = rtX . anhq5q4adb ; if (
ssIsModeUpdateTimeStep ( rtS ) ) { if ( rtDW . dol3avfpna != 0 ) { rtX .
h1nzd5kpbz = rtB . lv3onhl4x1 ; } rtB . ikogzwrhkz = rtX . h1nzd5kpbz ; if (
rtDW . mfw5n1ea1y != 0 ) { rtX . pmj0r0j4v2 = rtB . fmn4xw0byt ; } rtB .
fakrdhsmg5 = rtX . pmj0r0j4v2 ; } else { rtB . ikogzwrhkz = rtX . h1nzd5kpbz
; rtB . fakrdhsmg5 = rtX . pmj0r0j4v2 ; } rtDW . brdpyuefrq = oswwtj4xei ; Pg
= rtB . ikogzwrhkz - rtP . param . P_L ; P0 = rtP . param . P_H - rtP . param
. P_L ; rtB . dw11js5duv = ( ( muDoubleScalarExp ( Pg / rtP . param . beta )
- ( Pg / rtP . param . beta + 1.0 ) ) * rtP . param . beta - (
muDoubleScalarExp ( P0 / rtP . param . beta ) - ( P0 / rtP . param . beta +
1.0 ) ) * rtP . param . beta ) * rtP . param . V3_0 ; rtDW . bx2j5ck3fy =
oswwtj4xei ; Pg = rtB . fakrdhsmg5 - rtP . param . P_L ; P0 = rtP . param .
P_M - rtP . param . P_L ; rtB . lxsvfwftlj = ( ( muDoubleScalarExp ( Pg / rtP
. param . beta ) - ( Pg / rtP . param . beta + 1.0 ) ) * rtP . param . beta -
( muDoubleScalarExp ( P0 / rtP . param . beta ) - ( P0 / rtP . param . beta +
1.0 ) ) * rtP . param . beta ) * rtP . param . V4_0 ; rtB . cy5jq00nod = rtX
. dusu5owudw ; rtB . ht3vr3fcw1 = rtX . decqh12jsd ; rtB . axh1siq1cq = rtX .
f3sxyupuw4 ; rtB . fkszch0qwe = rtX . javki5lnyj ; rtB . idqoyikyv2 = rtX .
d1tt3gwruv ; rtB . hgyc45hfyg = rtX . erk0rjtwtw ; rtB . ku0vq41nvc = rtX .
m2od42nnln ; Pg = ( ( rtB . fkszch0qwe + rtB . idqoyikyv2 ) + rtB .
hgyc45hfyg ) + rtB . ku0vq41nvc ; if ( ssIsModeUpdateTimeStep ( rtS ) ) { if
( rtDW . oj5q0jsluk != 0 ) { rtX . ith00xbmqs = rtB . okpqgwd3h2 ; } rtB .
cjeggiimao = rtX . ith00xbmqs ; } else { rtB . cjeggiimao = rtX . ith00xbmqs
; } rtB . ab5hqtcfdy = ( rtP . param . J_elec + rtP . param . J_hyd ) * ( rtB
. cjeggiimao * rtB . cjeggiimao ) * rtP . Gain3_Gain ; rtB . ip00peawg3 = ( (
( ( rtB . dw11js5duv + rtB . lxsvfwftlj ) + rtB . ht3vr3fcw1 ) + rtB .
axh1siq1cq ) + Pg ) + rtB . ab5hqtcfdy ; rtB . iubu1mncgz = rtX . djreaptdva
; rtB . na1kfpe1mt = rtB . cbpdlgxiw5 + rtB . iubu1mncgz ; rtB . eee41maadj =
rtX . pnlrwk52wo ; rtB . pebffqxdnd = rtP . TransferFcn_C [ 0 ] * rtX .
lxsclehbdn [ 0 ] ; rtB . pebffqxdnd += rtP . TransferFcn_C [ 1 ] * rtX .
lxsclehbdn [ 1 ] ; if ( ssIsModeUpdateTimeStep ( rtS ) ) { if ( rtB .
pebffqxdnd >= rtP . param . max_Avt ) { rtDW . bapvsmhne1 = 1 ; } else if (
rtB . pebffqxdnd > rtP . Saturation_LowerSat ) { rtDW . bapvsmhne1 = 0 ; }
else { rtDW . bapvsmhne1 = - 1 ; } } if ( rtDW . bapvsmhne1 == 1 ) { rtB .
p2rr4uzpkw = rtP . param . max_Avt ; } else if ( rtDW . bapvsmhne1 == - 1 ) {
rtB . p2rr4uzpkw = rtP . Saturation_LowerSat ; } else { rtB . p2rr4uzpkw =
rtB . pebffqxdnd ; } rtDW . nfkviszx2j = oswwtj4xei ; P0 = rtP . param . P_H
- rtB . ikogzwrhkz ; rtB . fo05goe5wa = rtP . param . Cd * rtB . p2rr4uzpkw *
muDoubleScalarSqrt ( muDoubleScalarAbs ( P0 ) ) * muDoubleScalarSign ( P0 ) ;
rtB . fj3bxcs4ue = rtP . TransferFcn_C_nnurr522qa [ 0 ] * rtX . herd5wji3e [
0 ] ; rtB . fj3bxcs4ue += rtP . TransferFcn_C_nnurr522qa [ 1 ] * rtX .
herd5wji3e [ 1 ] ; if ( ssIsModeUpdateTimeStep ( rtS ) ) { if ( rtB .
fj3bxcs4ue >= rtP . param . max_Avt ) { rtDW . gmzfgv4p33 = 1 ; } else if (
rtB . fj3bxcs4ue > rtP . Saturation_LowerSat_l0ksmwxpmq ) { rtDW . gmzfgv4p33
= 0 ; } else { rtDW . gmzfgv4p33 = - 1 ; } } if ( rtDW . gmzfgv4p33 == 1 ) {
rtB . kkbkkd5z5m = rtP . param . max_Avt ; } else if ( rtDW . gmzfgv4p33 == -
1 ) { rtB . kkbkkd5z5m = rtP . Saturation_LowerSat_l0ksmwxpmq ; } else { rtB
. kkbkkd5z5m = rtB . fj3bxcs4ue ; } rtDW . ojcciswcd3 = oswwtj4xei ; P0 = rtB
. ikogzwrhkz - rtB . mcwb4pv1pd ; rtB . e0pcr2qr33 = rtP . param . Cd * rtB .
kkbkkd5z5m * muDoubleScalarSqrt ( muDoubleScalarAbs ( P0 ) ) *
muDoubleScalarSign ( P0 ) ; rtB . gztmpcdsdg = rtP . TransferFcn_C_j01bq5y01n
[ 0 ] * rtX . ekhpkuwiae [ 0 ] ; rtB . gztmpcdsdg += rtP .
TransferFcn_C_j01bq5y01n [ 1 ] * rtX . ekhpkuwiae [ 1 ] ; if (
ssIsModeUpdateTimeStep ( rtS ) ) { if ( rtB . gztmpcdsdg >= rtP . param .
max_Avt ) { rtDW . bxfqhgcloy = 1 ; } else if ( rtB . gztmpcdsdg > rtP .
Saturation_LowerSat_m0gx2ayxie ) { rtDW . bxfqhgcloy = 0 ; } else { rtDW .
bxfqhgcloy = - 1 ; } } if ( rtDW . bxfqhgcloy == 1 ) { rtB . mxnffxyshc = rtP
. param . max_Avt ; } else if ( rtDW . bxfqhgcloy == - 1 ) { rtB . mxnffxyshc
= rtP . Saturation_LowerSat_m0gx2ayxie ; } else { rtB . mxnffxyshc = rtB .
gztmpcdsdg ; } rtDW . brka54cpul = oswwtj4xei ; P0 = rtB . fakrdhsmg5 - rtB .
mcwb4pv1pd ; rtB . haqtqgqlej = rtP . param . Cd * rtB . mxnffxyshc *
muDoubleScalarSqrt ( muDoubleScalarAbs ( P0 ) ) * muDoubleScalarSign ( P0 ) ;
rtB . fr14joswus = rtP . TransferFcn_C_kavhf31tp1 [ 0 ] * rtX . k2eyyio13c [
0 ] ; rtB . fr14joswus += rtP . TransferFcn_C_kavhf31tp1 [ 1 ] * rtX .
k2eyyio13c [ 1 ] ; if ( ssIsModeUpdateTimeStep ( rtS ) ) { if ( rtB .
fr14joswus >= rtP . param . max_Avt ) { rtDW . gwwzsdaztz = 1 ; } else if (
rtB . fr14joswus > rtP . Saturation_LowerSat_l5sm1vdo10 ) { rtDW . gwwzsdaztz
= 0 ; } else { rtDW . gwwzsdaztz = - 1 ; } } if ( rtDW . gwwzsdaztz == 1 ) {
rtB . g3vckbtrp3 = rtP . param . max_Avt ; } else if ( rtDW . gwwzsdaztz == -
1 ) { rtB . g3vckbtrp3 = rtP . Saturation_LowerSat_l5sm1vdo10 ; } else { rtB
. g3vckbtrp3 = rtB . fr14joswus ; } rtDW . d1ijg3op2k = oswwtj4xei ; P0 = rtP
. param . P_H - rtB . fakrdhsmg5 ; rtB . idpdlwl5kt = rtP . param . Cd * rtB
. g3vckbtrp3 * muDoubleScalarSqrt ( muDoubleScalarAbs ( P0 ) ) *
muDoubleScalarSign ( P0 ) ; rtDW . kjn3oyfrac = oswwtj4xei ; rtB . pwpjs2gjry
= rtP . param . D / 6.2831853071795862 * rtB . cjeggiimao ; rtB . cqkdaw1iu3
= ssGetT ( rtS ) ; if ( ssIsSampleHit ( rtS , 1 , 0 ) &&
ssIsModeUpdateTimeStep ( rtS ) ) { zcEvent = rt_ZCFcn ( RISING_ZERO_CROSSING
, & rtPrevZCX . ii2vsk3a4a , ( rtB . a5ddnr2k0h ) ) ; if ( zcEvent !=
NO_ZCEVENT ) { rtB . a1drd1cazf = rtB . cqkdaw1iu3 ; rtDW . astl0qyquk = 4 ;
} } rtB . o123mdetox = rtX . mh2ntyj1vq ; { if ( rtDW . cpa2ih4ymh .
AQHandles && ssGetLogOutput ( rtS ) ) { sdiWriteSignal ( rtDW . cpa2ih4ymh .
AQHandles , ssGetTaskTime ( rtS , 0 ) , ( char * ) & rtB . mcwb4pv1pd + 0 ) ;
} } { if ( rtDW . igylchmu0t . AQHandles && ssGetLogOutput ( rtS ) ) {
sdiWriteSignal ( rtDW . igylchmu0t . AQHandles , ssGetTaskTime ( rtS , 0 ) ,
( char * ) & rtB . ht3vr3fcw1 + 0 ) ; } } if ( ssIsSampleHit ( rtS , 2 , 0 )
) { rtB . jsgd1kggkx = rtDW . mog1ypd0na ; } rtDW . d3q0cyta13 = oswwtj4xei ;
if ( rtB . a5ddnr2k0h == 1.0 ) { if ( rtB . cjeggiimao < 200.0 ) { rtB .
b3cdqkxcij = rtB . jsgd1kggkx ; } else { rtB . b3cdqkxcij = 5.0 ; } } else if
( rtB . cjeggiimao >= 300.0 ) { rtB . b3cdqkxcij = rtB . jsgd1kggkx ; } else
{ rtB . b3cdqkxcij = - 5.0 ; } { if ( rtDW . aqjm4jt1o5 . AQHandles &&
ssGetLogOutput ( rtS ) ) { sdiWriteSignal ( rtDW . aqjm4jt1o5 . AQHandles ,
ssGetTaskTime ( rtS , 0 ) , ( char * ) & rtB . b3cdqkxcij + 0 ) ; } } { if (
rtDW . gzhbrbv4c4 . AQHandles && ssGetLogOutput ( rtS ) ) { sdiWriteSignal (
rtDW . gzhbrbv4c4 . AQHandles , ssGetTaskTime ( rtS , 0 ) , ( char * ) & rtB
. iubu1mncgz + 0 ) ; } } { if ( rtDW . htile1lzbu . AQHandles &&
ssGetLogOutput ( rtS ) ) { sdiWriteSignal ( rtDW . htile1lzbu . AQHandles ,
ssGetTaskTime ( rtS , 0 ) , ( char * ) & rtB . ab5hqtcfdy + 0 ) ; } } if (
ssIsSampleHit ( rtS , 1 , 0 ) ) { { if ( rtDW . p1ayluai4b . AQHandles &&
ssGetLogOutput ( rtS ) ) { sdiWriteSignal ( rtDW . p1ayluai4b . AQHandles ,
ssGetTaskTime ( rtS , 1 ) , ( char * ) & rtB . a1drd1cazf + 0 ) ; } } } { if
( rtDW . gjtuyi3wgw . AQHandles && ssGetLogOutput ( rtS ) ) { sdiWriteSignal
( rtDW . gjtuyi3wgw . AQHandles , ssGetTaskTime ( rtS , 0 ) , ( char * ) &
rtB . cjeggiimao + 0 ) ; } } rtDW . ncj5ml0ihc = oswwtj4xei ; P0 = rtP .
param . D / 6.2831853071795862 * ( rtB . ikogzwrhkz - rtB . fakrdhsmg5 ) ;
rtB . iikt2yfv4y = 1.0 / ( rtP . param . J_elec + rtP . param . J_hyd ) * (
P0 - rtB . b3cdqkxcij ) ; rtB . ilwj52njki = P0 ; rtB . nc3ujqp1pg = rtB .
iikt2yfv4y * rtB . cjeggiimao ; rtB . iilt4zdhq4 = ( rtP . param . J_elec +
rtP . param . J_hyd ) * rtB . nc3ujqp1pg ; { if ( rtDW . jd4fubbir0 .
AQHandles && ssGetLogOutput ( rtS ) ) { sdiWriteSignal ( rtDW . jd4fubbir0 .
AQHandles , ssGetTaskTime ( rtS , 0 ) , ( char * ) & rtB . iilt4zdhq4 + 0 ) ;
} } rtB . foabxlso54 = rtB . b3cdqkxcij * rtB . cjeggiimao ; { if ( rtDW .
kbuynfluan . AQHandles && ssGetLogOutput ( rtS ) ) { sdiWriteSignal ( rtDW .
kbuynfluan . AQHandles , ssGetTaskTime ( rtS , 0 ) , ( char * ) & rtB .
foabxlso54 + 0 ) ; } } { if ( rtDW . j5pyh2o05c . AQHandles && ssGetLogOutput
( rtS ) ) { sdiWriteSignal ( rtDW . j5pyh2o05c . AQHandles , ssGetTaskTime (
rtS , 0 ) , ( char * ) & Pg + 0 ) ; } } if ( rtP .
ManualSwitch_CurrentSetting == 1 ) { rtB . mscp3mtqiu = muDoubleScalarSin (
rtP . Velocity_Freq * ssGetTaskTime ( rtS , 0 ) + rtP . Velocity_Phase ) *
rtP . Velocity_Amp + rtP . Velocity_Bias ; } else { rtB . mscp3mtqiu = rtP .
Constant_Value ; } { if ( rtDW . p0fxbp34np . AQHandles && ssGetLogOutput (
rtS ) ) { sdiWriteSignal ( rtDW . p0fxbp34np . AQHandles , ssGetTaskTime (
rtS , 0 ) , ( char * ) & rtB . mscp3mtqiu + 0 ) ; } } rtB . orhxqsmvfm = rtX
. kh4uwqwwas ; { if ( rtDW . n4qxc3kzkl . AQHandles && ssGetLogOutput ( rtS )
) { sdiWriteSignal ( rtDW . n4qxc3kzkl . AQHandles , ssGetTaskTime ( rtS , 0
) , ( char * ) & rtB . orhxqsmvfm + 0 ) ; } } rtB . pcceognrbb = rtB .
e0pcr2qr33 + rtB . haqtqgqlej ; rtB . hcfgbzdbcf = rtB . fo05goe5wa + rtB .
idpdlwl5kt ; rtDW . esun0jx3xw = oswwtj4xei ; rtB . jabohzh110 = rtP . param
. beta / ( rtP . param . Acap * rtB . orhxqsmvfm + rtP . param . V1_0 ) * ( (
rtB . e0pcr2qr33 + rtB . haqtqgqlej ) - rtP . param . Acap * rtB . mscp3mtqiu
) ; if ( ssIsSampleHit ( rtS , 1 , 0 ) ) { rtDW . fzispyfeqr = oswwtj4xei ;
if ( rtB . a5ddnr2k0h == 1.0 ) { rtB . jjwnm1umsr = rtP . param . max_Avt ;
rtB . jbbgrda5ov = rtP . param . max_Avt ; rtB . bojkujnbsq = 0.0 ; rtB .
bnvm1inada = rtP . param . max_Avt ; } else { rtB . jjwnm1umsr = rtP . param
. max_Avt ; rtB . jbbgrda5ov = 0.0 ; rtB . bojkujnbsq = rtP . param . max_Avt
; rtB . bnvm1inada = 0.0 ; } } rtDW . puivctwyq5 = oswwtj4xei ; Pg = rtB .
mcwb4pv1pd - rtP . param . P_L ; rtB . msv1wjj5cz = ( ( muDoubleScalarExp (
Pg / rtP . param . beta ) - ( Pg / rtP . param . beta + 1.0 ) ) * rtP . param
. beta + Pg ) * ( rtB . e0pcr2qr33 + rtB . haqtqgqlej ) ; rtDW . egiikq3opb =
oswwtj4xei ; rtB . fzo2unh14f = ( muDoubleScalarExp ( ( rtB . ikogzwrhkz -
rtP . param . P_L ) / rtP . param . beta ) - 1.0 ) * rtP . param . beta * ( (
rtB . fo05goe5wa - rtB . e0pcr2qr33 ) - rtB . pwpjs2gjry ) ; rtDW .
mevs4smziy = oswwtj4xei ; rtB . ollqjfiay1 = ( muDoubleScalarExp ( ( rtB .
fakrdhsmg5 - rtP . param . P_L ) / rtP . param . beta ) - 1.0 ) * rtP . param
. beta * ( ( rtB . pwpjs2gjry + rtB . idpdlwl5kt ) - rtB . haqtqgqlej ) ; rtB
. jost0dydiu = rtB . mcwb4pv1pd - rtP . param . P_L ; rtB . nodqar1hqn = rtB
. pcceognrbb * rtB . jost0dydiu ; rtB . oi30omicce = rtB . hcfgbzdbcf * rtP .
param . P_H ; rtB . nmkcufy13s = rtB . ikogzwrhkz - rtB . mcwb4pv1pd ; rtB .
bhydlzvdgu = rtB . nmkcufy13s * rtB . e0pcr2qr33 ; rtB . jnk41mcgi4 = rtB .
fakrdhsmg5 - rtB . mcwb4pv1pd ; rtB . nhtjb14mqm = rtB . jnk41mcgi4 * rtB .
haqtqgqlej ; rtB . becosfp3hi = rtP . param . P_H - rtB . ikogzwrhkz ; rtB .
lphlstzov0 = rtB . becosfp3hi * rtB . fo05goe5wa ; rtB . dcpeged5b2 = rtP .
param . P_H - rtB . fakrdhsmg5 ; rtB . eqjfsu4kur = rtB . dcpeged5b2 * rtB .
idpdlwl5kt ; rtDW . c2jtsctdjh = oswwtj4xei ; rtB . hrpirqqfri = ( ( rtB .
fo05goe5wa - rtB . e0pcr2qr33 ) - rtB . pwpjs2gjry ) * ( rtP . param . beta /
rtP . param . V3_0 ) ; rtDW . ndyrr3tae0 = oswwtj4xei ; rtB . lqpliiebzv = (
( rtB . idpdlwl5kt + rtB . pwpjs2gjry ) - rtB . haqtqgqlej ) * ( rtP . param
. beta / rtP . param . V4_0 ) ; rtB . cpcl24qdbb = rtP .
TransferFcn_C_lu1wzht1on [ 0 ] * rtX . l1m0lcksjp [ 0 ] ; rtB . cpcl24qdbb +=
rtP . TransferFcn_C_lu1wzht1on [ 1 ] * rtX . l1m0lcksjp [ 1 ] ; if (
ssIsModeUpdateTimeStep ( rtS ) ) { if ( rtB . cpcl24qdbb >= rtP . param .
max_Avt ) { rtDW . agicmdjqdi = 1 ; } else if ( rtB . cpcl24qdbb > rtP .
Saturation_LowerSat_mclkng5niv ) { rtDW . agicmdjqdi = 0 ; } else { rtDW .
agicmdjqdi = - 1 ; } } if ( rtDW . agicmdjqdi == 1 ) { rtB . hpzll4vow0 = rtP
. param . max_Avt ; } else if ( rtDW . agicmdjqdi == - 1 ) { rtB . hpzll4vow0
= rtP . Saturation_LowerSat_mclkng5niv ; } else { rtB . hpzll4vow0 = rtB .
cpcl24qdbb ; } UNUSED_PARAMETER ( tid ) ; } void MdlOutputsTID3 ( int_T tid )
{ rtB . okpqgwd3h2 = rtP . InitialSpeed_Value ; rtB . cbpdlgxiw5 = ( rtP .
param . J_elec + rtP . param . J_hyd ) * ( rtB . okpqgwd3h2 * rtB .
okpqgwd3h2 ) * rtP . Gain5_Gain ; rtB . lv3onhl4x1 = rtP . param . P_H ; rtB
. fmn4xw0byt = rtP . param . P_M ; rtB . aeifxkbw2r = rtP . param . P_M ;
UNUSED_PARAMETER ( tid ) ; } void MdlUpdate ( int_T tid ) { rtDW . di0xkhax5v
= 0 ; rtDW . dol3avfpna = 0 ; rtDW . mfw5n1ea1y = 0 ; rtDW . oj5q0jsluk = 0 ;
if ( ssIsSampleHit ( rtS , 2 , 0 ) ) { rtDW . mog1ypd0na = rtB . ilwj52njki ;
} UNUSED_PARAMETER ( tid ) ; } void MdlUpdateTID3 ( int_T tid ) {
UNUSED_PARAMETER ( tid ) ; } void MdlDerivatives ( void ) { XDot * _rtXdot ;
_rtXdot = ( ( XDot * ) ssGetdX ( rtS ) ) ; _rtXdot -> njnhyztmul = rtB .
cjeggiimao ; _rtXdot -> cxtidlrvgx = rtB . jabohzh110 ; _rtXdot -> anhq5q4adb
= rtB . ollqjfiay1 ; _rtXdot -> h1nzd5kpbz = rtB . hrpirqqfri ; _rtXdot ->
pmj0r0j4v2 = rtB . lqpliiebzv ; _rtXdot -> dusu5owudw = rtB . fzo2unh14f ;
_rtXdot -> decqh12jsd = rtB . foabxlso54 ; _rtXdot -> f3sxyupuw4 = rtB .
msv1wjj5cz ; _rtXdot -> javki5lnyj = rtB . lphlstzov0 ; _rtXdot -> d1tt3gwruv
= rtB . bhydlzvdgu ; _rtXdot -> erk0rjtwtw = rtB . nhtjb14mqm ; _rtXdot ->
m2od42nnln = rtB . eqjfsu4kur ; _rtXdot -> ith00xbmqs = rtB . iikt2yfv4y ;
_rtXdot -> djreaptdva = rtB . oi30omicce ; _rtXdot -> pnlrwk52wo = rtB .
iilt4zdhq4 ; _rtXdot -> lxsclehbdn [ 0 ] = rtP . TransferFcn_A [ 0 ] * rtX .
lxsclehbdn [ 0 ] ; _rtXdot -> lxsclehbdn [ 0 ] += rtP . TransferFcn_A [ 1 ] *
rtX . lxsclehbdn [ 1 ] ; _rtXdot -> lxsclehbdn [ 1 ] = rtX . lxsclehbdn [ 0 ]
; _rtXdot -> lxsclehbdn [ 0 ] += rtB . jjwnm1umsr ; _rtXdot -> herd5wji3e [ 0
] = rtP . TransferFcn_A_msxihde20w [ 0 ] * rtX . herd5wji3e [ 0 ] ; _rtXdot
-> herd5wji3e [ 0 ] += rtP . TransferFcn_A_msxihde20w [ 1 ] * rtX .
herd5wji3e [ 1 ] ; _rtXdot -> herd5wji3e [ 1 ] = rtX . herd5wji3e [ 0 ] ;
_rtXdot -> herd5wji3e [ 0 ] += rtB . jbbgrda5ov ; _rtXdot -> ekhpkuwiae [ 0 ]
= rtP . TransferFcn_A_bwuzoqq3iu [ 0 ] * rtX . ekhpkuwiae [ 0 ] ; _rtXdot ->
ekhpkuwiae [ 0 ] += rtP . TransferFcn_A_bwuzoqq3iu [ 1 ] * rtX . ekhpkuwiae [
1 ] ; _rtXdot -> ekhpkuwiae [ 1 ] = rtX . ekhpkuwiae [ 0 ] ; _rtXdot ->
ekhpkuwiae [ 0 ] += rtB . bojkujnbsq ; _rtXdot -> k2eyyio13c [ 0 ] = rtP .
TransferFcn_A_cqt2tjfue3 [ 0 ] * rtX . k2eyyio13c [ 0 ] ; _rtXdot ->
k2eyyio13c [ 0 ] += rtP . TransferFcn_A_cqt2tjfue3 [ 1 ] * rtX . k2eyyio13c [
1 ] ; _rtXdot -> k2eyyio13c [ 1 ] = rtX . k2eyyio13c [ 0 ] ; _rtXdot ->
k2eyyio13c [ 0 ] += rtB . bnvm1inada ; _rtXdot -> mh2ntyj1vq = rtB .
nodqar1hqn ; _rtXdot -> kh4uwqwwas = rtB . mscp3mtqiu ; _rtXdot -> l1m0lcksjp
[ 0 ] = rtP . TransferFcn_A_m2qs34lzon [ 0 ] * rtX . l1m0lcksjp [ 0 ] ;
_rtXdot -> l1m0lcksjp [ 0 ] += rtP . TransferFcn_A_m2qs34lzon [ 1 ] * rtX .
l1m0lcksjp [ 1 ] ; _rtXdot -> l1m0lcksjp [ 1 ] = rtX . l1m0lcksjp [ 0 ] ;
_rtXdot -> l1m0lcksjp [ 0 ] += rtP . param . max_Avt ; } void MdlProjection (
void ) { } void MdlZeroCrossings ( void ) { ZCV * _rtZCSV ; _rtZCSV = ( ( ZCV
* ) ssGetSolverZcSignalVector ( rtS ) ) ; _rtZCSV -> n52axficl2 = rtB .
cqrrgwzxyv ; if ( rtDW . e3b1amt0gc ) { _rtZCSV -> a15z3evljq = rtB .
dxeimajz1j - rtP . Relay_OffVal ; } else { _rtZCSV -> a15z3evljq = rtB .
dxeimajz1j - rtP . Relay_OnVal ; } _rtZCSV -> fbpmdji415 = rtB . pebffqxdnd -
rtP . param . max_Avt ; _rtZCSV -> odfwjcwoxs = rtB . pebffqxdnd - rtP .
Saturation_LowerSat ; _rtZCSV -> fdgm22ouln = rtB . fj3bxcs4ue - rtP . param
. max_Avt ; _rtZCSV -> haaywk4krz = rtB . fj3bxcs4ue - rtP .
Saturation_LowerSat_l0ksmwxpmq ; _rtZCSV -> l2jklhprzm = rtB . gztmpcdsdg -
rtP . param . max_Avt ; _rtZCSV -> jdcwcq0rv4 = rtB . gztmpcdsdg - rtP .
Saturation_LowerSat_m0gx2ayxie ; _rtZCSV -> dqejsrbctn = rtB . fr14joswus -
rtP . param . max_Avt ; _rtZCSV -> fs30lsvp00 = rtB . fr14joswus - rtP .
Saturation_LowerSat_l5sm1vdo10 ; _rtZCSV -> jy5xmnlgzp = rtB . cpcl24qdbb -
rtP . param . max_Avt ; _rtZCSV -> daixzoq4vk = rtB . cpcl24qdbb - rtP .
Saturation_LowerSat_mclkng5niv ; } void MdlTerminate ( void ) { { if ( rtDW .
cpa2ih4ymh . AQHandles ) { sdiTerminateStreaming ( & rtDW . cpa2ih4ymh .
AQHandles ) ; } } { if ( rtDW . igylchmu0t . AQHandles ) {
sdiTerminateStreaming ( & rtDW . igylchmu0t . AQHandles ) ; } } { if ( rtDW .
aqjm4jt1o5 . AQHandles ) { sdiTerminateStreaming ( & rtDW . aqjm4jt1o5 .
AQHandles ) ; } } { if ( rtDW . gzhbrbv4c4 . AQHandles ) {
sdiTerminateStreaming ( & rtDW . gzhbrbv4c4 . AQHandles ) ; } } { if ( rtDW .
htile1lzbu . AQHandles ) { sdiTerminateStreaming ( & rtDW . htile1lzbu .
AQHandles ) ; } } { if ( rtDW . p1ayluai4b . AQHandles ) {
sdiTerminateStreaming ( & rtDW . p1ayluai4b . AQHandles ) ; } } { if ( rtDW .
gjtuyi3wgw . AQHandles ) { sdiTerminateStreaming ( & rtDW . gjtuyi3wgw .
AQHandles ) ; } } { if ( rtDW . jd4fubbir0 . AQHandles ) {
sdiTerminateStreaming ( & rtDW . jd4fubbir0 . AQHandles ) ; } } { if ( rtDW .
kbuynfluan . AQHandles ) { sdiTerminateStreaming ( & rtDW . kbuynfluan .
AQHandles ) ; } } { if ( rtDW . j5pyh2o05c . AQHandles ) {
sdiTerminateStreaming ( & rtDW . j5pyh2o05c . AQHandles ) ; } } { if ( rtDW .
p0fxbp34np . AQHandles ) { sdiTerminateStreaming ( & rtDW . p0fxbp34np .
AQHandles ) ; } } { if ( rtDW . n4qxc3kzkl . AQHandles ) {
sdiTerminateStreaming ( & rtDW . n4qxc3kzkl . AQHandles ) ; } } } static void
mr_trial_with_intermediate_tanks_cacheDataAsMxArray ( mxArray * destArray ,
mwIndex i , int j , const void * srcData , size_t numBytes ) ; static void
mr_trial_with_intermediate_tanks_cacheDataAsMxArray ( mxArray * destArray ,
mwIndex i , int j , const void * srcData , size_t numBytes ) { mxArray *
newArray = mxCreateUninitNumericMatrix ( ( size_t ) 1 , numBytes ,
mxUINT8_CLASS , mxREAL ) ; memcpy ( ( uint8_T * ) mxGetData ( newArray ) , (
const uint8_T * ) srcData , numBytes ) ; mxSetFieldByNumber ( destArray , i ,
j , newArray ) ; } static void
mr_trial_with_intermediate_tanks_restoreDataFromMxArray ( void * destData ,
const mxArray * srcArray , mwIndex i , int j , size_t numBytes ) ; static
void mr_trial_with_intermediate_tanks_restoreDataFromMxArray ( void *
destData , const mxArray * srcArray , mwIndex i , int j , size_t numBytes ) {
memcpy ( ( uint8_T * ) destData , ( const uint8_T * ) mxGetData (
mxGetFieldByNumber ( srcArray , i , j ) ) , numBytes ) ; } static void
mr_trial_with_intermediate_tanks_cacheBitFieldToMxArray ( mxArray * destArray
, mwIndex i , int j , uint_T bitVal ) ; static void
mr_trial_with_intermediate_tanks_cacheBitFieldToMxArray ( mxArray * destArray
, mwIndex i , int j , uint_T bitVal ) { mxSetFieldByNumber ( destArray , i ,
j , mxCreateDoubleScalar ( ( real_T ) bitVal ) ) ; } static uint_T
mr_trial_with_intermediate_tanks_extractBitFieldFromMxArray ( const mxArray *
srcArray , mwIndex i , int j , uint_T numBits ) ; static uint_T
mr_trial_with_intermediate_tanks_extractBitFieldFromMxArray ( const mxArray *
srcArray , mwIndex i , int j , uint_T numBits ) { const uint_T varVal = (
uint_T ) mxGetScalar ( mxGetFieldByNumber ( srcArray , i , j ) ) ; return
varVal & ( ( 1u << numBits ) - 1u ) ; } static void
mr_trial_with_intermediate_tanks_cacheDataToMxArrayWithOffset ( mxArray *
destArray , mwIndex i , int j , mwIndex offset , const void * srcData ,
size_t numBytes ) ; static void
mr_trial_with_intermediate_tanks_cacheDataToMxArrayWithOffset ( mxArray *
destArray , mwIndex i , int j , mwIndex offset , const void * srcData ,
size_t numBytes ) { uint8_T * varData = ( uint8_T * ) mxGetData (
mxGetFieldByNumber ( destArray , i , j ) ) ; memcpy ( ( uint8_T * ) & varData
[ offset * numBytes ] , ( const uint8_T * ) srcData , numBytes ) ; } static
void mr_trial_with_intermediate_tanks_restoreDataFromMxArrayWithOffset ( void
* destData , const mxArray * srcArray , mwIndex i , int j , mwIndex offset ,
size_t numBytes ) ; static void
mr_trial_with_intermediate_tanks_restoreDataFromMxArrayWithOffset ( void *
destData , const mxArray * srcArray , mwIndex i , int j , mwIndex offset ,
size_t numBytes ) { const uint8_T * varData = ( const uint8_T * ) mxGetData (
mxGetFieldByNumber ( srcArray , i , j ) ) ; memcpy ( ( uint8_T * ) destData ,
( const uint8_T * ) & varData [ offset * numBytes ] , numBytes ) ; } static
void mr_trial_with_intermediate_tanks_cacheBitFieldToCellArrayWithOffset (
mxArray * destArray , mwIndex i , int j , mwIndex offset , uint_T fieldVal )
; static void
mr_trial_with_intermediate_tanks_cacheBitFieldToCellArrayWithOffset ( mxArray
* destArray , mwIndex i , int j , mwIndex offset , uint_T fieldVal ) {
mxSetCell ( mxGetFieldByNumber ( destArray , i , j ) , offset ,
mxCreateDoubleScalar ( ( real_T ) fieldVal ) ) ; } static uint_T
mr_trial_with_intermediate_tanks_extractBitFieldFromCellArrayWithOffset (
const mxArray * srcArray , mwIndex i , int j , mwIndex offset , uint_T
numBits ) ; static uint_T
mr_trial_with_intermediate_tanks_extractBitFieldFromCellArrayWithOffset (
const mxArray * srcArray , mwIndex i , int j , mwIndex offset , uint_T
numBits ) { const uint_T fieldVal = ( uint_T ) mxGetScalar ( mxGetCell (
mxGetFieldByNumber ( srcArray , i , j ) , offset ) ) ; return fieldVal & ( (
1u << numBits ) - 1u ) ; } mxArray *
mr_trial_with_intermediate_tanks_GetDWork ( ) { static const char_T *
ssDWFieldNames [ 3 ] = { "rtB" , "rtDW" , "rtPrevZCX" , } ; mxArray * ssDW =
mxCreateStructMatrix ( 1 , 1 , 3 , ssDWFieldNames ) ;
mr_trial_with_intermediate_tanks_cacheDataAsMxArray ( ssDW , 0 , 0 , ( const
void * ) & ( rtB ) , sizeof ( rtB ) ) ; { static const char_T *
rtdwDataFieldNames [ 61 ] = { "rtDW.mog1ypd0na" , "rtDW.d1ijg3op2k" ,
"rtDW.brka54cpul" , "rtDW.ojcciswcd3" , "rtDW.nfkviszx2j" , "rtDW.ndyrr3tae0"
, "rtDW.c2jtsctdjh" , "rtDW.mevs4smziy" , "rtDW.egiikq3opb" ,
"rtDW.puivctwyq5" , "rtDW.bx2j5ck3fy" , "rtDW.fzispyfeqr" , "rtDW.brdpyuefrq"
, "rtDW.esun0jx3xw" , "rtDW.ncj5ml0ihc" , "rtDW.kjn3oyfrac" ,
"rtDW.d3q0cyta13" , "rtDW.di0xkhax5v" , "rtDW.dol3avfpna" , "rtDW.mfw5n1ea1y"
, "rtDW.oj5q0jsluk" , "rtDW.dfi4io44l1" , "rtDW.bapvsmhne1" ,
"rtDW.gmzfgv4p33" , "rtDW.bxfqhgcloy" , "rtDW.gwwzsdaztz" , "rtDW.agicmdjqdi"
, "rtDW.astl0qyquk" , "rtDW.ekvnfwjtsb" , "rtDW.dql53suw4t" ,
"rtDW.gz0frs1wku" , "rtDW.iwrzpan53u" , "rtDW.ivbuycouap" , "rtDW.oatrzraqor"
, "rtDW.h3yqjmxcau" , "rtDW.c2gp1hcawr" , "rtDW.kzezuc5kmj" ,
"rtDW.baej0t4mz4" , "rtDW.e4o1f3himj" , "rtDW.etn0cabray" , "rtDW.nq5efgjr4u"
, "rtDW.ldwj0svuv5" , "rtDW.et5jt4roea" , "rtDW.jhgzfejyh1" ,
"rtDW.e3b1amt0gc" , "rtDW.g3swjv3kh4" , "rtDW.cifi3jtean" , "rtDW.mx1mzmmudm"
, "rtDW.nhmw0vweas" , "rtDW.enb1yfnvht" , "rtDW.au1n1jxucn" ,
"rtDW.cwu3pvv3pv" , "rtDW.dwzumrr4bn" , "rtDW.jelszbraal" , "rtDW.cbdxrybiau"
, "rtDW.dia2paiwqh" , "rtDW.f40yghkee4" , "rtDW.pke0ukp1xy" ,
"rtDW.ob1pklob1o" , "rtDW.jvf3pz3ask" , "rtDW.cy4ormxj51" , } ; mxArray *
rtdwData = mxCreateStructMatrix ( 1 , 1 , 61 , rtdwDataFieldNames ) ;
mr_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 , 0 , (
const void * ) & ( rtDW . mog1ypd0na ) , sizeof ( rtDW . mog1ypd0na ) ) ;
mr_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 , 1 , (
const void * ) & ( rtDW . d1ijg3op2k ) , sizeof ( rtDW . d1ijg3op2k ) ) ;
mr_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 , 2 , (
const void * ) & ( rtDW . brka54cpul ) , sizeof ( rtDW . brka54cpul ) ) ;
mr_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 , 3 , (
const void * ) & ( rtDW . ojcciswcd3 ) , sizeof ( rtDW . ojcciswcd3 ) ) ;
mr_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 , 4 , (
const void * ) & ( rtDW . nfkviszx2j ) , sizeof ( rtDW . nfkviszx2j ) ) ;
mr_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 , 5 , (
const void * ) & ( rtDW . ndyrr3tae0 ) , sizeof ( rtDW . ndyrr3tae0 ) ) ;
mr_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 , 6 , (
const void * ) & ( rtDW . c2jtsctdjh ) , sizeof ( rtDW . c2jtsctdjh ) ) ;
mr_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 , 7 , (
const void * ) & ( rtDW . mevs4smziy ) , sizeof ( rtDW . mevs4smziy ) ) ;
mr_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 , 8 , (
const void * ) & ( rtDW . egiikq3opb ) , sizeof ( rtDW . egiikq3opb ) ) ;
mr_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 , 9 , (
const void * ) & ( rtDW . puivctwyq5 ) , sizeof ( rtDW . puivctwyq5 ) ) ;
mr_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 , 10 , (
const void * ) & ( rtDW . bx2j5ck3fy ) , sizeof ( rtDW . bx2j5ck3fy ) ) ;
mr_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 , 11 , (
const void * ) & ( rtDW . fzispyfeqr ) , sizeof ( rtDW . fzispyfeqr ) ) ;
mr_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 , 12 , (
const void * ) & ( rtDW . brdpyuefrq ) , sizeof ( rtDW . brdpyuefrq ) ) ;
mr_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 , 13 , (
const void * ) & ( rtDW . esun0jx3xw ) , sizeof ( rtDW . esun0jx3xw ) ) ;
mr_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 , 14 , (
const void * ) & ( rtDW . ncj5ml0ihc ) , sizeof ( rtDW . ncj5ml0ihc ) ) ;
mr_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 , 15 , (
const void * ) & ( rtDW . kjn3oyfrac ) , sizeof ( rtDW . kjn3oyfrac ) ) ;
mr_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 , 16 , (
const void * ) & ( rtDW . d3q0cyta13 ) , sizeof ( rtDW . d3q0cyta13 ) ) ;
mr_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 , 17 , (
const void * ) & ( rtDW . di0xkhax5v ) , sizeof ( rtDW . di0xkhax5v ) ) ;
mr_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 , 18 , (
const void * ) & ( rtDW . dol3avfpna ) , sizeof ( rtDW . dol3avfpna ) ) ;
mr_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 , 19 , (
const void * ) & ( rtDW . mfw5n1ea1y ) , sizeof ( rtDW . mfw5n1ea1y ) ) ;
mr_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 , 20 , (
const void * ) & ( rtDW . oj5q0jsluk ) , sizeof ( rtDW . oj5q0jsluk ) ) ;
mr_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 , 21 , (
const void * ) & ( rtDW . dfi4io44l1 ) , sizeof ( rtDW . dfi4io44l1 ) ) ;
mr_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 , 22 , (
const void * ) & ( rtDW . bapvsmhne1 ) , sizeof ( rtDW . bapvsmhne1 ) ) ;
mr_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 , 23 , (
const void * ) & ( rtDW . gmzfgv4p33 ) , sizeof ( rtDW . gmzfgv4p33 ) ) ;
mr_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 , 24 , (
const void * ) & ( rtDW . bxfqhgcloy ) , sizeof ( rtDW . bxfqhgcloy ) ) ;
mr_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 , 25 , (
const void * ) & ( rtDW . gwwzsdaztz ) , sizeof ( rtDW . gwwzsdaztz ) ) ;
mr_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 , 26 , (
const void * ) & ( rtDW . agicmdjqdi ) , sizeof ( rtDW . agicmdjqdi ) ) ;
mr_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 , 27 , (
const void * ) & ( rtDW . astl0qyquk ) , sizeof ( rtDW . astl0qyquk ) ) ;
mr_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 , 28 , (
const void * ) & ( rtDW . ekvnfwjtsb ) , sizeof ( rtDW . ekvnfwjtsb ) ) ;
mr_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 , 29 , (
const void * ) & ( rtDW . dql53suw4t ) , sizeof ( rtDW . dql53suw4t ) ) ;
mr_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 , 30 , (
const void * ) & ( rtDW . gz0frs1wku ) , sizeof ( rtDW . gz0frs1wku ) ) ;
mr_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 , 31 , (
const void * ) & ( rtDW . iwrzpan53u ) , sizeof ( rtDW . iwrzpan53u ) ) ;
mr_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 , 32 , (
const void * ) & ( rtDW . ivbuycouap ) , sizeof ( rtDW . ivbuycouap ) ) ;
mr_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 , 33 , (
const void * ) & ( rtDW . oatrzraqor ) , sizeof ( rtDW . oatrzraqor ) ) ;
mr_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 , 34 , (
const void * ) & ( rtDW . h3yqjmxcau ) , sizeof ( rtDW . h3yqjmxcau ) ) ;
mr_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 , 35 , (
const void * ) & ( rtDW . c2gp1hcawr ) , sizeof ( rtDW . c2gp1hcawr ) ) ;
mr_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 , 36 , (
const void * ) & ( rtDW . kzezuc5kmj ) , sizeof ( rtDW . kzezuc5kmj ) ) ;
mr_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 , 37 , (
const void * ) & ( rtDW . baej0t4mz4 ) , sizeof ( rtDW . baej0t4mz4 ) ) ;
mr_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 , 38 , (
const void * ) & ( rtDW . e4o1f3himj ) , sizeof ( rtDW . e4o1f3himj ) ) ;
mr_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 , 39 , (
const void * ) & ( rtDW . etn0cabray ) , sizeof ( rtDW . etn0cabray ) ) ;
mr_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 , 40 , (
const void * ) & ( rtDW . nq5efgjr4u ) , sizeof ( rtDW . nq5efgjr4u ) ) ;
mr_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 , 41 , (
const void * ) & ( rtDW . ldwj0svuv5 ) , sizeof ( rtDW . ldwj0svuv5 ) ) ;
mr_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 , 42 , (
const void * ) & ( rtDW . et5jt4roea ) , sizeof ( rtDW . et5jt4roea ) ) ;
mr_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 , 43 , (
const void * ) & ( rtDW . jhgzfejyh1 ) , sizeof ( rtDW . jhgzfejyh1 ) ) ;
mr_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 , 44 , (
const void * ) & ( rtDW . e3b1amt0gc ) , sizeof ( rtDW . e3b1amt0gc ) ) ;
mr_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 , 45 , (
const void * ) & ( rtDW . g3swjv3kh4 ) , sizeof ( rtDW . g3swjv3kh4 ) ) ;
mr_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 , 46 , (
const void * ) & ( rtDW . cifi3jtean ) , sizeof ( rtDW . cifi3jtean ) ) ;
mr_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 , 47 , (
const void * ) & ( rtDW . mx1mzmmudm ) , sizeof ( rtDW . mx1mzmmudm ) ) ;
mr_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 , 48 , (
const void * ) & ( rtDW . nhmw0vweas ) , sizeof ( rtDW . nhmw0vweas ) ) ;
mr_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 , 49 , (
const void * ) & ( rtDW . enb1yfnvht ) , sizeof ( rtDW . enb1yfnvht ) ) ;
mr_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 , 50 , (
const void * ) & ( rtDW . au1n1jxucn ) , sizeof ( rtDW . au1n1jxucn ) ) ;
mr_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 , 51 , (
const void * ) & ( rtDW . cwu3pvv3pv ) , sizeof ( rtDW . cwu3pvv3pv ) ) ;
mr_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 , 52 , (
const void * ) & ( rtDW . dwzumrr4bn ) , sizeof ( rtDW . dwzumrr4bn ) ) ;
mr_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 , 53 , (
const void * ) & ( rtDW . jelszbraal ) , sizeof ( rtDW . jelszbraal ) ) ;
mr_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 , 54 , (
const void * ) & ( rtDW . cbdxrybiau ) , sizeof ( rtDW . cbdxrybiau ) ) ;
mr_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 , 55 , (
const void * ) & ( rtDW . dia2paiwqh ) , sizeof ( rtDW . dia2paiwqh ) ) ;
mr_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 , 56 , (
const void * ) & ( rtDW . f40yghkee4 ) , sizeof ( rtDW . f40yghkee4 ) ) ;
mr_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 , 57 , (
const void * ) & ( rtDW . pke0ukp1xy ) , sizeof ( rtDW . pke0ukp1xy ) ) ;
mr_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 , 58 , (
const void * ) & ( rtDW . ob1pklob1o ) , sizeof ( rtDW . ob1pklob1o ) ) ;
mr_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 , 59 , (
const void * ) & ( rtDW . jvf3pz3ask ) , sizeof ( rtDW . jvf3pz3ask ) ) ;
mr_trial_with_intermediate_tanks_cacheDataAsMxArray ( rtdwData , 0 , 60 , (
const void * ) & ( rtDW . cy4ormxj51 ) , sizeof ( rtDW . cy4ormxj51 ) ) ;
mxSetFieldByNumber ( ssDW , 0 , 1 , rtdwData ) ; }
mr_trial_with_intermediate_tanks_cacheDataAsMxArray ( ssDW , 0 , 2 , ( const
void * ) & ( rtPrevZCX ) , sizeof ( rtPrevZCX ) ) ; return ssDW ; } void
mr_trial_with_intermediate_tanks_SetDWork ( const mxArray * ssDW ) { ( void )
ssDW ; mr_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * ) &
( rtB ) , ssDW , 0 , 0 , sizeof ( rtB ) ) ; { const mxArray * rtdwData =
mxGetFieldByNumber ( ssDW , 0 , 1 ) ;
mr_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * ) & ( rtDW
. mog1ypd0na ) , rtdwData , 0 , 0 , sizeof ( rtDW . mog1ypd0na ) ) ;
mr_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * ) & ( rtDW
. d1ijg3op2k ) , rtdwData , 0 , 1 , sizeof ( rtDW . d1ijg3op2k ) ) ;
mr_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * ) & ( rtDW
. brka54cpul ) , rtdwData , 0 , 2 , sizeof ( rtDW . brka54cpul ) ) ;
mr_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * ) & ( rtDW
. ojcciswcd3 ) , rtdwData , 0 , 3 , sizeof ( rtDW . ojcciswcd3 ) ) ;
mr_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * ) & ( rtDW
. nfkviszx2j ) , rtdwData , 0 , 4 , sizeof ( rtDW . nfkviszx2j ) ) ;
mr_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * ) & ( rtDW
. ndyrr3tae0 ) , rtdwData , 0 , 5 , sizeof ( rtDW . ndyrr3tae0 ) ) ;
mr_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * ) & ( rtDW
. c2jtsctdjh ) , rtdwData , 0 , 6 , sizeof ( rtDW . c2jtsctdjh ) ) ;
mr_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * ) & ( rtDW
. mevs4smziy ) , rtdwData , 0 , 7 , sizeof ( rtDW . mevs4smziy ) ) ;
mr_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * ) & ( rtDW
. egiikq3opb ) , rtdwData , 0 , 8 , sizeof ( rtDW . egiikq3opb ) ) ;
mr_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * ) & ( rtDW
. puivctwyq5 ) , rtdwData , 0 , 9 , sizeof ( rtDW . puivctwyq5 ) ) ;
mr_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * ) & ( rtDW
. bx2j5ck3fy ) , rtdwData , 0 , 10 , sizeof ( rtDW . bx2j5ck3fy ) ) ;
mr_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * ) & ( rtDW
. fzispyfeqr ) , rtdwData , 0 , 11 , sizeof ( rtDW . fzispyfeqr ) ) ;
mr_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * ) & ( rtDW
. brdpyuefrq ) , rtdwData , 0 , 12 , sizeof ( rtDW . brdpyuefrq ) ) ;
mr_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * ) & ( rtDW
. esun0jx3xw ) , rtdwData , 0 , 13 , sizeof ( rtDW . esun0jx3xw ) ) ;
mr_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * ) & ( rtDW
. ncj5ml0ihc ) , rtdwData , 0 , 14 , sizeof ( rtDW . ncj5ml0ihc ) ) ;
mr_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * ) & ( rtDW
. kjn3oyfrac ) , rtdwData , 0 , 15 , sizeof ( rtDW . kjn3oyfrac ) ) ;
mr_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * ) & ( rtDW
. d3q0cyta13 ) , rtdwData , 0 , 16 , sizeof ( rtDW . d3q0cyta13 ) ) ;
mr_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * ) & ( rtDW
. di0xkhax5v ) , rtdwData , 0 , 17 , sizeof ( rtDW . di0xkhax5v ) ) ;
mr_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * ) & ( rtDW
. dol3avfpna ) , rtdwData , 0 , 18 , sizeof ( rtDW . dol3avfpna ) ) ;
mr_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * ) & ( rtDW
. mfw5n1ea1y ) , rtdwData , 0 , 19 , sizeof ( rtDW . mfw5n1ea1y ) ) ;
mr_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * ) & ( rtDW
. oj5q0jsluk ) , rtdwData , 0 , 20 , sizeof ( rtDW . oj5q0jsluk ) ) ;
mr_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * ) & ( rtDW
. dfi4io44l1 ) , rtdwData , 0 , 21 , sizeof ( rtDW . dfi4io44l1 ) ) ;
mr_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * ) & ( rtDW
. bapvsmhne1 ) , rtdwData , 0 , 22 , sizeof ( rtDW . bapvsmhne1 ) ) ;
mr_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * ) & ( rtDW
. gmzfgv4p33 ) , rtdwData , 0 , 23 , sizeof ( rtDW . gmzfgv4p33 ) ) ;
mr_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * ) & ( rtDW
. bxfqhgcloy ) , rtdwData , 0 , 24 , sizeof ( rtDW . bxfqhgcloy ) ) ;
mr_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * ) & ( rtDW
. gwwzsdaztz ) , rtdwData , 0 , 25 , sizeof ( rtDW . gwwzsdaztz ) ) ;
mr_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * ) & ( rtDW
. agicmdjqdi ) , rtdwData , 0 , 26 , sizeof ( rtDW . agicmdjqdi ) ) ;
mr_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * ) & ( rtDW
. astl0qyquk ) , rtdwData , 0 , 27 , sizeof ( rtDW . astl0qyquk ) ) ;
mr_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * ) & ( rtDW
. ekvnfwjtsb ) , rtdwData , 0 , 28 , sizeof ( rtDW . ekvnfwjtsb ) ) ;
mr_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * ) & ( rtDW
. dql53suw4t ) , rtdwData , 0 , 29 , sizeof ( rtDW . dql53suw4t ) ) ;
mr_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * ) & ( rtDW
. gz0frs1wku ) , rtdwData , 0 , 30 , sizeof ( rtDW . gz0frs1wku ) ) ;
mr_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * ) & ( rtDW
. iwrzpan53u ) , rtdwData , 0 , 31 , sizeof ( rtDW . iwrzpan53u ) ) ;
mr_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * ) & ( rtDW
. ivbuycouap ) , rtdwData , 0 , 32 , sizeof ( rtDW . ivbuycouap ) ) ;
mr_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * ) & ( rtDW
. oatrzraqor ) , rtdwData , 0 , 33 , sizeof ( rtDW . oatrzraqor ) ) ;
mr_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * ) & ( rtDW
. h3yqjmxcau ) , rtdwData , 0 , 34 , sizeof ( rtDW . h3yqjmxcau ) ) ;
mr_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * ) & ( rtDW
. c2gp1hcawr ) , rtdwData , 0 , 35 , sizeof ( rtDW . c2gp1hcawr ) ) ;
mr_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * ) & ( rtDW
. kzezuc5kmj ) , rtdwData , 0 , 36 , sizeof ( rtDW . kzezuc5kmj ) ) ;
mr_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * ) & ( rtDW
. baej0t4mz4 ) , rtdwData , 0 , 37 , sizeof ( rtDW . baej0t4mz4 ) ) ;
mr_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * ) & ( rtDW
. e4o1f3himj ) , rtdwData , 0 , 38 , sizeof ( rtDW . e4o1f3himj ) ) ;
mr_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * ) & ( rtDW
. etn0cabray ) , rtdwData , 0 , 39 , sizeof ( rtDW . etn0cabray ) ) ;
mr_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * ) & ( rtDW
. nq5efgjr4u ) , rtdwData , 0 , 40 , sizeof ( rtDW . nq5efgjr4u ) ) ;
mr_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * ) & ( rtDW
. ldwj0svuv5 ) , rtdwData , 0 , 41 , sizeof ( rtDW . ldwj0svuv5 ) ) ;
mr_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * ) & ( rtDW
. et5jt4roea ) , rtdwData , 0 , 42 , sizeof ( rtDW . et5jt4roea ) ) ;
mr_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * ) & ( rtDW
. jhgzfejyh1 ) , rtdwData , 0 , 43 , sizeof ( rtDW . jhgzfejyh1 ) ) ;
mr_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * ) & ( rtDW
. e3b1amt0gc ) , rtdwData , 0 , 44 , sizeof ( rtDW . e3b1amt0gc ) ) ;
mr_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * ) & ( rtDW
. g3swjv3kh4 ) , rtdwData , 0 , 45 , sizeof ( rtDW . g3swjv3kh4 ) ) ;
mr_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * ) & ( rtDW
. cifi3jtean ) , rtdwData , 0 , 46 , sizeof ( rtDW . cifi3jtean ) ) ;
mr_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * ) & ( rtDW
. mx1mzmmudm ) , rtdwData , 0 , 47 , sizeof ( rtDW . mx1mzmmudm ) ) ;
mr_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * ) & ( rtDW
. nhmw0vweas ) , rtdwData , 0 , 48 , sizeof ( rtDW . nhmw0vweas ) ) ;
mr_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * ) & ( rtDW
. enb1yfnvht ) , rtdwData , 0 , 49 , sizeof ( rtDW . enb1yfnvht ) ) ;
mr_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * ) & ( rtDW
. au1n1jxucn ) , rtdwData , 0 , 50 , sizeof ( rtDW . au1n1jxucn ) ) ;
mr_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * ) & ( rtDW
. cwu3pvv3pv ) , rtdwData , 0 , 51 , sizeof ( rtDW . cwu3pvv3pv ) ) ;
mr_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * ) & ( rtDW
. dwzumrr4bn ) , rtdwData , 0 , 52 , sizeof ( rtDW . dwzumrr4bn ) ) ;
mr_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * ) & ( rtDW
. jelszbraal ) , rtdwData , 0 , 53 , sizeof ( rtDW . jelszbraal ) ) ;
mr_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * ) & ( rtDW
. cbdxrybiau ) , rtdwData , 0 , 54 , sizeof ( rtDW . cbdxrybiau ) ) ;
mr_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * ) & ( rtDW
. dia2paiwqh ) , rtdwData , 0 , 55 , sizeof ( rtDW . dia2paiwqh ) ) ;
mr_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * ) & ( rtDW
. f40yghkee4 ) , rtdwData , 0 , 56 , sizeof ( rtDW . f40yghkee4 ) ) ;
mr_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * ) & ( rtDW
. pke0ukp1xy ) , rtdwData , 0 , 57 , sizeof ( rtDW . pke0ukp1xy ) ) ;
mr_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * ) & ( rtDW
. ob1pklob1o ) , rtdwData , 0 , 58 , sizeof ( rtDW . ob1pklob1o ) ) ;
mr_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * ) & ( rtDW
. jvf3pz3ask ) , rtdwData , 0 , 59 , sizeof ( rtDW . jvf3pz3ask ) ) ;
mr_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * ) & ( rtDW
. cy4ormxj51 ) , rtdwData , 0 , 60 , sizeof ( rtDW . cy4ormxj51 ) ) ; }
mr_trial_with_intermediate_tanks_restoreDataFromMxArray ( ( void * ) & (
rtPrevZCX ) , ssDW , 0 , 2 , sizeof ( rtPrevZCX ) ) ; } mxArray *
mr_trial_with_intermediate_tanks_GetSimStateDisallowedBlocks ( ) { mxArray *
data = mxCreateCellMatrix ( 8 , 3 ) ; mwIndex subs [ 2 ] , offset ; { static
const char_T * blockType [ 8 ] = { "Scope" , "Scope" , "Scope" , "Scope" ,
"Scope" , "Scope" , "Scope" , "Scope" , } ; static const char_T * blockPath [
8 ] = { "trial_with_intermediate_tanks/Energy Comparison" ,
"trial_with_intermediate_tanks/Scope" ,
"trial_with_intermediate_tanks/Scope1" ,
"trial_with_intermediate_tanks/Scope2" ,
"trial_with_intermediate_tanks/Scope3" ,
"trial_with_intermediate_tanks/Scope4" ,
"trial_with_intermediate_tanks/Scope5" ,
"trial_with_intermediate_tanks/Scope6" , } ; static const int reason [ 8 ] =
{ 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , } ; for ( subs [ 0 ] = 0 ; subs [ 0 ] < 8 ;
++ ( subs [ 0 ] ) ) { subs [ 1 ] = 0 ; offset = mxCalcSingleSubscript ( data
, 2 , subs ) ; mxSetCell ( data , offset , mxCreateString ( blockType [ subs
[ 0 ] ] ) ) ; subs [ 1 ] = 1 ; offset = mxCalcSingleSubscript ( data , 2 ,
subs ) ; mxSetCell ( data , offset , mxCreateString ( blockPath [ subs [ 0 ]
] ) ) ; subs [ 1 ] = 2 ; offset = mxCalcSingleSubscript ( data , 2 , subs ) ;
mxSetCell ( data , offset , mxCreateDoubleScalar ( ( real_T ) reason [ subs [
0 ] ] ) ) ; } } return data ; } void MdlInitializeSizes ( void ) {
ssSetNumContStates ( rtS , 27 ) ; ssSetNumPeriodicContStates ( rtS , 0 ) ;
ssSetNumY ( rtS , 0 ) ; ssSetNumU ( rtS , 0 ) ; ssSetDirectFeedThrough ( rtS
, 0 ) ; ssSetNumSampleTimes ( rtS , 3 ) ; ssSetNumBlocks ( rtS , 146 ) ;
ssSetNumBlockIO ( rtS , 79 ) ; ssSetNumBlockParams ( rtS , 56 ) ; } void
MdlInitializeSampleTimes ( void ) { ssSetSampleTime ( rtS , 0 , 0.0 ) ;
ssSetSampleTime ( rtS , 1 , 0.0 ) ; ssSetSampleTime ( rtS , 2 , 1.0E-6 ) ;
ssSetOffsetTime ( rtS , 0 , 0.0 ) ; ssSetOffsetTime ( rtS , 1 , 1.0 ) ;
ssSetOffsetTime ( rtS , 2 , 0.0 ) ; } void raccel_set_checksum ( ) {
ssSetChecksumVal ( rtS , 0 , 2806597705U ) ; ssSetChecksumVal ( rtS , 1 ,
165651654U ) ; ssSetChecksumVal ( rtS , 2 , 2791374129U ) ; ssSetChecksumVal
( rtS , 3 , 1143219654U ) ; }
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
trial_with_intermediate_tanks_InitializeDataMapInfo ( ) ;
ssSetIsRapidAcceleratorActive ( rtS , true ) ; ssSetRootSS ( rtS , rtS ) ;
ssSetVersion ( rtS , SIMSTRUCT_VERSION_LEVEL2 ) ; ssSetModelName ( rtS ,
"trial_with_intermediate_tanks" ) ; ssSetPath ( rtS ,
"trial_with_intermediate_tanks" ) ; ssSetTStart ( rtS , 0.0 ) ; ssSetTFinal (
rtS , 0.16 ) ; { static RTWLogInfo rt_DataLoggingInfo ; rt_DataLoggingInfo .
loggingInterval = ( NULL ) ; ssSetRTWLogInfo ( rtS , & rt_DataLoggingInfo ) ;
} { { static int_T rt_LoggedStateWidths [ ] = { 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1
, 1 , 1 , 1 , 1 , 1 , 1 , 1 , 2 , 2 , 2 , 2 , 1 , 1 , 2 , 1 } ; static int_T
rt_LoggedStateNumDimensions [ ] = { 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1
, 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 } ; static int_T
rt_LoggedStateDimensions [ ] = { 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 ,
1 , 1 , 1 , 1 , 2 , 2 , 2 , 2 , 1 , 1 , 2 , 1 } ; static boolean_T
rt_LoggedStateIsVarDims [ ] = { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0
, 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 } ; static BuiltInDTypeId
rt_LoggedStateDataTypeIds [ ] = { SS_DOUBLE , SS_DOUBLE , SS_DOUBLE ,
SS_DOUBLE , SS_DOUBLE , SS_DOUBLE , SS_DOUBLE , SS_DOUBLE , SS_DOUBLE ,
SS_DOUBLE , SS_DOUBLE , SS_DOUBLE , SS_DOUBLE , SS_DOUBLE , SS_DOUBLE ,
SS_DOUBLE , SS_DOUBLE , SS_DOUBLE , SS_DOUBLE , SS_DOUBLE , SS_DOUBLE ,
SS_DOUBLE , SS_DOUBLE } ; static int_T rt_LoggedStateComplexSignals [ ] = { 0
, 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ,
0 , 0 , 0 } ; static RTWPreprocessingFcnPtr
rt_LoggingStatePreprocessingFcnPtrs [ ] = { ( NULL ) , ( NULL ) , ( NULL ) ,
( NULL ) , ( NULL ) , ( NULL ) , ( NULL ) , ( NULL ) , ( NULL ) , ( NULL ) ,
( NULL ) , ( NULL ) , ( NULL ) , ( NULL ) , ( NULL ) , ( NULL ) , ( NULL ) ,
( NULL ) , ( NULL ) , ( NULL ) , ( NULL ) , ( NULL ) , ( NULL ) } ; static
const char_T * rt_LoggedStateLabels [ ] = { "CSTATE" , "CSTATE" , "CSTATE" ,
"CSTATE" , "CSTATE" , "CSTATE" , "CSTATE" , "CSTATE" , "CSTATE" , "CSTATE" ,
"CSTATE" , "CSTATE" , "CSTATE" , "CSTATE" , "CSTATE" , "CSTATE" , "CSTATE" ,
"CSTATE" , "CSTATE" , "CSTATE" , "CSTATE" , "CSTATE" , "DSTATE" } ; static
const char_T * rt_LoggedStateBlockNames [ ] = {
"trial_with_intermediate_tanks/Integrator4" ,
"trial_with_intermediate_tanks/Integrator2" ,
"trial_with_intermediate_tanks/Integrator13" ,
"trial_with_intermediate_tanks/Integrator" ,
"trial_with_intermediate_tanks/Integrator1" ,
"trial_with_intermediate_tanks/Integrator14" ,
"trial_with_intermediate_tanks/Integrator7" ,
"trial_with_intermediate_tanks/Integrator12" ,
"trial_with_intermediate_tanks/Integrator10" ,
"trial_with_intermediate_tanks/Integrator8" ,
"trial_with_intermediate_tanks/Integrator9" ,
"trial_with_intermediate_tanks/Integrator11" ,
"trial_with_intermediate_tanks/Integrator3" ,
"trial_with_intermediate_tanks/Integrator6" ,
"trial_with_intermediate_tanks/Integrator16" ,
"trial_with_intermediate_tanks/Dynamics of Valve1/Transfer Fcn" ,
"trial_with_intermediate_tanks/Dynamics of Valve2/Transfer Fcn" ,
"trial_with_intermediate_tanks/Dynamics of Valve3/Transfer Fcn" ,
"trial_with_intermediate_tanks/Dynamics of Valve4/Transfer Fcn" ,
"trial_with_intermediate_tanks/Integrator5" ,
"trial_with_intermediate_tanks/Integrator15" ,
"trial_with_intermediate_tanks/Dynamics of Valve/Transfer Fcn" ,
"trial_with_intermediate_tanks/Delay" } ; static const char_T *
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
) ; rt_LoggedStateSignalPtrs [ 0 ] = ( void * ) & rtX . njnhyztmul ;
rt_LoggedStateSignalPtrs [ 1 ] = ( void * ) & rtX . cxtidlrvgx ;
rt_LoggedStateSignalPtrs [ 2 ] = ( void * ) & rtX . anhq5q4adb ;
rt_LoggedStateSignalPtrs [ 3 ] = ( void * ) & rtX . h1nzd5kpbz ;
rt_LoggedStateSignalPtrs [ 4 ] = ( void * ) & rtX . pmj0r0j4v2 ;
rt_LoggedStateSignalPtrs [ 5 ] = ( void * ) & rtX . dusu5owudw ;
rt_LoggedStateSignalPtrs [ 6 ] = ( void * ) & rtX . decqh12jsd ;
rt_LoggedStateSignalPtrs [ 7 ] = ( void * ) & rtX . f3sxyupuw4 ;
rt_LoggedStateSignalPtrs [ 8 ] = ( void * ) & rtX . javki5lnyj ;
rt_LoggedStateSignalPtrs [ 9 ] = ( void * ) & rtX . d1tt3gwruv ;
rt_LoggedStateSignalPtrs [ 10 ] = ( void * ) & rtX . erk0rjtwtw ;
rt_LoggedStateSignalPtrs [ 11 ] = ( void * ) & rtX . m2od42nnln ;
rt_LoggedStateSignalPtrs [ 12 ] = ( void * ) & rtX . ith00xbmqs ;
rt_LoggedStateSignalPtrs [ 13 ] = ( void * ) & rtX . djreaptdva ;
rt_LoggedStateSignalPtrs [ 14 ] = ( void * ) & rtX . pnlrwk52wo ;
rt_LoggedStateSignalPtrs [ 15 ] = ( void * ) & rtX . lxsclehbdn [ 0 ] ;
rt_LoggedStateSignalPtrs [ 16 ] = ( void * ) & rtX . herd5wji3e [ 0 ] ;
rt_LoggedStateSignalPtrs [ 17 ] = ( void * ) & rtX . ekhpkuwiae [ 0 ] ;
rt_LoggedStateSignalPtrs [ 18 ] = ( void * ) & rtX . k2eyyio13c [ 0 ] ;
rt_LoggedStateSignalPtrs [ 19 ] = ( void * ) & rtX . mh2ntyj1vq ;
rt_LoggedStateSignalPtrs [ 20 ] = ( void * ) & rtX . kh4uwqwwas ;
rt_LoggedStateSignalPtrs [ 21 ] = ( void * ) & rtX . l1m0lcksjp [ 0 ] ;
rt_LoggedStateSignalPtrs [ 22 ] = ( void * ) & rtDW . mog1ypd0na ; }
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
[ 6 ] = { { 1 * sizeof ( real_T ) , ( char * ) ( & rtB . jsgd1kggkx ) , (
NULL ) } , { 1 * sizeof ( real_T ) , ( char * ) ( & rtB . a5ddnr2k0h ) , (
NULL ) } , { 1 * sizeof ( real_T ) , ( char * ) ( & rtB . bnvm1inada ) , (
NULL ) } , { 1 * sizeof ( real_T ) , ( char * ) ( & rtB . bojkujnbsq ) , (
NULL ) } , { 1 * sizeof ( real_T ) , ( char * ) ( & rtB . jbbgrda5ov ) , (
NULL ) } , { 1 * sizeof ( real_T ) , ( char * ) ( & rtB . jjwnm1umsr ) , (
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
zc ) ; } { rtPrevZCX . ii2vsk3a4a = UNINITIALIZED_ZCSIG ; } ssSetChecksumVal
( rtS , 0 , 2806597705U ) ; ssSetChecksumVal ( rtS , 1 , 165651654U ) ;
ssSetChecksumVal ( rtS , 2 , 2791374129U ) ; ssSetChecksumVal ( rtS , 3 ,
1143219654U ) ; { static const sysRanDType rtAlwaysEnabled =
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
rtDW . astl0qyquk ; systemRan [ 14 ] = & rtAlwaysEnabled ; systemRan [ 15 ] =
& rtAlwaysEnabled ; systemRan [ 16 ] = & rtAlwaysEnabled ; systemRan [ 17 ] =
& rtAlwaysEnabled ; rteiSetModelMappingInfoPtr ( ssGetRTWExtModeInfo ( rtS )
, & ssGetModelMappingInfo ( rtS ) ) ; rteiSetChecksumsPtr (
ssGetRTWExtModeInfo ( rtS ) , ssGetChecksums ( rtS ) ) ; rteiSetTPtr (
ssGetRTWExtModeInfo ( rtS ) , ssGetTPtr ( rtS ) ) ; }
slsaDisallowedBlocksForSimTargetOP ( rtS ,
mr_trial_with_intermediate_tanks_GetSimStateDisallowedBlocks ) ;
slsaGetWorkFcnForSimTargetOP ( rtS ,
mr_trial_with_intermediate_tanks_GetDWork ) ; slsaSetWorkFcnForSimTargetOP (
rtS , mr_trial_with_intermediate_tanks_SetDWork ) ;
rt_RapidReadMatFileAndUpdateParams ( rtS ) ; if ( ssGetErrorStatus ( rtS ) )
{ return rtS ; } return rtS ; }
#if defined(_MSC_VER)
#pragma optimize( "", on )
#endif
void MdlOutputsParameterSampleTime ( int_T tid ) { MdlOutputsTID3 ( tid ) ; }
