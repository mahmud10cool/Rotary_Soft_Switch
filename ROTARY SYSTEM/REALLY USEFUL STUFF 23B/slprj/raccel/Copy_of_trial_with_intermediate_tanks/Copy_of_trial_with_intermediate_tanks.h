#ifndef RTW_HEADER_Copy_of_trial_with_intermediate_tanks_h_
#define RTW_HEADER_Copy_of_trial_with_intermediate_tanks_h_
#ifndef Copy_of_trial_with_intermediate_tanks_COMMON_INCLUDES_
#define Copy_of_trial_with_intermediate_tanks_COMMON_INCLUDES_
#include <stdlib.h>
#include "sl_AsyncioQueue/AsyncioQueueCAPI.h"
#include "rtwtypes.h"
#include "sigstream_rtw.h"
#include "simtarget/slSimTgtSigstreamRTW.h"
#include "simtarget/slSimTgtSlioCoreRTW.h"
#include "simtarget/slSimTgtSlioClientsRTW.h"
#include "simtarget/slSimTgtSlioSdiRTW.h"
#include "simstruc.h"
#include "fixedpoint.h"
#include "raccel.h"
#include "slsv_diagnostic_codegen_c_api.h"
#include "rt_logging_simtarget.h"
#include "dt_info.h"
#include "ext_work.h"
#endif
#include "Copy_of_trial_with_intermediate_tanks_types.h"
#include "rt_zcfcn.h"
#include <stddef.h>
#include "rtw_modelmap_simtarget.h"
#include "rt_defines.h"
#include <string.h>
#include "rtGetInf.h"
#include "rt_nonfinite.h"
#include "zero_crossing_types.h"
#define MODEL_NAME Copy_of_trial_with_intermediate_tanks
#define NSAMPLE_TIMES (4) 
#define NINPUTS (0)       
#define NOUTPUTS (0)     
#define NBLOCKIO (79) 
#define NUM_ZC_EVENTS (1) 
#ifndef NCSTATES
#define NCSTATES (27)   
#elif NCSTATES != 27
#error Invalid specification of NCSTATES defined in compiler command
#endif
#ifndef rtmGetDataMapInfo
#define rtmGetDataMapInfo(rtm) (*rt_dataMapInfoPtr)
#endif
#ifndef rtmSetDataMapInfo
#define rtmSetDataMapInfo(rtm, val) (rt_dataMapInfoPtr = &val)
#endif
#ifndef IN_RACCEL_MAIN
#endif
typedef struct { real_T mzfn3lhia4 ; real_T jswgtkntlr ; real_T pyky0iaydb ;
real_T cwl31aims5 ; real_T dk50fnzsks ; real_T lr3evndp05 ; real_T gaymsnblvn
; real_T hfumqnwosi ; real_T jcnaw1aobb ; real_T c5jlsr0y2m ; real_T
gjajnpeavu ; real_T onszijeutj ; real_T o5jgaxymjd ; real_T aqvwe35511 ;
real_T ga4spwnzmb ; real_T iqcx4ji2gy ; real_T b5jzfqhswx ; real_T o5psjsgjbl
; real_T plcvs3tnfq ; real_T j0krl2hmjm ; real_T l0b3kthfn5 ; real_T
g5clq5nhx1 ; real_T mfk254enll ; real_T id3dy3zvhz ; real_T exiuxpqhzf ;
real_T jjqotddbo4 ; real_T hwx2zv05iw ; real_T n3lttcqcta ; real_T klsduu5tgr
; real_T huebprkurf ; real_T limtet2woe ; real_T ap3x1johma ; real_T
iiwlxh5nit ; real_T eiqmkhbtyt ; real_T btfjf0eyh0 ; real_T cbizxpnxe4 ;
real_T gzd0bvtkjr ; real_T m3tgvpb0je ; real_T kfmldcnjnq ; real_T hdgz1v45gg
; real_T lzxak0a220 ; real_T faofsw34yz ; real_T oz3rry1khh ; real_T
jhxf0vx510 ; real_T ahjuj0uqhz ; real_T eidzbo2axi ; real_T dzr21vacs5 ;
real_T g5dn334cgf ; real_T hllpacqnya ; real_T fad5yhgoij ; real_T nzchyezcdl
; real_T ouvddzwv14 ; real_T ixdckelvxz ; real_T edzhheuutq ; real_T
ci3hl3yrta ; real_T fgm5y4v5wu ; real_T ciwvymzmbd ; real_T bb2zjb32zc ;
real_T b0rncimetq ; real_T eks3p0husl ; real_T miudztndok ; real_T aufytnaglh
; real_T oguodtbfqn ; real_T djfhnx1214 ; real_T afrqwefb1y ; real_T
n2srwnob2z ; real_T folkmnz0t4 ; real_T k35lmpwnot ; real_T gnjhdq1qvh ;
real_T hmcqeqif2o ; real_T owtjgmweb5 ; real_T jpwjwhvvnr ; real_T pju2ibjy53
; real_T m0rkotwez0 ; real_T fsfzk25tf0 ; real_T jkrmpymttx ; real_T
mrxnnzezb3 ; real_T eqrybinqmz ; real_T hx5qshmezi ; } B ; typedef struct {
real_T j0cswhuhs0 ; struct { void * LoggedData [ 2 ] ; } mg3u3mjjuf ; struct
{ void * LoggedData ; } iljidhxeal ; struct { void * LoggedData [ 4 ] ; }
gfpi0a3qtl ; struct { void * LoggedData ; } p3eoc0sotd ; struct { void *
LoggedData [ 3 ] ; } naarlwo2lz ; struct { void * LoggedData ; } km0sbvoafb ;
struct { void * LoggedData ; } lrvek4v1dp ; struct { void * AQHandles ; }
homjuevopz ; struct { void * AQHandles ; } h1nrd1cc1r ; struct { void *
AQHandles ; } l444m1qlra ; struct { void * AQHandles ; } fcjl52ojvi ; struct
{ void * AQHandles ; } cpiqtbmrnh ; struct { void * AQHandles ; } bya1bx33so
; struct { void * AQHandles ; } h1qcylur1s ; struct { void * AQHandles ; }
ohqxmkyead ; struct { void * AQHandles ; } jnhubbfo5j ; struct { void *
AQHandles ; } fokz15rol1 ; struct { void * AQHandles ; } lj0dkmhahi ; struct
{ void * AQHandles ; } oqdb4h0bik ; struct { void * LoggedData ; } gkxaf1jyuy
; int32_T atxpvhc5n0 ; int32_T inlqdh4kh1 ; int32_T glwrm0q1go ; int32_T
apkag20puf ; int32_T mpaqmqlwie ; int32_T a4pqll03x1 ; int32_T p1ctkshy11 ;
int32_T d34s0zbduf ; int32_T epcaebihmu ; int32_T kvbq5jwbvw ; int32_T
audz3xgmq3 ; int32_T loaiy2cacx ; int32_T e0qqpzfruh ; int32_T hamy5qs5qe ;
int32_T je5l4cyavv ; int32_T ebwtwofxpd ; int_T aku4xs4fdc ; int_T aigzt14viy
; int_T jwuiy2atz3 ; int_T isc2qgvxqr ; int_T awwuosv0rp ; int_T g40nougyde ;
int_T cangmlit2r ; int_T od3s0wtt1r ; int_T ppnlaclhj1 ; int_T exujqdty5l ;
int8_T iln55qrgxk ; uint8_T i213auut4h ; uint8_T e13oplclke ; uint8_T
laquyfl1kp ; uint8_T f45mx2s03m ; uint8_T gstdmwepbc ; uint8_T ikvzyivi0j ;
uint8_T libi2dsvky ; uint8_T ga2xty0o4o ; uint8_T e1b1tx1d4l ; uint8_T
p0eq4zkt3w ; uint8_T oahmxad2ft ; uint8_T ksvy5jh2vc ; uint8_T fyaulsukmr ;
uint8_T hdus1cnzep ; uint8_T amiin4kdmj ; uint8_T nwxy1ani3r ; boolean_T
j5zzepflri ; boolean_T hk25kkxqwm ; boolean_T ffd0olavlx ; boolean_T
eiohhcwuqa ; boolean_T eyyp3zs3le ; boolean_T nakeaw2rna ; boolean_T
n4ne2vzdw3 ; boolean_T ljfwvqlys1 ; boolean_T ct4cduudvv ; boolean_T
jkvzly4a1o ; boolean_T fxawe3htlm ; boolean_T irq32iym0t ; boolean_T
jz3v2t25p2 ; boolean_T m0i52xetkk ; boolean_T dfvcjzwy4h ; boolean_T
nas3vtbw5t ; boolean_T awqg33it32 ; } DW ; typedef struct { real_T eer5ulgpq1
; real_T gvqywnxquo ; real_T mui3hk2kdn ; real_T kllv1nez4n ; real_T
orxgnzk15p ; real_T apttuwn22n ; real_T avebadcay4 ; real_T o2kfetkxay ;
real_T fvfvpkek2n ; real_T aqqzk12fdv ; real_T e30fqjkckb ; real_T mizgi5ztop
; real_T iikwtjyqxf ; real_T c5bbd2cwaf ; real_T loqzkyynb5 ; real_T
klye04tpmo [ 2 ] ; real_T ht5behjiqe [ 2 ] ; real_T aa2dkshjtd [ 2 ] ; real_T
kco1vq4cxg [ 2 ] ; real_T fklbo31ygf ; real_T kwoozlzy2w ; real_T mvzjw32hje
[ 2 ] ; } X ; typedef struct { real_T eer5ulgpq1 ; real_T gvqywnxquo ; real_T
mui3hk2kdn ; real_T kllv1nez4n ; real_T orxgnzk15p ; real_T apttuwn22n ;
real_T avebadcay4 ; real_T o2kfetkxay ; real_T fvfvpkek2n ; real_T aqqzk12fdv
; real_T e30fqjkckb ; real_T mizgi5ztop ; real_T iikwtjyqxf ; real_T
c5bbd2cwaf ; real_T loqzkyynb5 ; real_T klye04tpmo [ 2 ] ; real_T ht5behjiqe
[ 2 ] ; real_T aa2dkshjtd [ 2 ] ; real_T kco1vq4cxg [ 2 ] ; real_T fklbo31ygf
; real_T kwoozlzy2w ; real_T mvzjw32hje [ 2 ] ; } XDot ; typedef struct {
boolean_T eer5ulgpq1 ; boolean_T gvqywnxquo ; boolean_T mui3hk2kdn ;
boolean_T kllv1nez4n ; boolean_T orxgnzk15p ; boolean_T apttuwn22n ;
boolean_T avebadcay4 ; boolean_T o2kfetkxay ; boolean_T fvfvpkek2n ;
boolean_T aqqzk12fdv ; boolean_T e30fqjkckb ; boolean_T mizgi5ztop ;
boolean_T iikwtjyqxf ; boolean_T c5bbd2cwaf ; boolean_T loqzkyynb5 ;
boolean_T klye04tpmo [ 2 ] ; boolean_T ht5behjiqe [ 2 ] ; boolean_T
aa2dkshjtd [ 2 ] ; boolean_T kco1vq4cxg [ 2 ] ; boolean_T fklbo31ygf ;
boolean_T kwoozlzy2w ; boolean_T mvzjw32hje [ 2 ] ; } XDis ; typedef struct {
real_T eer5ulgpq1 ; real_T gvqywnxquo ; real_T mui3hk2kdn ; real_T kllv1nez4n
; real_T orxgnzk15p ; real_T apttuwn22n ; real_T avebadcay4 ; real_T
o2kfetkxay ; real_T fvfvpkek2n ; real_T aqqzk12fdv ; real_T e30fqjkckb ;
real_T mizgi5ztop ; real_T iikwtjyqxf ; real_T c5bbd2cwaf ; real_T loqzkyynb5
; real_T klye04tpmo [ 2 ] ; real_T ht5behjiqe [ 2 ] ; real_T aa2dkshjtd [ 2 ]
; real_T kco1vq4cxg [ 2 ] ; real_T fklbo31ygf ; real_T kwoozlzy2w ; real_T
mvzjw32hje [ 2 ] ; } CStateAbsTol ; typedef struct { real_T eer5ulgpq1 ;
real_T gvqywnxquo ; real_T mui3hk2kdn ; real_T kllv1nez4n ; real_T orxgnzk15p
; real_T apttuwn22n ; real_T avebadcay4 ; real_T o2kfetkxay ; real_T
fvfvpkek2n ; real_T aqqzk12fdv ; real_T e30fqjkckb ; real_T mizgi5ztop ;
real_T iikwtjyqxf ; real_T c5bbd2cwaf ; real_T loqzkyynb5 ; real_T klye04tpmo
[ 2 ] ; real_T ht5behjiqe [ 2 ] ; real_T aa2dkshjtd [ 2 ] ; real_T kco1vq4cxg
[ 2 ] ; real_T fklbo31ygf ; real_T kwoozlzy2w ; real_T mvzjw32hje [ 2 ] ; }
CXPtMin ; typedef struct { real_T eer5ulgpq1 ; real_T gvqywnxquo ; real_T
mui3hk2kdn ; real_T kllv1nez4n ; real_T orxgnzk15p ; real_T apttuwn22n ;
real_T avebadcay4 ; real_T o2kfetkxay ; real_T fvfvpkek2n ; real_T aqqzk12fdv
; real_T e30fqjkckb ; real_T mizgi5ztop ; real_T iikwtjyqxf ; real_T
c5bbd2cwaf ; real_T loqzkyynb5 ; real_T klye04tpmo [ 2 ] ; real_T ht5behjiqe
[ 2 ] ; real_T aa2dkshjtd [ 2 ] ; real_T kco1vq4cxg [ 2 ] ; real_T fklbo31ygf
; real_T kwoozlzy2w ; real_T mvzjw32hje [ 2 ] ; } CXPtMax ; typedef struct {
real_T krptktw1n2 ; real_T py1rahqile ; real_T a0e1ektqfm ; real_T lonlr0pekt
; real_T ghimdihoai ; real_T hen3ks1v4p ; real_T kij5fz4qdb ; real_T
l4cui55vk5 ; real_T h2qdorbmuq ; real_T pseezcu2gs ; real_T i3mlogp5q3 ;
real_T azuoncxh5k ; real_T nmdinmuvcs ; } ZCV ; typedef struct { ZCSigState
ah03fydvt5 ; } PrevZCX ; typedef struct { rtwCAPI_ModelMappingInfo mmi ; }
DataMapInfo ; struct P_ { struct_BLEksg6c1ggfsVfFA9HVZ param ; real_T Out1_Y0
; real_T Integrator4_IC ; real_T Gain_Gain ; real_T Relay_OnVal ; real_T
Relay_OffVal ; real_T Relay_YOn ; real_T Relay_YOff ; real_T Integrator13_IC
; real_T Integrator14_IC ; real_T Integrator7_IC ; real_T Integrator12_IC ;
real_T Integrator10_IC ; real_T Integrator8_IC ; real_T Integrator9_IC ;
real_T Integrator11_IC ; real_T Gain3_Gain ; real_T Integrator6_IC ; real_T
Integrator16_IC ; real_T TransferFcn_A [ 2 ] ; real_T TransferFcn_C [ 2 ] ;
real_T Saturation_LowerSat ; real_T TransferFcn_A_msxihde20w [ 2 ] ; real_T
TransferFcn_C_nnurr522qa [ 2 ] ; real_T Saturation_LowerSat_l0ksmwxpmq ;
real_T TransferFcn_A_bwuzoqq3iu [ 2 ] ; real_T TransferFcn_C_j01bq5y01n [ 2 ]
; real_T Saturation_LowerSat_m0gx2ayxie ; real_T TransferFcn_A_cqt2tjfue3 [ 2
] ; real_T TransferFcn_C_kavhf31tp1 [ 2 ] ; real_T
Saturation_LowerSat_l5sm1vdo10 ; real_T Integrator5_IC ; real_T
Delay_InitialCondition ; real_T Velocity_Amp ; real_T Velocity_Bias ; real_T
Velocity_Freq ; real_T Velocity_Phase ; real_T Integrator15_IC ; real_T
TransferFcn_A_m2qs34lzon [ 2 ] ; real_T TransferFcn_C_lu1wzht1on [ 2 ] ;
real_T Saturation_LowerSat_mclkng5niv ; real_T InitialSpeed_Value ; real_T
Gain5_Gain ; real_T Constant_Value ; real_T Constant5_Value ; uint8_T
ManualSwitch_CurrentSetting ; } ; extern const char_T *
RT_MEMORY_ALLOCATION_ERROR ; extern B rtB ; extern X rtX ; extern DW rtDW ;
extern PrevZCX rtPrevZCX ; extern P rtP ; extern mxArray *
mr_Copy_of_trial_with_intermediate_tanks_GetDWork ( ) ; extern void
mr_Copy_of_trial_with_intermediate_tanks_SetDWork ( const mxArray * ssDW ) ;
extern mxArray *
mr_Copy_of_trial_with_intermediate_tanks_GetSimStateDisallowedBlocks ( ) ;
extern const rtwCAPI_ModelMappingStaticInfo *
Copy_of_trial_with_intermediate_tanks_GetCAPIStaticMap ( void ) ; extern
SimStruct * const rtS ; extern DataMapInfo * rt_dataMapInfoPtr ; extern
rtwCAPI_ModelMappingInfo * rt_modelMapInfoPtr ; void MdlOutputs ( int_T tid )
; void MdlOutputsParameterSampleTime ( int_T tid ) ; void MdlUpdate ( int_T
tid ) ; void MdlTerminate ( void ) ; void MdlInitializeSizes ( void ) ; void
MdlInitializeSampleTimes ( void ) ; SimStruct * raccel_register_model (
ssExecutionInfo * executionInfo ) ;
#endif
