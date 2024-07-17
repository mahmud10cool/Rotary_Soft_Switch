#ifndef RTW_HEADER_trial_with_intermediate_tanks_h_
#define RTW_HEADER_trial_with_intermediate_tanks_h_
#ifndef trial_with_intermediate_tanks_COMMON_INCLUDES_
#define trial_with_intermediate_tanks_COMMON_INCLUDES_
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
#include "trial_with_intermediate_tanks_types.h"
#include "rt_zcfcn.h"
#include <stddef.h>
#include "rtw_modelmap_simtarget.h"
#include "rt_defines.h"
#include <string.h>
#include "rtGetInf.h"
#include "rt_nonfinite.h"
#include "zero_crossing_types.h"
#define MODEL_NAME trial_with_intermediate_tanks
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
typedef struct { real_T du3hz2fyzs ; real_T frh4wmn3vf ; real_T mcwb4pv1pd ;
real_T cqrrgwzxyv ; real_T dxeimajz1j ; real_T a5ddnr2k0h ; real_T blk35wfutl
; real_T ikogzwrhkz ; real_T fakrdhsmg5 ; real_T cy5jq00nod ; real_T
ht3vr3fcw1 ; real_T axh1siq1cq ; real_T fkszch0qwe ; real_T idqoyikyv2 ;
real_T hgyc45hfyg ; real_T ku0vq41nvc ; real_T cjeggiimao ; real_T ab5hqtcfdy
; real_T ip00peawg3 ; real_T iubu1mncgz ; real_T na1kfpe1mt ; real_T
eee41maadj ; real_T pebffqxdnd ; real_T p2rr4uzpkw ; real_T fj3bxcs4ue ;
real_T kkbkkd5z5m ; real_T gztmpcdsdg ; real_T mxnffxyshc ; real_T fr14joswus
; real_T g3vckbtrp3 ; real_T cqkdaw1iu3 ; real_T o123mdetox ; real_T
jsgd1kggkx ; real_T nc3ujqp1pg ; real_T iilt4zdhq4 ; real_T foabxlso54 ;
real_T mscp3mtqiu ; real_T orhxqsmvfm ; real_T pcceognrbb ; real_T hcfgbzdbcf
; real_T jost0dydiu ; real_T nodqar1hqn ; real_T oi30omicce ; real_T
nmkcufy13s ; real_T bhydlzvdgu ; real_T jnk41mcgi4 ; real_T nhtjb14mqm ;
real_T becosfp3hi ; real_T lphlstzov0 ; real_T dcpeged5b2 ; real_T eqjfsu4kur
; real_T cpcl24qdbb ; real_T hpzll4vow0 ; real_T okpqgwd3h2 ; real_T
cbpdlgxiw5 ; real_T lv3onhl4x1 ; real_T fmn4xw0byt ; real_T aeifxkbw2r ;
real_T idpdlwl5kt ; real_T haqtqgqlej ; real_T e0pcr2qr33 ; real_T fo05goe5wa
; real_T a1drd1cazf ; real_T lqpliiebzv ; real_T hrpirqqfri ; real_T
ollqjfiay1 ; real_T fzo2unh14f ; real_T msv1wjj5cz ; real_T lxsvfwftlj ;
real_T jjwnm1umsr ; real_T jbbgrda5ov ; real_T bojkujnbsq ; real_T bnvm1inada
; real_T dw11js5duv ; real_T jabohzh110 ; real_T iikt2yfv4y ; real_T
ilwj52njki ; real_T pwpjs2gjry ; real_T b3cdqkxcij ; } B ; typedef struct {
real_T mog1ypd0na ; struct { void * LoggedData [ 2 ] ; } mjjnlupcqc ; struct
{ void * LoggedData ; } m5bl3o4gs5 ; struct { void * LoggedData [ 4 ] ; }
pojnsaofwq ; struct { void * LoggedData ; } ajsjglxjta ; struct { void *
LoggedData [ 3 ] ; } ahihqvv0qg ; struct { void * LoggedData ; } lftwng3afb ;
struct { void * LoggedData ; } dgbo00d2tw ; struct { void * AQHandles ; }
cpa2ih4ymh ; struct { void * AQHandles ; } igylchmu0t ; struct { void *
AQHandles ; } aqjm4jt1o5 ; struct { void * AQHandles ; } gzhbrbv4c4 ; struct
{ void * AQHandles ; } htile1lzbu ; struct { void * AQHandles ; } p1ayluai4b
; struct { void * AQHandles ; } gjtuyi3wgw ; struct { void * AQHandles ; }
jd4fubbir0 ; struct { void * AQHandles ; } kbuynfluan ; struct { void *
AQHandles ; } j5pyh2o05c ; struct { void * AQHandles ; } p0fxbp34np ; struct
{ void * AQHandles ; } n4qxc3kzkl ; struct { void * LoggedData ; } mma0ea0nbs
; int32_T d1ijg3op2k ; int32_T brka54cpul ; int32_T ojcciswcd3 ; int32_T
nfkviszx2j ; int32_T ndyrr3tae0 ; int32_T c2jtsctdjh ; int32_T mevs4smziy ;
int32_T egiikq3opb ; int32_T puivctwyq5 ; int32_T bx2j5ck3fy ; int32_T
fzispyfeqr ; int32_T brdpyuefrq ; int32_T esun0jx3xw ; int32_T ncj5ml0ihc ;
int32_T kjn3oyfrac ; int32_T d3q0cyta13 ; int_T di0xkhax5v ; int_T dol3avfpna
; int_T mfw5n1ea1y ; int_T oj5q0jsluk ; int_T dfi4io44l1 ; int_T bapvsmhne1 ;
int_T gmzfgv4p33 ; int_T bxfqhgcloy ; int_T gwwzsdaztz ; int_T agicmdjqdi ;
int8_T astl0qyquk ; uint8_T ekvnfwjtsb ; uint8_T dql53suw4t ; uint8_T
gz0frs1wku ; uint8_T iwrzpan53u ; uint8_T ivbuycouap ; uint8_T oatrzraqor ;
uint8_T h3yqjmxcau ; uint8_T c2gp1hcawr ; uint8_T kzezuc5kmj ; uint8_T
baej0t4mz4 ; uint8_T e4o1f3himj ; uint8_T etn0cabray ; uint8_T nq5efgjr4u ;
uint8_T ldwj0svuv5 ; uint8_T et5jt4roea ; uint8_T jhgzfejyh1 ; boolean_T
e3b1amt0gc ; boolean_T g3swjv3kh4 ; boolean_T cifi3jtean ; boolean_T
mx1mzmmudm ; boolean_T nhmw0vweas ; boolean_T enb1yfnvht ; boolean_T
au1n1jxucn ; boolean_T cwu3pvv3pv ; boolean_T dwzumrr4bn ; boolean_T
jelszbraal ; boolean_T cbdxrybiau ; boolean_T dia2paiwqh ; boolean_T
f40yghkee4 ; boolean_T pke0ukp1xy ; boolean_T ob1pklob1o ; boolean_T
jvf3pz3ask ; boolean_T cy4ormxj51 ; } DW ; typedef struct { real_T njnhyztmul
; real_T cxtidlrvgx ; real_T anhq5q4adb ; real_T h1nzd5kpbz ; real_T
pmj0r0j4v2 ; real_T dusu5owudw ; real_T decqh12jsd ; real_T f3sxyupuw4 ;
real_T javki5lnyj ; real_T d1tt3gwruv ; real_T erk0rjtwtw ; real_T m2od42nnln
; real_T ith00xbmqs ; real_T djreaptdva ; real_T pnlrwk52wo ; real_T
lxsclehbdn [ 2 ] ; real_T herd5wji3e [ 2 ] ; real_T ekhpkuwiae [ 2 ] ; real_T
k2eyyio13c [ 2 ] ; real_T mh2ntyj1vq ; real_T kh4uwqwwas ; real_T l1m0lcksjp
[ 2 ] ; } X ; typedef struct { real_T njnhyztmul ; real_T cxtidlrvgx ; real_T
anhq5q4adb ; real_T h1nzd5kpbz ; real_T pmj0r0j4v2 ; real_T dusu5owudw ;
real_T decqh12jsd ; real_T f3sxyupuw4 ; real_T javki5lnyj ; real_T d1tt3gwruv
; real_T erk0rjtwtw ; real_T m2od42nnln ; real_T ith00xbmqs ; real_T
djreaptdva ; real_T pnlrwk52wo ; real_T lxsclehbdn [ 2 ] ; real_T herd5wji3e
[ 2 ] ; real_T ekhpkuwiae [ 2 ] ; real_T k2eyyio13c [ 2 ] ; real_T mh2ntyj1vq
; real_T kh4uwqwwas ; real_T l1m0lcksjp [ 2 ] ; } XDot ; typedef struct {
boolean_T njnhyztmul ; boolean_T cxtidlrvgx ; boolean_T anhq5q4adb ;
boolean_T h1nzd5kpbz ; boolean_T pmj0r0j4v2 ; boolean_T dusu5owudw ;
boolean_T decqh12jsd ; boolean_T f3sxyupuw4 ; boolean_T javki5lnyj ;
boolean_T d1tt3gwruv ; boolean_T erk0rjtwtw ; boolean_T m2od42nnln ;
boolean_T ith00xbmqs ; boolean_T djreaptdva ; boolean_T pnlrwk52wo ;
boolean_T lxsclehbdn [ 2 ] ; boolean_T herd5wji3e [ 2 ] ; boolean_T
ekhpkuwiae [ 2 ] ; boolean_T k2eyyio13c [ 2 ] ; boolean_T mh2ntyj1vq ;
boolean_T kh4uwqwwas ; boolean_T l1m0lcksjp [ 2 ] ; } XDis ; typedef struct {
real_T njnhyztmul ; real_T cxtidlrvgx ; real_T anhq5q4adb ; real_T h1nzd5kpbz
; real_T pmj0r0j4v2 ; real_T dusu5owudw ; real_T decqh12jsd ; real_T
f3sxyupuw4 ; real_T javki5lnyj ; real_T d1tt3gwruv ; real_T erk0rjtwtw ;
real_T m2od42nnln ; real_T ith00xbmqs ; real_T djreaptdva ; real_T pnlrwk52wo
; real_T lxsclehbdn [ 2 ] ; real_T herd5wji3e [ 2 ] ; real_T ekhpkuwiae [ 2 ]
; real_T k2eyyio13c [ 2 ] ; real_T mh2ntyj1vq ; real_T kh4uwqwwas ; real_T
l1m0lcksjp [ 2 ] ; } CStateAbsTol ; typedef struct { real_T njnhyztmul ;
real_T cxtidlrvgx ; real_T anhq5q4adb ; real_T h1nzd5kpbz ; real_T pmj0r0j4v2
; real_T dusu5owudw ; real_T decqh12jsd ; real_T f3sxyupuw4 ; real_T
javki5lnyj ; real_T d1tt3gwruv ; real_T erk0rjtwtw ; real_T m2od42nnln ;
real_T ith00xbmqs ; real_T djreaptdva ; real_T pnlrwk52wo ; real_T lxsclehbdn
[ 2 ] ; real_T herd5wji3e [ 2 ] ; real_T ekhpkuwiae [ 2 ] ; real_T k2eyyio13c
[ 2 ] ; real_T mh2ntyj1vq ; real_T kh4uwqwwas ; real_T l1m0lcksjp [ 2 ] ; }
CXPtMin ; typedef struct { real_T njnhyztmul ; real_T cxtidlrvgx ; real_T
anhq5q4adb ; real_T h1nzd5kpbz ; real_T pmj0r0j4v2 ; real_T dusu5owudw ;
real_T decqh12jsd ; real_T f3sxyupuw4 ; real_T javki5lnyj ; real_T d1tt3gwruv
; real_T erk0rjtwtw ; real_T m2od42nnln ; real_T ith00xbmqs ; real_T
djreaptdva ; real_T pnlrwk52wo ; real_T lxsclehbdn [ 2 ] ; real_T herd5wji3e
[ 2 ] ; real_T ekhpkuwiae [ 2 ] ; real_T k2eyyio13c [ 2 ] ; real_T mh2ntyj1vq
; real_T kh4uwqwwas ; real_T l1m0lcksjp [ 2 ] ; } CXPtMax ; typedef struct {
real_T n52axficl2 ; real_T a15z3evljq ; real_T fbpmdji415 ; real_T odfwjcwoxs
; real_T fdgm22ouln ; real_T haaywk4krz ; real_T l2jklhprzm ; real_T
jdcwcq0rv4 ; real_T dqejsrbctn ; real_T fs30lsvp00 ; real_T jy5xmnlgzp ;
real_T daixzoq4vk ; real_T cueld2yjzo ; } ZCV ; typedef struct { ZCSigState
ii2vsk3a4a ; } PrevZCX ; typedef struct { rtwCAPI_ModelMappingInfo mmi ; }
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
mr_trial_with_intermediate_tanks_GetDWork ( ) ; extern void
mr_trial_with_intermediate_tanks_SetDWork ( const mxArray * ssDW ) ; extern
mxArray * mr_trial_with_intermediate_tanks_GetSimStateDisallowedBlocks ( ) ;
extern const rtwCAPI_ModelMappingStaticInfo *
trial_with_intermediate_tanks_GetCAPIStaticMap ( void ) ; extern SimStruct *
const rtS ; extern DataMapInfo * rt_dataMapInfoPtr ; extern
rtwCAPI_ModelMappingInfo * rt_modelMapInfoPtr ; void MdlOutputs ( int_T tid )
; void MdlOutputsParameterSampleTime ( int_T tid ) ; void MdlUpdate ( int_T
tid ) ; void MdlTerminate ( void ) ; void MdlInitializeSizes ( void ) ; void
MdlInitializeSampleTimes ( void ) ; SimStruct * raccel_register_model (
ssExecutionInfo * executionInfo ) ;
#endif
