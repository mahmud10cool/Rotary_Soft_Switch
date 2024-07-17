#include "rtw_capi.h"
#ifdef HOST_CAPI_BUILD
#include "trial_with_intermediate_tanks_capi_host.h"
#define sizeof(s) ((size_t)(0xFFFF))
#undef rt_offsetof
#define rt_offsetof(s,el) ((uint16_T)(0xFFFF))
#define TARGET_CONST
#define TARGET_STRING(s) (s)
#ifndef SS_UINT64
#define SS_UINT64 18
#endif
#ifndef SS_INT64
#define SS_INT64 19
#endif
#else
#include "builtin_typeid_types.h"
#include "trial_with_intermediate_tanks.h"
#include "trial_with_intermediate_tanks_capi.h"
#include "trial_with_intermediate_tanks_private.h"
#ifdef LIGHT_WEIGHT_CAPI
#define TARGET_CONST
#define TARGET_STRING(s)               ((NULL))
#else
#define TARGET_CONST                   const
#define TARGET_STRING(s)               (s)
#endif
#endif
static const rtwCAPI_Signals rtBlockSignals [ ] = { { 0 , 13 , TARGET_STRING
( "trial_with_intermediate_tanks/Triggered Subsystem" ) , TARGET_STRING (
"Event Time" ) , 0 , 0 , 0 , 0 , 0 } , { 1 , 1 , TARGET_STRING (
"trial_with_intermediate_tanks/Electric Torque Controller" ) , TARGET_STRING
( "" ) , 0 , 0 , 0 , 0 , 1 } , { 2 , 0 , TARGET_STRING (
 "trial_with_intermediate_tanks/Electric Torque Controller/is_active_c11_trial_with_intermediate_tanks"
) , TARGET_STRING ( "is_active_c11_trial_with_intermediate_tanks" ) , 0 , 1 ,
0 , 0 , 1 } , { 3 , 2 , TARGET_STRING (
"trial_with_intermediate_tanks/Hydraulic Pump//Motor" ) , TARGET_STRING ( ""
) , 0 , 0 , 0 , 0 , 1 } , { 4 , 0 , TARGET_STRING (
 "trial_with_intermediate_tanks/Hydraulic Pump//Motor/is_active_c9_trial_with_intermediate_tanks"
) , TARGET_STRING ( "is_active_c9_trial_with_intermediate_tanks" ) , 0 , 1 ,
0 , 0 , 1 } , { 5 , 3 , TARGET_STRING (
"trial_with_intermediate_tanks/Hydraulic Pump//Motor Dynamics" ) ,
TARGET_STRING ( "omega_dot" ) , 0 , 0 , 0 , 0 , 1 } , { 6 , 3 , TARGET_STRING
( "trial_with_intermediate_tanks/Hydraulic Pump//Motor Dynamics" ) ,
TARGET_STRING ( "" ) , 1 , 0 , 0 , 0 , 1 } , { 7 , 0 , TARGET_STRING (
 "trial_with_intermediate_tanks/Hydraulic Pump//Motor Dynamics/is_active_c10_trial_with_intermediate_tanks"
) , TARGET_STRING ( "is_active_c10_trial_with_intermediate_tanks" ) , 0 , 1 ,
0 , 0 , 1 } , { 8 , 4 , TARGET_STRING (
"trial_with_intermediate_tanks/Linear Actuator Chamber" ) , TARGET_STRING (
"" ) , 0 , 0 , 0 , 0 , 1 } , { 9 , 0 , TARGET_STRING (
 "trial_with_intermediate_tanks/Linear Actuator Chamber/is_active_c8_trial_with_intermediate_tanks"
) , TARGET_STRING ( "is_active_c8_trial_with_intermediate_tanks" ) , 0 , 1 ,
0 , 0 , 1 } , { 10 , 5 , TARGET_STRING (
"trial_with_intermediate_tanks/MATLAB Function1" ) , TARGET_STRING (
"Stored Energy 3" ) , 0 , 0 , 0 , 0 , 1 } , { 11 , 0 , TARGET_STRING (
 "trial_with_intermediate_tanks/MATLAB Function1/is_active_c13_trial_with_intermediate_tanks"
) , TARGET_STRING ( "is_active_c13_trial_with_intermediate_tanks" ) , 0 , 1 ,
0 , 0 , 1 } , { 12 , 6 , TARGET_STRING (
"trial_with_intermediate_tanks/MATLAB Function2" ) , TARGET_STRING ( "" ) , 0
, 0 , 0 , 0 , 2 } , { 13 , 6 , TARGET_STRING (
"trial_with_intermediate_tanks/MATLAB Function2" ) , TARGET_STRING ( "" ) , 1
, 0 , 0 , 0 , 2 } , { 14 , 6 , TARGET_STRING (
"trial_with_intermediate_tanks/MATLAB Function2" ) , TARGET_STRING ( "" ) , 2
, 0 , 0 , 0 , 2 } , { 15 , 6 , TARGET_STRING (
"trial_with_intermediate_tanks/MATLAB Function2" ) , TARGET_STRING ( "" ) , 3
, 0 , 0 , 0 , 2 } , { 16 , 0 , TARGET_STRING (
 "trial_with_intermediate_tanks/MATLAB Function2/is_active_c3_trial_with_intermediate_tanks"
) , TARGET_STRING ( "is_active_c3_trial_with_intermediate_tanks" ) , 0 , 1 ,
0 , 0 , 2 } , { 17 , 7 , TARGET_STRING (
"trial_with_intermediate_tanks/MATLAB Function3" ) , TARGET_STRING (
"Stored Energy 4" ) , 0 , 0 , 0 , 0 , 1 } , { 18 , 0 , TARGET_STRING (
 "trial_with_intermediate_tanks/MATLAB Function3/is_active_c12_trial_with_intermediate_tanks"
) , TARGET_STRING ( "is_active_c12_trial_with_intermediate_tanks" ) , 0 , 1 ,
0 , 0 , 1 } , { 19 , 8 , TARGET_STRING (
"trial_with_intermediate_tanks/MATLAB Function4" ) , TARGET_STRING (
"Power_out" ) , 0 , 0 , 0 , 0 , 1 } , { 20 , 0 , TARGET_STRING (
 "trial_with_intermediate_tanks/MATLAB Function4/is_active_c14_trial_with_intermediate_tanks"
) , TARGET_STRING ( "is_active_c14_trial_with_intermediate_tanks" ) , 0 , 1 ,
0 , 0 , 1 } , { 21 , 9 , TARGET_STRING (
"trial_with_intermediate_tanks/MATLAB Function5" ) , TARGET_STRING (
"Stored Energy 3" ) , 0 , 0 , 0 , 0 , 1 } , { 22 , 0 , TARGET_STRING (
 "trial_with_intermediate_tanks/MATLAB Function5/is_active_c15_trial_with_intermediate_tanks"
) , TARGET_STRING ( "is_active_c15_trial_with_intermediate_tanks" ) , 0 , 1 ,
0 , 0 , 1 } , { 23 , 10 , TARGET_STRING (
"trial_with_intermediate_tanks/MATLAB Function6" ) , TARGET_STRING (
"Stored Energy 4" ) , 0 , 0 , 0 , 0 , 1 } , { 24 , 0 , TARGET_STRING (
 "trial_with_intermediate_tanks/MATLAB Function6/is_active_c16_trial_with_intermediate_tanks"
) , TARGET_STRING ( "is_active_c16_trial_with_intermediate_tanks" ) , 0 , 1 ,
0 , 0 , 1 } , { 25 , 11 , TARGET_STRING (
"trial_with_intermediate_tanks/Small Chamber 3" ) , TARGET_STRING ( "" ) , 0
, 0 , 0 , 0 , 1 } , { 26 , 0 , TARGET_STRING (
 "trial_with_intermediate_tanks/Small Chamber 3/is_active_c1_trial_with_intermediate_tanks"
) , TARGET_STRING ( "is_active_c1_trial_with_intermediate_tanks" ) , 0 , 1 ,
0 , 0 , 1 } , { 27 , 12 , TARGET_STRING (
"trial_with_intermediate_tanks/Small Chamber 4" ) , TARGET_STRING ( "" ) , 0
, 0 , 0 , 0 , 1 } , { 28 , 0 , TARGET_STRING (
 "trial_with_intermediate_tanks/Small Chamber 4/is_active_c7_trial_with_intermediate_tanks"
) , TARGET_STRING ( "is_active_c7_trial_with_intermediate_tanks" ) , 0 , 1 ,
0 , 0 , 1 } , { 29 , 14 , TARGET_STRING (
"trial_with_intermediate_tanks/Valve A" ) , TARGET_STRING ( "" ) , 0 , 0 , 0
, 0 , 1 } , { 30 , 0 , TARGET_STRING (
 "trial_with_intermediate_tanks/Valve A/is_active_c2_trial_with_intermediate_tanks"
) , TARGET_STRING ( "is_active_c2_trial_with_intermediate_tanks" ) , 0 , 1 ,
0 , 0 , 1 } , { 31 , 15 , TARGET_STRING (
"trial_with_intermediate_tanks/Valve B" ) , TARGET_STRING ( "" ) , 0 , 0 , 0
, 0 , 1 } , { 32 , 0 , TARGET_STRING (
 "trial_with_intermediate_tanks/Valve B/is_active_c4_trial_with_intermediate_tanks"
) , TARGET_STRING ( "is_active_c4_trial_with_intermediate_tanks" ) , 0 , 1 ,
0 , 0 , 1 } , { 33 , 16 , TARGET_STRING (
"trial_with_intermediate_tanks/Valve C" ) , TARGET_STRING ( "" ) , 0 , 0 , 0
, 0 , 1 } , { 34 , 0 , TARGET_STRING (
 "trial_with_intermediate_tanks/Valve C/is_active_c6_trial_with_intermediate_tanks"
) , TARGET_STRING ( "is_active_c6_trial_with_intermediate_tanks" ) , 0 , 1 ,
0 , 0 , 1 } , { 35 , 17 , TARGET_STRING (
"trial_with_intermediate_tanks/Valve D" ) , TARGET_STRING ( "" ) , 0 , 0 , 0
, 0 , 1 } , { 36 , 0 , TARGET_STRING (
 "trial_with_intermediate_tanks/Valve D/is_active_c5_trial_with_intermediate_tanks"
) , TARGET_STRING ( "is_active_c5_trial_with_intermediate_tanks" ) , 0 , 1 ,
0 , 0 , 1 } , { 37 , 0 , TARGET_STRING (
"trial_with_intermediate_tanks/Clock" ) , TARGET_STRING ( "Simulation Time" )
, 0 , 0 , 0 , 0 , 1 } , { 38 , 0 , TARGET_STRING (
"trial_with_intermediate_tanks/Constant1" ) , TARGET_STRING ( "" ) , 0 , 0 ,
0 , 0 , 3 } , { 39 , 0 , TARGET_STRING (
"trial_with_intermediate_tanks/Constant3" ) , TARGET_STRING ( "" ) , 0 , 0 ,
0 , 0 , 3 } , { 40 , 0 , TARGET_STRING (
"trial_with_intermediate_tanks/Constant4" ) , TARGET_STRING ( "" ) , 0 , 0 ,
0 , 0 , 3 } , { 41 , 0 , TARGET_STRING (
"trial_with_intermediate_tanks/Initial Speed" ) , TARGET_STRING ( "" ) , 0 ,
0 , 0 , 0 , 3 } , { 42 , 0 , TARGET_STRING (
"trial_with_intermediate_tanks/Gain" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0
, 1 } , { 43 , 0 , TARGET_STRING ( "trial_with_intermediate_tanks/Gain1" ) ,
TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 1 } , { 44 , 0 , TARGET_STRING (
"trial_with_intermediate_tanks/Gain3" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 ,
0 , 1 } , { 45 , 0 , TARGET_STRING ( "trial_with_intermediate_tanks/Gain5" )
, TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 3 } , { 46 , 0 , TARGET_STRING (
"trial_with_intermediate_tanks/Integrator" ) , TARGET_STRING ( "P3" ) , 0 , 0
, 0 , 0 , 1 } , { 47 , 0 , TARGET_STRING (
"trial_with_intermediate_tanks/Integrator1" ) , TARGET_STRING ( "P3" ) , 0 ,
0 , 0 , 0 , 1 } , { 48 , 0 , TARGET_STRING (
"trial_with_intermediate_tanks/Integrator10" ) , TARGET_STRING ( "" ) , 0 , 0
, 0 , 0 , 1 } , { 49 , 0 , TARGET_STRING (
"trial_with_intermediate_tanks/Integrator11" ) , TARGET_STRING ( "" ) , 0 , 0
, 0 , 0 , 1 } , { 50 , 0 , TARGET_STRING (
"trial_with_intermediate_tanks/Integrator12" ) , TARGET_STRING ( "" ) , 0 , 0
, 0 , 0 , 1 } , { 51 , 0 , TARGET_STRING (
"trial_with_intermediate_tanks/Integrator13" ) , TARGET_STRING ( "" ) , 0 , 0
, 0 , 0 , 1 } , { 52 , 0 , TARGET_STRING (
"trial_with_intermediate_tanks/Integrator14" ) , TARGET_STRING ( "" ) , 0 , 0
, 0 , 0 , 1 } , { 53 , 0 , TARGET_STRING (
"trial_with_intermediate_tanks/Integrator15" ) , TARGET_STRING (
"Displacement" ) , 0 , 0 , 0 , 0 , 1 } , { 54 , 0 , TARGET_STRING (
"trial_with_intermediate_tanks/Integrator16" ) , TARGET_STRING ( "" ) , 0 , 0
, 0 , 0 , 1 } , { 55 , 0 , TARGET_STRING (
"trial_with_intermediate_tanks/Integrator2" ) , TARGET_STRING ( "P3" ) , 0 ,
0 , 0 , 0 , 1 } , { 56 , 0 , TARGET_STRING (
"trial_with_intermediate_tanks/Integrator3" ) , TARGET_STRING ( "omega" ) , 0
, 0 , 0 , 0 , 1 } , { 57 , 0 , TARGET_STRING (
"trial_with_intermediate_tanks/Integrator4" ) , TARGET_STRING ( "" ) , 0 , 0
, 0 , 0 , 1 } , { 58 , 0 , TARGET_STRING (
"trial_with_intermediate_tanks/Integrator5" ) , TARGET_STRING ( "" ) , 0 , 0
, 0 , 0 , 1 } , { 59 , 0 , TARGET_STRING (
"trial_with_intermediate_tanks/Integrator6" ) , TARGET_STRING ( "" ) , 0 , 0
, 0 , 0 , 1 } , { 60 , 0 , TARGET_STRING (
"trial_with_intermediate_tanks/Integrator7" ) , TARGET_STRING ( "" ) , 0 , 0
, 0 , 0 , 1 } , { 61 , 0 , TARGET_STRING (
"trial_with_intermediate_tanks/Integrator8" ) , TARGET_STRING ( "" ) , 0 , 0
, 0 , 0 , 1 } , { 62 , 0 , TARGET_STRING (
"trial_with_intermediate_tanks/Integrator9" ) , TARGET_STRING ( "" ) , 0 , 0
, 0 , 0 , 1 } , { 63 , 0 , TARGET_STRING (
"trial_with_intermediate_tanks/Product" ) , TARGET_STRING ( "" ) , 0 , 0 , 0
, 0 , 1 } , { 64 , 0 , TARGET_STRING (
"trial_with_intermediate_tanks/Product1" ) , TARGET_STRING ( "" ) , 0 , 0 , 0
, 0 , 1 } , { 65 , 0 , TARGET_STRING (
"trial_with_intermediate_tanks/Product2" ) , TARGET_STRING ( "" ) , 0 , 0 , 0
, 0 , 1 } , { 66 , 0 , TARGET_STRING (
"trial_with_intermediate_tanks/Product3" ) , TARGET_STRING ( "" ) , 0 , 0 , 0
, 0 , 1 } , { 67 , 0 , TARGET_STRING (
"trial_with_intermediate_tanks/Product4" ) , TARGET_STRING ( "" ) , 0 , 0 , 0
, 0 , 1 } , { 68 , 0 , TARGET_STRING (
"trial_with_intermediate_tanks/Product5" ) , TARGET_STRING ( "" ) , 0 , 0 , 0
, 0 , 1 } , { 69 , 0 , TARGET_STRING (
"trial_with_intermediate_tanks/Product6" ) , TARGET_STRING ( "" ) , 0 , 0 , 0
, 0 , 1 } , { 70 , 0 , TARGET_STRING (
"trial_with_intermediate_tanks/Product7" ) , TARGET_STRING ( "" ) , 0 , 0 , 0
, 0 , 1 } , { 71 , 0 , TARGET_STRING ( "trial_with_intermediate_tanks/Add" )
, TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 1 } , { 72 , 0 , TARGET_STRING (
"trial_with_intermediate_tanks/Add1" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0
, 1 } , { 73 , 0 , TARGET_STRING ( "trial_with_intermediate_tanks/Add3" ) ,
TARGET_STRING ( "Out" ) , 0 , 0 , 0 , 0 , 1 } , { 74 , 0 , TARGET_STRING (
"trial_with_intermediate_tanks/Add4" ) , TARGET_STRING ( "In" ) , 0 , 0 , 0 ,
0 , 1 } , { 75 , 0 , TARGET_STRING ( "trial_with_intermediate_tanks/Subtract"
) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 1 } , { 76 , 0 , TARGET_STRING (
"trial_with_intermediate_tanks/Subtract1" ) , TARGET_STRING ( "" ) , 0 , 0 ,
0 , 0 , 1 } , { 77 , 0 , TARGET_STRING (
"trial_with_intermediate_tanks/Subtract2" ) , TARGET_STRING ( "" ) , 0 , 0 ,
0 , 0 , 1 } , { 78 , 0 , TARGET_STRING (
"trial_with_intermediate_tanks/Subtract3" ) , TARGET_STRING ( "" ) , 0 , 0 ,
0 , 0 , 1 } , { 79 , 0 , TARGET_STRING (
"trial_with_intermediate_tanks/Subtract4" ) , TARGET_STRING ( "" ) , 0 , 0 ,
0 , 0 , 1 } , { 80 , 0 , TARGET_STRING (
"trial_with_intermediate_tanks/Manual Switch" ) , TARGET_STRING ( "" ) , 0 ,
0 , 0 , 0 , 1 } , { 81 , 0 , TARGET_STRING (
"trial_with_intermediate_tanks/Delay" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 ,
0 , 4 } , { 82 , 0 , TARGET_STRING (
"trial_with_intermediate_tanks/Dynamics of Valve/Saturation" ) ,
TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 1 } , { 83 , 0 , TARGET_STRING (
"trial_with_intermediate_tanks/Dynamics of Valve/Transfer Fcn" ) ,
TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 1 } , { 84 , 0 , TARGET_STRING (
"trial_with_intermediate_tanks/Dynamics of Valve1/Saturation" ) ,
TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 1 } , { 85 , 0 , TARGET_STRING (
"trial_with_intermediate_tanks/Dynamics of Valve1/Transfer Fcn" ) ,
TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 1 } , { 86 , 0 , TARGET_STRING (
"trial_with_intermediate_tanks/Dynamics of Valve2/Saturation" ) ,
TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 1 } , { 87 , 0 , TARGET_STRING (
"trial_with_intermediate_tanks/Dynamics of Valve2/Transfer Fcn" ) ,
TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 1 } , { 88 , 0 , TARGET_STRING (
"trial_with_intermediate_tanks/Dynamics of Valve3/Saturation" ) ,
TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 1 } , { 89 , 0 , TARGET_STRING (
"trial_with_intermediate_tanks/Dynamics of Valve3/Transfer Fcn" ) ,
TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 1 } , { 90 , 0 , TARGET_STRING (
"trial_with_intermediate_tanks/Dynamics of Valve4/Saturation" ) ,
TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 1 } , { 91 , 0 , TARGET_STRING (
"trial_with_intermediate_tanks/Dynamics of Valve4/Transfer Fcn" ) ,
TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 1 } , { 92 , 0 , TARGET_STRING (
"trial_with_intermediate_tanks/Subsystem/Abs" ) , TARGET_STRING ( "" ) , 0 ,
0 , 0 , 0 , 1 } , { 93 , 0 , TARGET_STRING (
"trial_with_intermediate_tanks/Subsystem/Relay" ) , TARGET_STRING ( "" ) , 0
, 0 , 0 , 0 , 2 } , { 94 , 0 , TARGET_STRING (
"trial_with_intermediate_tanks/Subsystem/Subtract" ) , TARGET_STRING ( "" ) ,
0 , 0 , 0 , 0 , 1 } , { 95 , 13 , TARGET_STRING (
"trial_with_intermediate_tanks/Triggered Subsystem/In1" ) , TARGET_STRING (
"" ) , 0 , 0 , 0 , 0 , 0 } , { 0 , 0 , ( NULL ) , ( NULL ) , 0 , 0 , 0 , 0 ,
0 } } ; static const rtwCAPI_BlockParameters rtBlockParameters [ ] = { { 96 ,
TARGET_STRING ( "trial_with_intermediate_tanks/Constant" ) , TARGET_STRING (
"Value" ) , 0 , 0 , 0 } , { 97 , TARGET_STRING (
"trial_with_intermediate_tanks/Constant5" ) , TARGET_STRING ( "Value" ) , 0 ,
0 , 0 } , { 98 , TARGET_STRING (
"trial_with_intermediate_tanks/Initial Speed" ) , TARGET_STRING ( "Value" ) ,
0 , 0 , 0 } , { 99 , TARGET_STRING ( "trial_with_intermediate_tanks/Gain" ) ,
TARGET_STRING ( "Gain" ) , 0 , 0 , 0 } , { 100 , TARGET_STRING (
"trial_with_intermediate_tanks/Gain3" ) , TARGET_STRING ( "Gain" ) , 0 , 0 ,
0 } , { 101 , TARGET_STRING ( "trial_with_intermediate_tanks/Gain5" ) ,
TARGET_STRING ( "Gain" ) , 0 , 0 , 0 } , { 102 , TARGET_STRING (
"trial_with_intermediate_tanks/Integrator10" ) , TARGET_STRING (
"InitialCondition" ) , 0 , 0 , 0 } , { 103 , TARGET_STRING (
"trial_with_intermediate_tanks/Integrator11" ) , TARGET_STRING (
"InitialCondition" ) , 0 , 0 , 0 } , { 104 , TARGET_STRING (
"trial_with_intermediate_tanks/Integrator12" ) , TARGET_STRING (
"InitialCondition" ) , 0 , 0 , 0 } , { 105 , TARGET_STRING (
"trial_with_intermediate_tanks/Integrator13" ) , TARGET_STRING (
"InitialCondition" ) , 0 , 0 , 0 } , { 106 , TARGET_STRING (
"trial_with_intermediate_tanks/Integrator14" ) , TARGET_STRING (
"InitialCondition" ) , 0 , 0 , 0 } , { 107 , TARGET_STRING (
"trial_with_intermediate_tanks/Integrator15" ) , TARGET_STRING (
"InitialCondition" ) , 0 , 0 , 0 } , { 108 , TARGET_STRING (
"trial_with_intermediate_tanks/Integrator16" ) , TARGET_STRING (
"InitialCondition" ) , 0 , 0 , 0 } , { 109 , TARGET_STRING (
"trial_with_intermediate_tanks/Integrator4" ) , TARGET_STRING (
"InitialCondition" ) , 0 , 0 , 0 } , { 110 , TARGET_STRING (
"trial_with_intermediate_tanks/Integrator5" ) , TARGET_STRING (
"InitialCondition" ) , 0 , 0 , 0 } , { 111 , TARGET_STRING (
"trial_with_intermediate_tanks/Integrator6" ) , TARGET_STRING (
"InitialCondition" ) , 0 , 0 , 0 } , { 112 , TARGET_STRING (
"trial_with_intermediate_tanks/Integrator7" ) , TARGET_STRING (
"InitialCondition" ) , 0 , 0 , 0 } , { 113 , TARGET_STRING (
"trial_with_intermediate_tanks/Integrator8" ) , TARGET_STRING (
"InitialCondition" ) , 0 , 0 , 0 } , { 114 , TARGET_STRING (
"trial_with_intermediate_tanks/Integrator9" ) , TARGET_STRING (
"InitialCondition" ) , 0 , 0 , 0 } , { 115 , TARGET_STRING (
"trial_with_intermediate_tanks/Velocity" ) , TARGET_STRING ( "Amplitude" ) ,
0 , 0 , 0 } , { 116 , TARGET_STRING (
"trial_with_intermediate_tanks/Velocity" ) , TARGET_STRING ( "Bias" ) , 0 , 0
, 0 } , { 117 , TARGET_STRING ( "trial_with_intermediate_tanks/Velocity" ) ,
TARGET_STRING ( "Frequency" ) , 0 , 0 , 0 } , { 118 , TARGET_STRING (
"trial_with_intermediate_tanks/Velocity" ) , TARGET_STRING ( "Phase" ) , 0 ,
0 , 0 } , { 119 , TARGET_STRING (
"trial_with_intermediate_tanks/Manual Switch" ) , TARGET_STRING (
"CurrentSetting" ) , 1 , 0 , 0 } , { 120 , TARGET_STRING (
"trial_with_intermediate_tanks/Delay" ) , TARGET_STRING ( "InitialCondition"
) , 0 , 0 , 0 } , { 121 , TARGET_STRING (
"trial_with_intermediate_tanks/Dynamics of Valve/Saturation" ) ,
TARGET_STRING ( "LowerLimit" ) , 0 , 0 , 0 } , { 122 , TARGET_STRING (
"trial_with_intermediate_tanks/Dynamics of Valve/Transfer Fcn" ) ,
TARGET_STRING ( "A" ) , 0 , 1 , 0 } , { 123 , TARGET_STRING (
"trial_with_intermediate_tanks/Dynamics of Valve/Transfer Fcn" ) ,
TARGET_STRING ( "C" ) , 0 , 2 , 0 } , { 124 , TARGET_STRING (
"trial_with_intermediate_tanks/Dynamics of Valve1/Saturation" ) ,
TARGET_STRING ( "LowerLimit" ) , 0 , 0 , 0 } , { 125 , TARGET_STRING (
"trial_with_intermediate_tanks/Dynamics of Valve1/Transfer Fcn" ) ,
TARGET_STRING ( "A" ) , 0 , 1 , 0 } , { 126 , TARGET_STRING (
"trial_with_intermediate_tanks/Dynamics of Valve1/Transfer Fcn" ) ,
TARGET_STRING ( "C" ) , 0 , 2 , 0 } , { 127 , TARGET_STRING (
"trial_with_intermediate_tanks/Dynamics of Valve2/Saturation" ) ,
TARGET_STRING ( "LowerLimit" ) , 0 , 0 , 0 } , { 128 , TARGET_STRING (
"trial_with_intermediate_tanks/Dynamics of Valve2/Transfer Fcn" ) ,
TARGET_STRING ( "A" ) , 0 , 1 , 0 } , { 129 , TARGET_STRING (
"trial_with_intermediate_tanks/Dynamics of Valve2/Transfer Fcn" ) ,
TARGET_STRING ( "C" ) , 0 , 2 , 0 } , { 130 , TARGET_STRING (
"trial_with_intermediate_tanks/Dynamics of Valve3/Saturation" ) ,
TARGET_STRING ( "LowerLimit" ) , 0 , 0 , 0 } , { 131 , TARGET_STRING (
"trial_with_intermediate_tanks/Dynamics of Valve3/Transfer Fcn" ) ,
TARGET_STRING ( "A" ) , 0 , 1 , 0 } , { 132 , TARGET_STRING (
"trial_with_intermediate_tanks/Dynamics of Valve3/Transfer Fcn" ) ,
TARGET_STRING ( "C" ) , 0 , 2 , 0 } , { 133 , TARGET_STRING (
"trial_with_intermediate_tanks/Dynamics of Valve4/Saturation" ) ,
TARGET_STRING ( "LowerLimit" ) , 0 , 0 , 0 } , { 134 , TARGET_STRING (
"trial_with_intermediate_tanks/Dynamics of Valve4/Transfer Fcn" ) ,
TARGET_STRING ( "A" ) , 0 , 1 , 0 } , { 135 , TARGET_STRING (
"trial_with_intermediate_tanks/Dynamics of Valve4/Transfer Fcn" ) ,
TARGET_STRING ( "C" ) , 0 , 2 , 0 } , { 136 , TARGET_STRING (
"trial_with_intermediate_tanks/Subsystem/Relay" ) , TARGET_STRING (
"OnSwitchValue" ) , 0 , 0 , 0 } , { 137 , TARGET_STRING (
"trial_with_intermediate_tanks/Subsystem/Relay" ) , TARGET_STRING (
"OffSwitchValue" ) , 0 , 0 , 0 } , { 138 , TARGET_STRING (
"trial_with_intermediate_tanks/Subsystem/Relay" ) , TARGET_STRING (
"OnOutputValue" ) , 0 , 0 , 0 } , { 139 , TARGET_STRING (
"trial_with_intermediate_tanks/Subsystem/Relay" ) , TARGET_STRING (
"OffOutputValue" ) , 0 , 0 , 0 } , { 140 , TARGET_STRING (
"trial_with_intermediate_tanks/Triggered Subsystem/Out1" ) , TARGET_STRING (
"InitialOutput" ) , 0 , 0 , 0 } , { 0 , ( NULL ) , ( NULL ) , 0 , 0 , 0 } } ;
static int_T rt_LoggedStateIdxList [ ] = { - 1 } ; static const
rtwCAPI_Signals rtRootInputs [ ] = { { 0 , 0 , ( NULL ) , ( NULL ) , 0 , 0 ,
0 , 0 , 0 } } ; static const rtwCAPI_Signals rtRootOutputs [ ] = { { 0 , 0 ,
( NULL ) , ( NULL ) , 0 , 0 , 0 , 0 , 0 } } ; static const
rtwCAPI_ModelParameters rtModelParameters [ ] = { { 141 , TARGET_STRING (
"param" ) , 2 , 0 , 0 } , { 0 , ( NULL ) , 0 , 0 , 0 } } ;
#ifndef HOST_CAPI_BUILD
static void * rtDataAddrMap [ ] = { & rtB . a1drd1cazf , & rtB . b3cdqkxcij ,
& rtDW . jhgzfejyh1 , & rtB . pwpjs2gjry , & rtDW . et5jt4roea , & rtB .
iikt2yfv4y , & rtB . ilwj52njki , & rtDW . ldwj0svuv5 , & rtB . jabohzh110 ,
& rtDW . nq5efgjr4u , & rtB . dw11js5duv , & rtDW . etn0cabray , & rtB .
jjwnm1umsr , & rtB . jbbgrda5ov , & rtB . bojkujnbsq , & rtB . bnvm1inada , &
rtDW . e4o1f3himj , & rtB . lxsvfwftlj , & rtDW . baej0t4mz4 , & rtB .
msv1wjj5cz , & rtDW . kzezuc5kmj , & rtB . fzo2unh14f , & rtDW . c2gp1hcawr ,
& rtB . ollqjfiay1 , & rtDW . h3yqjmxcau , & rtB . hrpirqqfri , & rtDW .
oatrzraqor , & rtB . lqpliiebzv , & rtDW . ivbuycouap , & rtB . fo05goe5wa ,
& rtDW . iwrzpan53u , & rtB . e0pcr2qr33 , & rtDW . gz0frs1wku , & rtB .
haqtqgqlej , & rtDW . dql53suw4t , & rtB . idpdlwl5kt , & rtDW . ekvnfwjtsb ,
& rtB . cqkdaw1iu3 , & rtB . lv3onhl4x1 , & rtB . fmn4xw0byt , & rtB .
aeifxkbw2r , & rtB . okpqgwd3h2 , & rtB . frh4wmn3vf , & rtB . iilt4zdhq4 , &
rtB . ab5hqtcfdy , & rtB . cbpdlgxiw5 , & rtB . ikogzwrhkz , & rtB .
fakrdhsmg5 , & rtB . fkszch0qwe , & rtB . ku0vq41nvc , & rtB . axh1siq1cq , &
rtB . blk35wfutl , & rtB . cy5jq00nod , & rtB . orhxqsmvfm , & rtB .
eee41maadj , & rtB . mcwb4pv1pd , & rtB . cjeggiimao , & rtB . du3hz2fyzs , &
rtB . o123mdetox , & rtB . iubu1mncgz , & rtB . ht3vr3fcw1 , & rtB .
idqoyikyv2 , & rtB . hgyc45hfyg , & rtB . nodqar1hqn , & rtB . oi30omicce , &
rtB . foabxlso54 , & rtB . bhydlzvdgu , & rtB . nhtjb14mqm , & rtB .
lphlstzov0 , & rtB . nc3ujqp1pg , & rtB . eqjfsu4kur , & rtB . pcceognrbb , &
rtB . hcfgbzdbcf , & rtB . ip00peawg3 , & rtB . na1kfpe1mt , & rtB .
jost0dydiu , & rtB . nmkcufy13s , & rtB . jnk41mcgi4 , & rtB . becosfp3hi , &
rtB . dcpeged5b2 , & rtB . mscp3mtqiu , & rtB . jsgd1kggkx , & rtB .
hpzll4vow0 , & rtB . cpcl24qdbb , & rtB . p2rr4uzpkw , & rtB . pebffqxdnd , &
rtB . kkbkkd5z5m , & rtB . fj3bxcs4ue , & rtB . mxnffxyshc , & rtB .
gztmpcdsdg , & rtB . g3vckbtrp3 , & rtB . fr14joswus , & rtB . dxeimajz1j , &
rtB . a5ddnr2k0h , & rtB . cqrrgwzxyv , & rtB . a1drd1cazf , & rtP .
Constant_Value , & rtP . Constant5_Value , & rtP . InitialSpeed_Value , & rtP
. Gain_Gain , & rtP . Gain3_Gain , & rtP . Gain5_Gain , & rtP .
Integrator10_IC , & rtP . Integrator11_IC , & rtP . Integrator12_IC , & rtP .
Integrator13_IC , & rtP . Integrator14_IC , & rtP . Integrator15_IC , & rtP .
Integrator16_IC , & rtP . Integrator4_IC , & rtP . Integrator5_IC , & rtP .
Integrator6_IC , & rtP . Integrator7_IC , & rtP . Integrator8_IC , & rtP .
Integrator9_IC , & rtP . Velocity_Amp , & rtP . Velocity_Bias , & rtP .
Velocity_Freq , & rtP . Velocity_Phase , & rtP . ManualSwitch_CurrentSetting
, & rtP . Delay_InitialCondition , & rtP . Saturation_LowerSat_mclkng5niv , &
rtP . TransferFcn_A_m2qs34lzon [ 0 ] , & rtP . TransferFcn_C_lu1wzht1on [ 0 ]
, & rtP . Saturation_LowerSat , & rtP . TransferFcn_A [ 0 ] , & rtP .
TransferFcn_C [ 0 ] , & rtP . Saturation_LowerSat_l0ksmwxpmq , & rtP .
TransferFcn_A_msxihde20w [ 0 ] , & rtP . TransferFcn_C_nnurr522qa [ 0 ] , &
rtP . Saturation_LowerSat_m0gx2ayxie , & rtP . TransferFcn_A_bwuzoqq3iu [ 0 ]
, & rtP . TransferFcn_C_j01bq5y01n [ 0 ] , & rtP .
Saturation_LowerSat_l5sm1vdo10 , & rtP . TransferFcn_A_cqt2tjfue3 [ 0 ] , &
rtP . TransferFcn_C_kavhf31tp1 [ 0 ] , & rtP . Relay_OnVal , & rtP .
Relay_OffVal , & rtP . Relay_YOn , & rtP . Relay_YOff , & rtP . Out1_Y0 , &
rtP . param , } ; static int32_T * rtVarDimsAddrMap [ ] = { ( NULL ) } ;
#endif
static TARGET_CONST rtwCAPI_DataTypeMap rtDataTypeMap [ ] = { { "double" ,
"real_T" , 0 , 0 , sizeof ( real_T ) , ( uint8_T ) SS_DOUBLE , 0 , 0 , 0 } ,
{ "unsigned char" , "uint8_T" , 0 , 0 , sizeof ( uint8_T ) , ( uint8_T )
SS_UINT8 , 0 , 0 , 0 } , { "struct" , "struct_BLEksg6c1ggfsVfFA9HVZ" , 15 , 1
, sizeof ( struct_BLEksg6c1ggfsVfFA9HVZ ) , ( uint8_T ) SS_STRUCT , 0 , 0 , 0
} } ;
#ifdef HOST_CAPI_BUILD
#undef sizeof
#endif
static TARGET_CONST rtwCAPI_ElementMap rtElementMap [ ] = { { ( NULL ) , 0 ,
0 , 0 , 0 } , { "P_H" , rt_offsetof ( struct_BLEksg6c1ggfsVfFA9HVZ , P_H ) ,
0 , 3 , 0 } , { "P_M" , rt_offsetof ( struct_BLEksg6c1ggfsVfFA9HVZ , P_M ) ,
0 , 3 , 0 } , { "P_L" , rt_offsetof ( struct_BLEksg6c1ggfsVfFA9HVZ , P_L ) ,
0 , 3 , 0 } , { "V1_0" , rt_offsetof ( struct_BLEksg6c1ggfsVfFA9HVZ , V1_0 )
, 0 , 3 , 0 } , { "Acap" , rt_offsetof ( struct_BLEksg6c1ggfsVfFA9HVZ , Acap
) , 0 , 3 , 0 } , { "V3_0" , rt_offsetof ( struct_BLEksg6c1ggfsVfFA9HVZ ,
V3_0 ) , 0 , 3 , 0 } , { "V4_0" , rt_offsetof ( struct_BLEksg6c1ggfsVfFA9HVZ
, V4_0 ) , 0 , 3 , 0 } , { "beta" , rt_offsetof (
struct_BLEksg6c1ggfsVfFA9HVZ , beta ) , 0 , 3 , 0 } , { "max_Avt" ,
rt_offsetof ( struct_BLEksg6c1ggfsVfFA9HVZ , max_Avt ) , 0 , 3 , 0 } , { "Cd"
, rt_offsetof ( struct_BLEksg6c1ggfsVfFA9HVZ , Cd ) , 0 , 3 , 0 } , { "wn" ,
rt_offsetof ( struct_BLEksg6c1ggfsVfFA9HVZ , wn ) , 0 , 3 , 0 } , { "zeta" ,
rt_offsetof ( struct_BLEksg6c1ggfsVfFA9HVZ , zeta ) , 0 , 3 , 0 } , { "J_hyd"
, rt_offsetof ( struct_BLEksg6c1ggfsVfFA9HVZ , J_hyd ) , 0 , 3 , 0 } , { "D"
, rt_offsetof ( struct_BLEksg6c1ggfsVfFA9HVZ , D ) , 0 , 3 , 0 } , { "J_elec"
, rt_offsetof ( struct_BLEksg6c1ggfsVfFA9HVZ , J_elec ) , 0 , 3 , 0 } } ;
static const rtwCAPI_DimensionMap rtDimensionMap [ ] = { { rtwCAPI_SCALAR , 0
, 2 , 0 } , { rtwCAPI_VECTOR , 2 , 2 , 0 } , { rtwCAPI_VECTOR , 4 , 2 , 0 } ,
{ rtwCAPI_MATRIX_COL_MAJOR , 0 , 2 , 0 } } ; static const uint_T
rtDimensionArray [ ] = { 1 , 1 , 2 , 1 , 1 , 2 } ; static const real_T
rtcapiStoredFloats [ ] = { 0.0 , 1.0 , 1.0E-6 } ; static const
rtwCAPI_FixPtMap rtFixPtMap [ ] = { { ( NULL ) , ( NULL ) ,
rtwCAPI_FIX_RESERVED , 0 , 0 , ( boolean_T ) 0 } , } ; static const
rtwCAPI_SampleTimeMap rtSampleTimeMap [ ] = { { ( NULL ) , ( NULL ) , - 1 , 0
} , { ( const void * ) & rtcapiStoredFloats [ 0 ] , ( const void * ) &
rtcapiStoredFloats [ 0 ] , ( int8_T ) 0 , ( uint8_T ) 0 } , { ( const void *
) & rtcapiStoredFloats [ 0 ] , ( const void * ) & rtcapiStoredFloats [ 1 ] ,
( int8_T ) 1 , ( uint8_T ) 0 } , { ( NULL ) , ( NULL ) , 3 , 0 } , { ( const
void * ) & rtcapiStoredFloats [ 2 ] , ( const void * ) & rtcapiStoredFloats [
0 ] , ( int8_T ) 2 , ( uint8_T ) 0 } } ; static
rtwCAPI_ModelMappingStaticInfo mmiStatic = { { rtBlockSignals , 96 ,
rtRootInputs , 0 , rtRootOutputs , 0 } , { rtBlockParameters , 45 ,
rtModelParameters , 1 } , { ( NULL ) , 0 } , { rtDataTypeMap , rtDimensionMap
, rtFixPtMap , rtElementMap , rtSampleTimeMap , rtDimensionArray } , "float"
, { 2806597705U , 165651654U , 2791374129U , 1143219654U } , ( NULL ) , 0 , (
boolean_T ) 0 , rt_LoggedStateIdxList } ; const
rtwCAPI_ModelMappingStaticInfo *
trial_with_intermediate_tanks_GetCAPIStaticMap ( void ) { return & mmiStatic
; }
#ifndef HOST_CAPI_BUILD
void trial_with_intermediate_tanks_InitializeDataMapInfo ( void ) {
rtwCAPI_SetVersion ( ( * rt_dataMapInfoPtr ) . mmi , 1 ) ;
rtwCAPI_SetStaticMap ( ( * rt_dataMapInfoPtr ) . mmi , & mmiStatic ) ;
rtwCAPI_SetLoggingStaticMap ( ( * rt_dataMapInfoPtr ) . mmi , ( NULL ) ) ;
rtwCAPI_SetDataAddressMap ( ( * rt_dataMapInfoPtr ) . mmi , rtDataAddrMap ) ;
rtwCAPI_SetVarDimsAddressMap ( ( * rt_dataMapInfoPtr ) . mmi ,
rtVarDimsAddrMap ) ; rtwCAPI_SetInstanceLoggingInfo ( ( * rt_dataMapInfoPtr )
. mmi , ( NULL ) ) ; rtwCAPI_SetChildMMIArray ( ( * rt_dataMapInfoPtr ) . mmi
, ( NULL ) ) ; rtwCAPI_SetChildMMIArrayLen ( ( * rt_dataMapInfoPtr ) . mmi ,
0 ) ; }
#else
#ifdef __cplusplus
extern "C" {
#endif
void trial_with_intermediate_tanks_host_InitializeDataMapInfo (
trial_with_intermediate_tanks_host_DataMapInfo_T * dataMap , const char *
path ) { rtwCAPI_SetVersion ( dataMap -> mmi , 1 ) ; rtwCAPI_SetStaticMap (
dataMap -> mmi , & mmiStatic ) ; rtwCAPI_SetDataAddressMap ( dataMap -> mmi ,
( NULL ) ) ; rtwCAPI_SetVarDimsAddressMap ( dataMap -> mmi , ( NULL ) ) ;
rtwCAPI_SetPath ( dataMap -> mmi , path ) ; rtwCAPI_SetFullPath ( dataMap ->
mmi , ( NULL ) ) ; rtwCAPI_SetChildMMIArray ( dataMap -> mmi , ( NULL ) ) ;
rtwCAPI_SetChildMMIArrayLen ( dataMap -> mmi , 0 ) ; }
#ifdef __cplusplus
}
#endif
#endif
