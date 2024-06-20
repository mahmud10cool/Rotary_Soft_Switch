/*
 * Copy_of_rotary_system_23a.h
 *
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Code generation for model "Copy_of_rotary_system_23a".
 *
 * Model version              : 1.8
 * Simulink Coder version : 9.9 (R2023a) 19-Nov-2022
 * C source code generated on : Thu Jun  6 14:54:00 2024
 *
 * Target selection: grt.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: Intel->x86-64 (Windows64)
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#ifndef RTW_HEADER_Copy_of_rotary_system_23a_h_
#define RTW_HEADER_Copy_of_rotary_system_23a_h_
#ifndef Copy_of_rotary_system_23a_COMMON_INCLUDES_
#define Copy_of_rotary_system_23a_COMMON_INCLUDES_
#include "rtwtypes.h"
#include "rtw_continuous.h"
#include "rtw_solver.h"
#include "rt_logging.h"
#endif                          /* Copy_of_rotary_system_23a_COMMON_INCLUDES_ */

#include "Copy_of_rotary_system_23a_types.h"
#include "rt_nonfinite.h"
#include "rtGetNaN.h"
#include <float.h>
#include <string.h>
#include <stddef.h>

/* Macros for accessing real-time model data structure */
#ifndef rtmGetContStateDisabled
#define rtmGetContStateDisabled(rtm)   ((rtm)->contStateDisabled)
#endif

#ifndef rtmSetContStateDisabled
#define rtmSetContStateDisabled(rtm, val) ((rtm)->contStateDisabled = (val))
#endif

#ifndef rtmGetContStates
#define rtmGetContStates(rtm)          ((rtm)->contStates)
#endif

#ifndef rtmSetContStates
#define rtmSetContStates(rtm, val)     ((rtm)->contStates = (val))
#endif

#ifndef rtmGetContTimeOutputInconsistentWithStateAtMajorStepFlag
#define rtmGetContTimeOutputInconsistentWithStateAtMajorStepFlag(rtm) ((rtm)->CTOutputIncnstWithState)
#endif

#ifndef rtmSetContTimeOutputInconsistentWithStateAtMajorStepFlag
#define rtmSetContTimeOutputInconsistentWithStateAtMajorStepFlag(rtm, val) ((rtm)->CTOutputIncnstWithState = (val))
#endif

#ifndef rtmGetDerivCacheNeedsReset
#define rtmGetDerivCacheNeedsReset(rtm) ((rtm)->derivCacheNeedsReset)
#endif

#ifndef rtmSetDerivCacheNeedsReset
#define rtmSetDerivCacheNeedsReset(rtm, val) ((rtm)->derivCacheNeedsReset = (val))
#endif

#ifndef rtmGetFinalTime
#define rtmGetFinalTime(rtm)           ((rtm)->Timing.tFinal)
#endif

#ifndef rtmGetIntgData
#define rtmGetIntgData(rtm)            ((rtm)->intgData)
#endif

#ifndef rtmSetIntgData
#define rtmSetIntgData(rtm, val)       ((rtm)->intgData = (val))
#endif

#ifndef rtmGetOdeF
#define rtmGetOdeF(rtm)                ((rtm)->odeF)
#endif

#ifndef rtmSetOdeF
#define rtmSetOdeF(rtm, val)           ((rtm)->odeF = (val))
#endif

#ifndef rtmGetOdeY
#define rtmGetOdeY(rtm)                ((rtm)->odeY)
#endif

#ifndef rtmSetOdeY
#define rtmSetOdeY(rtm, val)           ((rtm)->odeY = (val))
#endif

#ifndef rtmGetPeriodicContStateIndices
#define rtmGetPeriodicContStateIndices(rtm) ((rtm)->periodicContStateIndices)
#endif

#ifndef rtmSetPeriodicContStateIndices
#define rtmSetPeriodicContStateIndices(rtm, val) ((rtm)->periodicContStateIndices = (val))
#endif

#ifndef rtmGetPeriodicContStateRanges
#define rtmGetPeriodicContStateRanges(rtm) ((rtm)->periodicContStateRanges)
#endif

#ifndef rtmSetPeriodicContStateRanges
#define rtmSetPeriodicContStateRanges(rtm, val) ((rtm)->periodicContStateRanges = (val))
#endif

#ifndef rtmGetRTWLogInfo
#define rtmGetRTWLogInfo(rtm)          ((rtm)->rtwLogInfo)
#endif

#ifndef rtmGetZCCacheNeedsReset
#define rtmGetZCCacheNeedsReset(rtm)   ((rtm)->zCCacheNeedsReset)
#endif

#ifndef rtmSetZCCacheNeedsReset
#define rtmSetZCCacheNeedsReset(rtm, val) ((rtm)->zCCacheNeedsReset = (val))
#endif

#ifndef rtmGetdX
#define rtmGetdX(rtm)                  ((rtm)->derivs)
#endif

#ifndef rtmSetdX
#define rtmSetdX(rtm, val)             ((rtm)->derivs = (val))
#endif

#ifndef rtmGetErrorStatus
#define rtmGetErrorStatus(rtm)         ((rtm)->errorStatus)
#endif

#ifndef rtmSetErrorStatus
#define rtmSetErrorStatus(rtm, val)    ((rtm)->errorStatus = (val))
#endif

#ifndef rtmGetStopRequested
#define rtmGetStopRequested(rtm)       ((rtm)->Timing.stopRequestedFlag)
#endif

#ifndef rtmSetStopRequested
#define rtmSetStopRequested(rtm, val)  ((rtm)->Timing.stopRequestedFlag = (val))
#endif

#ifndef rtmGetStopRequestedPtr
#define rtmGetStopRequestedPtr(rtm)    (&((rtm)->Timing.stopRequestedFlag))
#endif

#ifndef rtmGetT
#define rtmGetT(rtm)                   (rtmGetTPtr((rtm))[0])
#endif

#ifndef rtmGetTFinal
#define rtmGetTFinal(rtm)              ((rtm)->Timing.tFinal)
#endif

#ifndef rtmGetTPtr
#define rtmGetTPtr(rtm)                ((rtm)->Timing.t)
#endif

/* Block signals (default storage) */
typedef struct {
  real_T Integrator1;                  /* '<Root>/Integrator1' */
  real_T Integrator2;                  /* '<Root>/Integrator2' */
  real_T RPM;                          /* '<Root>/Gain2' */
  real_T Integrator;                   /* '<Root>/Integrator' */
  real_T ManualSwitch;                 /* '<Root>/Manual Switch' */
  real_T KineticEnergy;                /* '<Root>/Gain1' */
  real_T Product1;                     /* '<Root>/Product1' */
  real_T Product3;                     /* '<Root>/Product3' */
  real_T Product4;                     /* '<Root>/Product4' */
  real_T Product5;                     /* '<Root>/Product5' */
  real_T Q_valve;                      /* '<Root>/Valve' */
  real_T P2dot;                        /* '<Root>/Small Chamber' */
  real_T shaft_acceleration;           /* '<Root>/Motor//Pump Dynamics' */
  real_T P1dot;                        /* '<Root>/Linear Actuator Chamber' */
  real_T Q_hyd;                        /* '<Root>/Hydraulic Pump//Motor Flow' */
} B_Copy_of_rotary_system_23a_T;

/* Continuous states (default storage) */
typedef struct {
  real_T Integrator1_CSTATE;           /* '<Root>/Integrator1' */
  real_T Integrator2_CSTATE;           /* '<Root>/Integrator2' */
  real_T Integrator_CSTATE;            /* '<Root>/Integrator' */
  real_T Integrator3_CSTATE;           /* '<Root>/Integrator3' */
  real_T Integrator4_CSTATE;           /* '<Root>/Integrator4' */
  real_T Integrator5_CSTATE;           /* '<Root>/Integrator5' */
  real_T Integrator6_CSTATE;           /* '<Root>/Integrator6' */
  real_T Integrator8_CSTATE;           /* '<Root>/Integrator8' */
} X_Copy_of_rotary_system_23a_T;

/* State derivatives (default storage) */
typedef struct {
  real_T Integrator1_CSTATE;           /* '<Root>/Integrator1' */
  real_T Integrator2_CSTATE;           /* '<Root>/Integrator2' */
  real_T Integrator_CSTATE;            /* '<Root>/Integrator' */
  real_T Integrator3_CSTATE;           /* '<Root>/Integrator3' */
  real_T Integrator4_CSTATE;           /* '<Root>/Integrator4' */
  real_T Integrator5_CSTATE;           /* '<Root>/Integrator5' */
  real_T Integrator6_CSTATE;           /* '<Root>/Integrator6' */
  real_T Integrator8_CSTATE;           /* '<Root>/Integrator8' */
} XDot_Copy_of_rotary_system_23_T;

/* State disabled  */
typedef struct {
  boolean_T Integrator1_CSTATE;        /* '<Root>/Integrator1' */
  boolean_T Integrator2_CSTATE;        /* '<Root>/Integrator2' */
  boolean_T Integrator_CSTATE;         /* '<Root>/Integrator' */
  boolean_T Integrator3_CSTATE;        /* '<Root>/Integrator3' */
  boolean_T Integrator4_CSTATE;        /* '<Root>/Integrator4' */
  boolean_T Integrator5_CSTATE;        /* '<Root>/Integrator5' */
  boolean_T Integrator6_CSTATE;        /* '<Root>/Integrator6' */
  boolean_T Integrator8_CSTATE;        /* '<Root>/Integrator8' */
} XDis_Copy_of_rotary_system_23_T;

#ifndef ODE4_INTG
#define ODE4_INTG

/* ODE4 Integration Data */
typedef struct {
  real_T *y;                           /* output */
  real_T *f[4];                        /* derivatives */
} ODE4_IntgData;

#endif

/* Parameters (default storage) */
struct P_Copy_of_rotary_system_23a_T_ {
  struct_5GcCYBJb8cnOhiSwmyHs9D param; /* Variable: param
                                        * Referenced by:
                                        *   '<Root>/Hydraulic Pump//Motor Flow'
                                        *   '<Root>/Linear Actuator Chamber'
                                        *   '<Root>/Motor//Pump Dynamics'
                                        *   '<Root>/Small Chamber'
                                        *   '<Root>/Valve'
                                        *   '<Root>/Constant'
                                        *   '<Root>/Constant1'
                                        *   '<Root>/Gain1'
                                        *   '<Root>/Integrator'
                                        *   '<Root>/Integrator1'
                                        */
  real_T Integrator2_IC;               /* Expression: 0
                                        * Referenced by: '<Root>/Integrator2'
                                        */
  real_T Gain2_Gain;                   /* Expression: 60/(2*pi)
                                        * Referenced by: '<Root>/Gain2'
                                        */
  real_T Constant3_Value;              /* Expression: 0
                                        * Referenced by: '<Root>/Constant3'
                                        */
  real_T Integrator3_IC;               /* Expression: 0
                                        * Referenced by: '<Root>/Integrator3'
                                        */
  real_T Integrator4_IC;               /* Expression: 0
                                        * Referenced by: '<Root>/Integrator4'
                                        */
  real_T Integrator5_IC;               /* Expression: 0
                                        * Referenced by: '<Root>/Integrator5'
                                        */
  real_T Integrator6_IC;               /* Expression: 0
                                        * Referenced by: '<Root>/Integrator6'
                                        */
  real_T Integrator8_IC;               /* Expression: 0
                                        * Referenced by: '<Root>/Integrator8'
                                        */
  uint8_T ManualSwitch_CurrentSetting;
                              /* Computed Parameter: ManualSwitch_CurrentSetting
                               * Referenced by: '<Root>/Manual Switch'
                               */
};

/* Real-time Model Data Structure */
struct tag_RTM_Copy_of_rotary_system_T {
  const char_T *errorStatus;
  RTWLogInfo *rtwLogInfo;
  RTWSolverInfo solverInfo;
  X_Copy_of_rotary_system_23a_T *contStates;
  int_T *periodicContStateIndices;
  real_T *periodicContStateRanges;
  real_T *derivs;
  XDis_Copy_of_rotary_system_23_T *contStateDisabled;
  boolean_T zCCacheNeedsReset;
  boolean_T derivCacheNeedsReset;
  boolean_T CTOutputIncnstWithState;
  real_T odeY[8];
  real_T odeF[4][8];
  ODE4_IntgData intgData;

  /*
   * Sizes:
   * The following substructure contains sizes information
   * for many of the model attributes such as inputs, outputs,
   * dwork, sample times, etc.
   */
  struct {
    int_T numContStates;
    int_T numPeriodicContStates;
    int_T numSampTimes;
  } Sizes;

  /*
   * Timing:
   * The following substructure contains information regarding
   * the timing information for the model.
   */
  struct {
    uint32_T clockTick0;
    uint32_T clockTickH0;
    time_T stepSize0;
    uint32_T clockTick1;
    uint32_T clockTickH1;
    time_T tFinal;
    SimTimeStep simTimeStep;
    boolean_T stopRequestedFlag;
    time_T *t;
    time_T tArray[2];
  } Timing;
};

/* Block parameters (default storage) */
extern P_Copy_of_rotary_system_23a_T Copy_of_rotary_system_23a_P;

/* Block signals (default storage) */
extern B_Copy_of_rotary_system_23a_T Copy_of_rotary_system_23a_B;

/* Continuous states (default storage) */
extern X_Copy_of_rotary_system_23a_T Copy_of_rotary_system_23a_X;

/* Model entry point functions */
extern void Copy_of_rotary_system_23a_initialize(void);
extern void Copy_of_rotary_system_23a_step(void);
extern void Copy_of_rotary_system_23a_terminate(void);

/* Real-time Model object */
extern RT_MODEL_Copy_of_rotary_syste_T *const Copy_of_rotary_system_23a_M;

/*-
 * The generated code includes comments that allow you to trace directly
 * back to the appropriate location in the model.  The basic format
 * is <system>/block_name, where system is the system number (uniquely
 * assigned by Simulink) and block_name is the name of the block.
 *
 * Use the MATLAB hilite_system command to trace the generated code back
 * to the model.  For example,
 *
 * hilite_system('<S3>')    - opens system 3
 * hilite_system('<S3>/Kp') - opens and selects block Kp which resides in S3
 *
 * Here is the system hierarchy for this model
 *
 * '<Root>' : 'Copy_of_rotary_system_23a'
 * '<S1>'   : 'Copy_of_rotary_system_23a/Hydraulic Pump//Motor Flow'
 * '<S2>'   : 'Copy_of_rotary_system_23a/Linear Actuator Chamber'
 * '<S3>'   : 'Copy_of_rotary_system_23a/Motor//Pump Dynamics'
 * '<S4>'   : 'Copy_of_rotary_system_23a/Small Chamber'
 * '<S5>'   : 'Copy_of_rotary_system_23a/Valve'
 */
#endif                             /* RTW_HEADER_Copy_of_rotary_system_23a_h_ */
