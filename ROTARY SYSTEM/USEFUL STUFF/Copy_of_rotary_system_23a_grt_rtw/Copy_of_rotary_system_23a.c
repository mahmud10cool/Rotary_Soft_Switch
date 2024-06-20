/*
 * Copy_of_rotary_system_23a.c
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

#include "Copy_of_rotary_system_23a.h"
#include "rtwtypes.h"
#include "rt_nonfinite.h"
#include <math.h>
#include "Copy_of_rotary_system_23a_private.h"

/* Block signals (default storage) */
B_Copy_of_rotary_system_23a_T Copy_of_rotary_system_23a_B;

/* Continuous states */
X_Copy_of_rotary_system_23a_T Copy_of_rotary_system_23a_X;

/* Real-time model */
static RT_MODEL_Copy_of_rotary_syste_T Copy_of_rotary_system_23a_M_;
RT_MODEL_Copy_of_rotary_syste_T *const Copy_of_rotary_system_23a_M =
  &Copy_of_rotary_system_23a_M_;

/*
 * This function updates continuous states using the ODE4 fixed-step
 * solver algorithm
 */
static void rt_ertODEUpdateContinuousStates(RTWSolverInfo *si )
{
  time_T t = rtsiGetT(si);
  time_T tnew = rtsiGetSolverStopTime(si);
  time_T h = rtsiGetStepSize(si);
  real_T *x = rtsiGetContStates(si);
  ODE4_IntgData *id = (ODE4_IntgData *)rtsiGetSolverData(si);
  real_T *y = id->y;
  real_T *f0 = id->f[0];
  real_T *f1 = id->f[1];
  real_T *f2 = id->f[2];
  real_T *f3 = id->f[3];
  real_T temp;
  int_T i;
  int_T nXc = 8;
  rtsiSetSimTimeStep(si,MINOR_TIME_STEP);

  /* Save the state values at time t in y, we'll use x as ynew. */
  (void) memcpy(y, x,
                (uint_T)nXc*sizeof(real_T));

  /* Assumes that rtsiSetT and ModelOutputs are up-to-date */
  /* f0 = f(t,y) */
  rtsiSetdX(si, f0);
  Copy_of_rotary_system_23a_derivatives();

  /* f1 = f(t + (h/2), y + (h/2)*f0) */
  temp = 0.5 * h;
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (temp*f0[i]);
  }

  rtsiSetT(si, t + temp);
  rtsiSetdX(si, f1);
  Copy_of_rotary_system_23a_step();
  Copy_of_rotary_system_23a_derivatives();

  /* f2 = f(t + (h/2), y + (h/2)*f1) */
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (temp*f1[i]);
  }

  rtsiSetdX(si, f2);
  Copy_of_rotary_system_23a_step();
  Copy_of_rotary_system_23a_derivatives();

  /* f3 = f(t + h, y + h*f2) */
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (h*f2[i]);
  }

  rtsiSetT(si, tnew);
  rtsiSetdX(si, f3);
  Copy_of_rotary_system_23a_step();
  Copy_of_rotary_system_23a_derivatives();

  /* tnew = t + h
     ynew = y + (h/6)*(f0 + 2*f1 + 2*f2 + 2*f3) */
  temp = h / 6.0;
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + temp*(f0[i] + 2.0*f1[i] + 2.0*f2[i] + f3[i]);
  }

  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

/* Model step function */
void Copy_of_rotary_system_23a_step(void)
{
  real_T KineticEnergy_tmp;
  real_T Q_hyd_tmp;
  real_T u_tmp;
  if (rtmIsMajorTimeStep(Copy_of_rotary_system_23a_M)) {
    /* set solver stop time */
    if (!(Copy_of_rotary_system_23a_M->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(&Copy_of_rotary_system_23a_M->solverInfo,
                            ((Copy_of_rotary_system_23a_M->Timing.clockTickH0 +
        1) * Copy_of_rotary_system_23a_M->Timing.stepSize0 * 4294967296.0));
    } else {
      rtsiSetSolverStopTime(&Copy_of_rotary_system_23a_M->solverInfo,
                            ((Copy_of_rotary_system_23a_M->Timing.clockTick0 + 1)
        * Copy_of_rotary_system_23a_M->Timing.stepSize0 +
        Copy_of_rotary_system_23a_M->Timing.clockTickH0 *
        Copy_of_rotary_system_23a_M->Timing.stepSize0 * 4294967296.0));
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(Copy_of_rotary_system_23a_M)) {
    Copy_of_rotary_system_23a_M->Timing.t[0] = rtsiGetT
      (&Copy_of_rotary_system_23a_M->solverInfo);
  }

  /* Integrator: '<Root>/Integrator1' */
  Copy_of_rotary_system_23a_B.Integrator1 =
    Copy_of_rotary_system_23a_X.Integrator1_CSTATE;
  if (rtmIsMajorTimeStep(Copy_of_rotary_system_23a_M)) {
  }

  /* Integrator: '<Root>/Integrator2' */
  Copy_of_rotary_system_23a_B.Integrator2 =
    Copy_of_rotary_system_23a_X.Integrator2_CSTATE;

  /* Gain: '<Root>/Gain2' */
  Copy_of_rotary_system_23a_B.RPM = Copy_of_rotary_system_23a_P.Gain2_Gain *
    Copy_of_rotary_system_23a_B.Integrator2;
  if (rtmIsMajorTimeStep(Copy_of_rotary_system_23a_M)) {
  }

  /* Integrator: '<Root>/Integrator' */
  Copy_of_rotary_system_23a_B.Integrator =
    Copy_of_rotary_system_23a_X.Integrator_CSTATE;
  if (rtmIsMajorTimeStep(Copy_of_rotary_system_23a_M)) {
    /* ManualSwitch: '<Root>/Manual Switch' */
    if (Copy_of_rotary_system_23a_P.ManualSwitch_CurrentSetting == 1) {
      /* ManualSwitch: '<Root>/Manual Switch' */
      Copy_of_rotary_system_23a_B.ManualSwitch = 0.0;
    } else {
      /* ManualSwitch: '<Root>/Manual Switch' incorporates:
       *  Constant: '<Root>/Constant'
       */
      Copy_of_rotary_system_23a_B.ManualSwitch =
        Copy_of_rotary_system_23a_P.param.P_H;
    }

    /* End of ManualSwitch: '<Root>/Manual Switch' */
  }

  /* MATLAB Function: '<Root>/Valve' incorporates:
   *  Constant: '<Root>/Constant1'
   *  Sum: '<Root>/Subtract'
   */
  u_tmp = Copy_of_rotary_system_23a_B.ManualSwitch -
    Copy_of_rotary_system_23a_B.Integrator;
  if (rtIsNaN(u_tmp)) {
    Q_hyd_tmp = (rtNaN);
  } else if (u_tmp < 0.0) {
    Q_hyd_tmp = -1.0;
  } else {
    Q_hyd_tmp = (u_tmp > 0.0);
  }

  Copy_of_rotary_system_23a_B.Q_valve = sqrt(2.0 /
    Copy_of_rotary_system_23a_P.param.rho * fabs(u_tmp)) *
    (Copy_of_rotary_system_23a_P.param.Cd *
     Copy_of_rotary_system_23a_P.param.max_Avt) * Q_hyd_tmp;

  /* End of MATLAB Function: '<Root>/Valve' */

  /* MATLAB Function: '<Root>/Hydraulic Pump//Motor Flow' */
  Q_hyd_tmp = Copy_of_rotary_system_23a_P.param.hyd_D / 6.2831853071795862;
  Copy_of_rotary_system_23a_B.Q_hyd = Q_hyd_tmp *
    Copy_of_rotary_system_23a_B.Integrator2;
  if (rtmIsMajorTimeStep(Copy_of_rotary_system_23a_M)) {
  }

  /* Gain: '<Root>/Gain1' incorporates:
   *  MATLAB Function: '<Root>/Motor//Pump Dynamics'
   */
  KineticEnergy_tmp = Copy_of_rotary_system_23a_P.param.J_elec +
    Copy_of_rotary_system_23a_P.param.J_hyd;

  /* Gain: '<Root>/Gain1' incorporates:
   *  Math: '<Root>/Square'
   */
  Copy_of_rotary_system_23a_B.KineticEnergy = KineticEnergy_tmp * 0.5 *
    (Copy_of_rotary_system_23a_B.Integrator2 *
     Copy_of_rotary_system_23a_B.Integrator2);
  if (rtmIsMajorTimeStep(Copy_of_rotary_system_23a_M)) {
  }

  /* MATLAB Function: '<Root>/Linear Actuator Chamber' */
  Copy_of_rotary_system_23a_B.P1dot = Copy_of_rotary_system_23a_P.param.beta /
    Copy_of_rotary_system_23a_P.param.V1_0 * Copy_of_rotary_system_23a_B.Q_hyd;

  /* MATLAB Function: '<Root>/Motor//Pump Dynamics' incorporates:
   *  Constant: '<Root>/Constant3'
   *  MATLAB Function: '<Root>/Hydraulic Pump//Motor Flow'
   */
  Copy_of_rotary_system_23a_B.shaft_acceleration =
    ((Copy_of_rotary_system_23a_B.Integrator -
      Copy_of_rotary_system_23a_B.Integrator1) * Q_hyd_tmp -
     Copy_of_rotary_system_23a_P.Constant3_Value) * (1.0 / KineticEnergy_tmp);

  /* Product: '<Root>/Product1' incorporates:
   *  Constant: '<Root>/Constant3'
   */
  Copy_of_rotary_system_23a_B.Product1 =
    Copy_of_rotary_system_23a_P.Constant3_Value *
    Copy_of_rotary_system_23a_B.Integrator2;

  /* Product: '<Root>/Product3' */
  Copy_of_rotary_system_23a_B.Product3 =
    Copy_of_rotary_system_23a_B.ManualSwitch *
    Copy_of_rotary_system_23a_B.Q_valve;

  /* Product: '<Root>/Product4' */
  Copy_of_rotary_system_23a_B.Product4 = Copy_of_rotary_system_23a_B.Integrator1
    * Copy_of_rotary_system_23a_B.Q_hyd;

  /* Product: '<Root>/Product5' */
  Copy_of_rotary_system_23a_B.Product5 = u_tmp *
    Copy_of_rotary_system_23a_B.Q_valve;

  /* MATLAB Function: '<Root>/Small Chamber' */
  Copy_of_rotary_system_23a_B.P2dot = Copy_of_rotary_system_23a_P.param.beta /
    Copy_of_rotary_system_23a_P.param.V2_0 *
    (Copy_of_rotary_system_23a_B.Q_valve - Copy_of_rotary_system_23a_B.Q_hyd);
  if (rtmIsMajorTimeStep(Copy_of_rotary_system_23a_M)) {
  }

  if (rtmIsMajorTimeStep(Copy_of_rotary_system_23a_M)) {
    /* Matfile logging */
    rt_UpdateTXYLogVars(Copy_of_rotary_system_23a_M->rtwLogInfo,
                        (Copy_of_rotary_system_23a_M->Timing.t));
  }                                    /* end MajorTimeStep */

  if (rtmIsMajorTimeStep(Copy_of_rotary_system_23a_M)) {
    /* signal main to stop simulation */
    {                                  /* Sample time: [0.0s, 0.0s] */
      if ((rtmGetTFinal(Copy_of_rotary_system_23a_M)!=-1) &&
          !((rtmGetTFinal(Copy_of_rotary_system_23a_M)-
             (((Copy_of_rotary_system_23a_M->Timing.clockTick1+
                Copy_of_rotary_system_23a_M->Timing.clockTickH1* 4294967296.0)) *
              1.0E-6)) > (((Copy_of_rotary_system_23a_M->Timing.clockTick1+
                            Copy_of_rotary_system_23a_M->Timing.clockTickH1*
                            4294967296.0)) * 1.0E-6) * (DBL_EPSILON))) {
        rtmSetErrorStatus(Copy_of_rotary_system_23a_M, "Simulation finished");
      }
    }

    rt_ertODEUpdateContinuousStates(&Copy_of_rotary_system_23a_M->solverInfo);

    /* Update absolute time for base rate */
    /* The "clockTick0" counts the number of times the code of this task has
     * been executed. The absolute time is the multiplication of "clockTick0"
     * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
     * overflow during the application lifespan selected.
     * Timer of this task consists of two 32 bit unsigned integers.
     * The two integers represent the low bits Timing.clockTick0 and the high bits
     * Timing.clockTickH0. When the low bit overflows to 0, the high bits increment.
     */
    if (!(++Copy_of_rotary_system_23a_M->Timing.clockTick0)) {
      ++Copy_of_rotary_system_23a_M->Timing.clockTickH0;
    }

    Copy_of_rotary_system_23a_M->Timing.t[0] = rtsiGetSolverStopTime
      (&Copy_of_rotary_system_23a_M->solverInfo);

    {
      /* Update absolute timer for sample time: [1.0E-6s, 0.0s] */
      /* The "clockTick1" counts the number of times the code of this task has
       * been executed. The resolution of this integer timer is 1.0E-6, which is the step size
       * of the task. Size of "clockTick1" ensures timer will not overflow during the
       * application lifespan selected.
       * Timer of this task consists of two 32 bit unsigned integers.
       * The two integers represent the low bits Timing.clockTick1 and the high bits
       * Timing.clockTickH1. When the low bit overflows to 0, the high bits increment.
       */
      Copy_of_rotary_system_23a_M->Timing.clockTick1++;
      if (!Copy_of_rotary_system_23a_M->Timing.clockTick1) {
        Copy_of_rotary_system_23a_M->Timing.clockTickH1++;
      }
    }
  }                                    /* end MajorTimeStep */
}

/* Derivatives for root system: '<Root>' */
void Copy_of_rotary_system_23a_derivatives(void)
{
  XDot_Copy_of_rotary_system_23_T *_rtXdot;
  _rtXdot = ((XDot_Copy_of_rotary_system_23_T *)
             Copy_of_rotary_system_23a_M->derivs);

  /* Derivatives for Integrator: '<Root>/Integrator1' */
  _rtXdot->Integrator1_CSTATE = Copy_of_rotary_system_23a_B.P1dot;

  /* Derivatives for Integrator: '<Root>/Integrator2' */
  _rtXdot->Integrator2_CSTATE = Copy_of_rotary_system_23a_B.shaft_acceleration;

  /* Derivatives for Integrator: '<Root>/Integrator' */
  _rtXdot->Integrator_CSTATE = Copy_of_rotary_system_23a_B.P2dot;

  /* Derivatives for Integrator: '<Root>/Integrator3' */
  _rtXdot->Integrator3_CSTATE = Copy_of_rotary_system_23a_B.Integrator2;

  /* Derivatives for Integrator: '<Root>/Integrator4' */
  _rtXdot->Integrator4_CSTATE = Copy_of_rotary_system_23a_B.Product5;

  /* Derivatives for Integrator: '<Root>/Integrator5' */
  _rtXdot->Integrator5_CSTATE = Copy_of_rotary_system_23a_B.Product1;

  /* Derivatives for Integrator: '<Root>/Integrator6' */
  _rtXdot->Integrator6_CSTATE = Copy_of_rotary_system_23a_B.Product3;

  /* Derivatives for Integrator: '<Root>/Integrator8' */
  _rtXdot->Integrator8_CSTATE = Copy_of_rotary_system_23a_B.Product4;
}

/* Model initialize function */
void Copy_of_rotary_system_23a_initialize(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* initialize real-time model */
  (void) memset((void *)Copy_of_rotary_system_23a_M, 0,
                sizeof(RT_MODEL_Copy_of_rotary_syste_T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&Copy_of_rotary_system_23a_M->solverInfo,
                          &Copy_of_rotary_system_23a_M->Timing.simTimeStep);
    rtsiSetTPtr(&Copy_of_rotary_system_23a_M->solverInfo, &rtmGetTPtr
                (Copy_of_rotary_system_23a_M));
    rtsiSetStepSizePtr(&Copy_of_rotary_system_23a_M->solverInfo,
                       &Copy_of_rotary_system_23a_M->Timing.stepSize0);
    rtsiSetdXPtr(&Copy_of_rotary_system_23a_M->solverInfo,
                 &Copy_of_rotary_system_23a_M->derivs);
    rtsiSetContStatesPtr(&Copy_of_rotary_system_23a_M->solverInfo, (real_T **)
                         &Copy_of_rotary_system_23a_M->contStates);
    rtsiSetNumContStatesPtr(&Copy_of_rotary_system_23a_M->solverInfo,
      &Copy_of_rotary_system_23a_M->Sizes.numContStates);
    rtsiSetNumPeriodicContStatesPtr(&Copy_of_rotary_system_23a_M->solverInfo,
      &Copy_of_rotary_system_23a_M->Sizes.numPeriodicContStates);
    rtsiSetPeriodicContStateIndicesPtr(&Copy_of_rotary_system_23a_M->solverInfo,
      &Copy_of_rotary_system_23a_M->periodicContStateIndices);
    rtsiSetPeriodicContStateRangesPtr(&Copy_of_rotary_system_23a_M->solverInfo,
      &Copy_of_rotary_system_23a_M->periodicContStateRanges);
    rtsiSetErrorStatusPtr(&Copy_of_rotary_system_23a_M->solverInfo,
                          (&rtmGetErrorStatus(Copy_of_rotary_system_23a_M)));
    rtsiSetRTModelPtr(&Copy_of_rotary_system_23a_M->solverInfo,
                      Copy_of_rotary_system_23a_M);
  }

  rtsiSetSimTimeStep(&Copy_of_rotary_system_23a_M->solverInfo, MAJOR_TIME_STEP);
  Copy_of_rotary_system_23a_M->intgData.y = Copy_of_rotary_system_23a_M->odeY;
  Copy_of_rotary_system_23a_M->intgData.f[0] = Copy_of_rotary_system_23a_M->
    odeF[0];
  Copy_of_rotary_system_23a_M->intgData.f[1] = Copy_of_rotary_system_23a_M->
    odeF[1];
  Copy_of_rotary_system_23a_M->intgData.f[2] = Copy_of_rotary_system_23a_M->
    odeF[2];
  Copy_of_rotary_system_23a_M->intgData.f[3] = Copy_of_rotary_system_23a_M->
    odeF[3];
  Copy_of_rotary_system_23a_M->contStates = ((X_Copy_of_rotary_system_23a_T *)
    &Copy_of_rotary_system_23a_X);
  rtsiSetSolverData(&Copy_of_rotary_system_23a_M->solverInfo, (void *)
                    &Copy_of_rotary_system_23a_M->intgData);
  rtsiSetIsMinorTimeStepWithModeChange(&Copy_of_rotary_system_23a_M->solverInfo,
    false);
  rtsiSetSolverName(&Copy_of_rotary_system_23a_M->solverInfo,"ode4");
  rtmSetTPtr(Copy_of_rotary_system_23a_M,
             &Copy_of_rotary_system_23a_M->Timing.tArray[0]);
  rtmSetTFinal(Copy_of_rotary_system_23a_M, 1.0);
  Copy_of_rotary_system_23a_M->Timing.stepSize0 = 1.0E-6;

  /* Setup for data logging */
  {
    static RTWLogInfo rt_DataLoggingInfo;
    rt_DataLoggingInfo.loggingInterval = (NULL);
    Copy_of_rotary_system_23a_M->rtwLogInfo = &rt_DataLoggingInfo;
  }

  /* Setup for data logging */
  {
    rtliSetLogXSignalInfo(Copy_of_rotary_system_23a_M->rtwLogInfo, (NULL));
    rtliSetLogXSignalPtrs(Copy_of_rotary_system_23a_M->rtwLogInfo, (NULL));
    rtliSetLogT(Copy_of_rotary_system_23a_M->rtwLogInfo, "tout");
    rtliSetLogX(Copy_of_rotary_system_23a_M->rtwLogInfo, "");
    rtliSetLogXFinal(Copy_of_rotary_system_23a_M->rtwLogInfo, "");
    rtliSetLogVarNameModifier(Copy_of_rotary_system_23a_M->rtwLogInfo, "rt_");
    rtliSetLogFormat(Copy_of_rotary_system_23a_M->rtwLogInfo, 4);
    rtliSetLogMaxRows(Copy_of_rotary_system_23a_M->rtwLogInfo, 0);
    rtliSetLogDecimation(Copy_of_rotary_system_23a_M->rtwLogInfo, 1);
    rtliSetLogY(Copy_of_rotary_system_23a_M->rtwLogInfo, "");
    rtliSetLogYSignalInfo(Copy_of_rotary_system_23a_M->rtwLogInfo, (NULL));
    rtliSetLogYSignalPtrs(Copy_of_rotary_system_23a_M->rtwLogInfo, (NULL));
  }

  /* block I/O */
  (void) memset(((void *) &Copy_of_rotary_system_23a_B), 0,
                sizeof(B_Copy_of_rotary_system_23a_T));

  /* states (continuous) */
  {
    (void) memset((void *)&Copy_of_rotary_system_23a_X, 0,
                  sizeof(X_Copy_of_rotary_system_23a_T));
  }

  /* Matfile logging */
  rt_StartDataLoggingWithStartTime(Copy_of_rotary_system_23a_M->rtwLogInfo, 0.0,
    rtmGetTFinal(Copy_of_rotary_system_23a_M),
    Copy_of_rotary_system_23a_M->Timing.stepSize0, (&rtmGetErrorStatus
    (Copy_of_rotary_system_23a_M)));

  /* InitializeConditions for Integrator: '<Root>/Integrator1' */
  Copy_of_rotary_system_23a_X.Integrator1_CSTATE =
    Copy_of_rotary_system_23a_P.param.P_M;

  /* InitializeConditions for Integrator: '<Root>/Integrator2' */
  Copy_of_rotary_system_23a_X.Integrator2_CSTATE =
    Copy_of_rotary_system_23a_P.Integrator2_IC;

  /* InitializeConditions for Integrator: '<Root>/Integrator' incorporates:
   *  Integrator: '<Root>/Integrator1'
   */
  Copy_of_rotary_system_23a_X.Integrator_CSTATE =
    Copy_of_rotary_system_23a_P.param.P_M;

  /* InitializeConditions for Integrator: '<Root>/Integrator3' */
  Copy_of_rotary_system_23a_X.Integrator3_CSTATE =
    Copy_of_rotary_system_23a_P.Integrator3_IC;

  /* InitializeConditions for Integrator: '<Root>/Integrator4' */
  Copy_of_rotary_system_23a_X.Integrator4_CSTATE =
    Copy_of_rotary_system_23a_P.Integrator4_IC;

  /* InitializeConditions for Integrator: '<Root>/Integrator5' */
  Copy_of_rotary_system_23a_X.Integrator5_CSTATE =
    Copy_of_rotary_system_23a_P.Integrator5_IC;

  /* InitializeConditions for Integrator: '<Root>/Integrator6' */
  Copy_of_rotary_system_23a_X.Integrator6_CSTATE =
    Copy_of_rotary_system_23a_P.Integrator6_IC;

  /* InitializeConditions for Integrator: '<Root>/Integrator8' */
  Copy_of_rotary_system_23a_X.Integrator8_CSTATE =
    Copy_of_rotary_system_23a_P.Integrator8_IC;
}

/* Model terminate function */
void Copy_of_rotary_system_23a_terminate(void)
{
  /* (no terminate code required) */
}
