/*
 * Copy_of_rotary_system_23a_data.c
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

/* Block parameters (default storage) */
P_Copy_of_rotary_system_23a_T Copy_of_rotary_system_23a_P = {
  /* Variable: param
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
  {
    0.001,
    2.0e-5,
    20000000.0,
    10000000.0,
    0.00015707963267948965,
    0.6,
    1800000000.0,
    870.0,
    0.0705,
    0.0705,
    1.4e-5,
    0.0003,
    1.6e-6
  },

  /* Expression: 0
   * Referenced by: '<Root>/Integrator2'
   */
  0.0,

  /* Expression: 60/(2*pi)
   * Referenced by: '<Root>/Gain2'
   */
  9.5492965855137211,

  /* Expression: 0
   * Referenced by: '<Root>/Constant3'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<Root>/Integrator3'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<Root>/Integrator4'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<Root>/Integrator5'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<Root>/Integrator6'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<Root>/Integrator8'
   */
  0.0,

  /* Computed Parameter: ManualSwitch_CurrentSetting
   * Referenced by: '<Root>/Manual Switch'
   */
  0U
};
