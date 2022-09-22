#pragma once

#include <math.h>

/**
 * Standard Atmosphere Computations: 
 * computes properties related to the 1976 standard atmosphere up to 84852 meters
 * 
 * input: 
 * \param height height in [m]
 * \param g Acceleration of gravity  in (m/s/s)
 * \param gamma Ratio of specific heats
 * \param R Gas constant for air in (J/kg-K)
 * 
 * output: 
 * \param temp temperature 
 * \param press pressure 
 * \param rho density 
 */

inline void IsoThermalAtmosphericLayer(double &T, double &P, double dh, double L, double g_R)
{ P *= exp(-dh*g_R/T); }

inline void GradientAtmosphericLayer(double &T, double &P, double dh, double L, double g_R)
{ auto ratio = 1 + L*dh/T; T *= ratio; P *= pow(ratio, -g_R/L); }

inline int getAtmosphericConditions(double &T, double &P, double &Rho, double h,
                                    double g=9.81, double R=287.058)
{
  // Altitudes (m)
  double H0=0.0, H1=11e3, H2=20e3, H3=32e3, H4=47e3, H5=51e3, H6=71e3, H7=84852;
  //Lapse Rates (K/m)
  double L1=-65e-4, L2=0, L3=1e-3, L4=28e-4, L5=0, L6=-28e-4, L7=-2e-3;
  // Ratio of g and R
  double g_R=g/R;
  // Initial Values
  T = 288.16; P = 1.01325e5;
  // Check if the calculation is out of bounds of the model's range of applicability
  if (h>H7) return 0;
  // Perform the calculation
  // !!! NESTED IF STATEMENTS -- INDENTATION IS WRONG BUT CODE IS READABLE !!!
  if (h>H0) {   GradientAtmosphericLayer(T, P, ((h>H1)?H1:h)-H0, L1, g_R);
  if (h>H1) { IsoThermalAtmosphericLayer(T, P, ((h>H2)?H2:h)-H1, L2, g_R);
  if (h>H2) {   GradientAtmosphericLayer(T, P, ((h>H3)?H3:h)-H2, L3, g_R);
  if (h>H3) {   GradientAtmosphericLayer(T, P, ((h>H4)?H4:h)-H3, L4, g_R);
  if (h>H4) { IsoThermalAtmosphericLayer(T, P, ((h>H5)?H5:h)-H4, L5, g_R);
  if (h>H5) {   GradientAtmosphericLayer(T, P, ((h>H6)?H6:h)-H5, L6, g_R);
  if (h>H6) {   GradientAtmosphericLayer(T, P, (          h)-H6, L7, g_R); }}}}}}}
  //Calculate density and return
  Rho = P/R/T; return 1;
}
