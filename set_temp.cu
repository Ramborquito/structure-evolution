#include "encabezados.h"

// ===================================================================================
//   ONLY HOST
// ===================================================================================

float calculate_temp_linear(float T0, float Tf, float t0, float tf, float t){
    if(t < t0) return T0;
    if(t > tf) return Tf;
    float m = (Tf-T0)/(tf-t0);
    return T0 + m * (t - t0);
}

float calculate_temp_sine(float T0, float Tf, float t0, float tf, float period, float t){
    // a + b sin [2pi ( t / tau - 1 / 4 )]
    float a, b, argument;
    if(t < t0) return T0;
    if(t > tf) return Tf;
    a = 0.5 * (T0 + Tf);
    b = 0.5 * (Tf - T0);
    argument = 2 * PI * (t / period - 0.25);

    return a + b * sin(argument);
}


void set_temp(float3 *vv_big_vec, float3 *vv_sml_vec, float energs[2],
         parametros pars)
{
  static float3 vv;
  static int mm, ngrain_big, ngrain_sml;
  static float ene_kin_big, ene_kin_sml, mass_big, mass_sml, temp_set,
                temp_present_big, temp_present_sml, temp_big, temp_sml, coeff, ran1,
                ran2, gauss1, gauss2, rr2, aux, xngrain_big, xngrain_sml;
  static int not_init = 1;
  static long idum;

  if (not_init)
  {
    ngrain_big = pars.ngrain_big;
    ngrain_sml = pars.ngrain_sml;
    mass_big = pars.mass_big;
    mass_sml = pars.mass_sml;
    idum = pars.idum;
    xngrain_big = (float) ngrain_big;
    xngrain_sml = (float) ngrain_sml;
    not_init = 0;
  }

  // obtain kinetic energies and instantaneous temperature
  temp_set = pars.temp_set;
  ene_kin_big = 0.0;
  for (mm = 0; mm < ngrain_big; mm++)
  {
    vv = vv_big_vec[mm];
    ene_kin_big += 0.5*mass_big*(vv.x*vv.x + vv.y*vv.y + vv.z*vv.z);
  }
  ene_kin_big /= xngrain_big;
  energs[0] = ene_kin_big;
  temp_present_big = (2.0/3.0)*ene_kin_big;

  ene_kin_sml = 0.0;
  for (mm = 0; mm < ngrain_sml; mm++)
  {
    vv = vv_sml_vec[mm];
    ene_kin_sml += 0.5*mass_sml*(vv.x*vv.x + vv.y*vv.y + vv.z*vv.z);
  }
  if (ngrain_sml) ene_kin_sml /= xngrain_sml;
  energs[1] = ene_kin_sml;
  temp_present_sml = (2.0/3.0)*ene_kin_sml;

  // get gaussian random for new temps assignment

  while (1)
  {
    while (1)
    {
      ran1 = 2.0*RAN(&idum) - 1.0;
      ran2 = 2.0*RAN(&idum) - 1.0;
      rr2 = ran1*ran1 + ran2*ran2;
      if (rr2 > 0.0 && rr2 < 1.0)
      {
        aux = sqrt(-2.0*log(rr2)/rr2);
        gauss1 = ran1*aux;
        gauss2 = ran2*aux;
        break;
      }
    }
    temp_big = temp_set*(1.0 + gauss1/sqrt(xngrain_big));
    if (temp_big >= 0.0) break;
  }

  while (1)
  {
    while (1)
    {
      ran1 = 2.0*RAN(&idum) - 1.0;
      ran2 = 2.0*RAN(&idum) - 1.0;
      rr2 = ran1*ran1 + ran2*ran2;
      if (rr2 > 0.0 && rr2 < 1.0)
      {
        aux = sqrt(-2.0*log(rr2)/rr2);
        gauss1 = ran1*aux;
        gauss2 = ran2*aux;
        break;
      }
    }
    temp_sml = temp_set*(1.0 + gauss2/sqrt(xngrain_sml));
    if (temp_sml >= 0.0) break;
  }

  // now rescales velocities for grains
      
  coeff = sqrt(temp_big/temp_present_big);
  for (mm = 0; mm < ngrain_big; mm++)
  {
    vv_big_vec[mm].x *= coeff; 
    vv_big_vec[mm].y *= coeff; 
    vv_big_vec[mm].z *= coeff; 
  }
  coeff = sqrt(temp_sml/temp_present_sml);
  for (mm = 0; mm < ngrain_sml; mm++)
  {
    vv_sml_vec[mm].x *= coeff; 
    vv_sml_vec[mm].y *= coeff; 
    vv_sml_vec[mm].z *= coeff; 
  }

  return;
}
