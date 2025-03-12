#include "encabezados.h"

// ================================================== 2023-05-22 =====================
// HOST			(VELOCITY VERLET)
// ===================================================================================

// initialize velocity verlet

void update_verlet_init_hst(char type, float3 *rr_vec, float3 *rr_raw_vec, 
         float3 *vv_vec, float3 *ff_vec, parametros pars)
{
  static float3 rr, rr_raw, drr, vv, ff;
  static float side, side_inv, dt, hdt_omass;
  static int ngrain, mm;

  // inicializa

  side = pars.side;
  dt = pars.dt;
  side_inv = 1.0/side;

  if (type == 'b') 
  {
    ngrain = pars.ngrain_big;
    hdt_omass = 0.5*dt/pars.mass_big;
  }
  else
  {
    ngrain = pars.ngrain_sml;
    hdt_omass = 0.5*dt/pars.mass_sml;
  }
  
  // run over all particles

  for (mm = 0; mm < ngrain; mm++)
  {
    // fetch

    rr = rr_vec[mm];
    rr_raw = rr_raw_vec[mm];
    vv = vv_vec[mm];
    ff = ff_vec[mm];

    // half-update velocities

    vv.x += hdt_omass*ff.x; 
    vv.y += hdt_omass*ff.y; 
    vv.z += hdt_omass*ff.z; 

    // move position

    drr.x = dt*vv.x;
    drr.y = dt*vv.y;
    drr.z = dt*vv.z;

    rr.x += drr.x;
    rr.y += drr.y;
    rr.z += drr.z;

    rr_raw.x += drr.x;
    rr_raw.y += drr.y;
    rr_raw.z += drr.z;

    // peridodic boundary conditions
  
    rr.x -= side*floor(side_inv*rr.x);
    rr.y -= side*floor(side_inv*rr.y);
    rr.z -= side*floor(side_inv*rr.z);

    // save

    rr_vec[mm] = rr;
    rr_raw_vec[mm] = rr_raw;
    vv_vec[mm] = vv;
  }

  return;
}

// ===================================================================================

void update_verlet_finish_hst(char type, float3 *vv_vec, float3 *ff_vec, 
         parametros pars)
{
  static float3 vv, ff;
  static float dt, hdt_omass;
  static int ngrain, mm;

  dt = pars.dt;

  if (type == 'b') 
  {
    ngrain = pars.ngrain_big;
    hdt_omass = 0.5*dt/pars.mass_big;
  }
  else
  {
    ngrain = pars.ngrain_sml;
    hdt_omass = 0.5*dt/pars.mass_sml;
  }

  // another half-update velocities

  for (mm = 0; mm < ngrain; mm++)
  {
    // fetch

    vv = vv_vec[mm];
    ff = ff_vec[mm];

    // increment velocity

    vv.x += hdt_omass*ff.x; 
    vv.y += hdt_omass*ff.y; 
    vv.z += hdt_omass*ff.z; 

    // save

    vv_vec[mm] = vv;
  }

  return;
}

// ===================================================================================
// DEVICE 		(VELOCITY VERLET)
// ===================================================================================

__global__ void update_verlet_init_dev(char type, float3 *rr_vec, 
	            float3 *rr_raw_vec, float3 *vv_vec, float3 *ff_vec, 
		    parametros pars)
{
  float3 rr, rr_raw, drr, vv, ff;
  float side, side_inv, dt, hdt_omass;
  int ngrain, mm;

  // inicializa

  side = pars.side;
  dt = pars.dt;
  side_inv = 1.0/side;

  if (type == 'b') 
  {
    ngrain = pars.ngrain_big;
    hdt_omass = 0.5*dt/pars.mass_big;
  }
  else
  {
    ngrain = pars.ngrain_sml;
    hdt_omass = 0.5*dt/pars.mass_sml;
  }
  
  mm = threadIdx.x + blockIdx.x*blockDim.x;

  if (mm < ngrain)
  {
    // fetch

    rr = rr_vec[mm];
    rr_raw = rr_raw_vec[mm];
    vv = vv_vec[mm];
    ff = ff_vec[mm];
    
    // half-update velocities

    vv.x += hdt_omass*ff.x; 
    vv.y += hdt_omass*ff.y; 
    vv.z += hdt_omass*ff.z; 

    // move position

    drr.x = dt*vv.x;
    drr.y = dt*vv.y;
    drr.z = dt*vv.z;

    rr.x += drr.x;
    rr.y += drr.y;
    rr.z += drr.z;

    rr_raw.x += drr.x;
    rr_raw.y += drr.y;
    rr_raw.z += drr.z;

    // periodic boundary conditions

    rr.x -= side*floor(side_inv*rr.x);
    rr.y -= side*floor(side_inv*rr.y);
    rr.z -= side*floor(side_inv*rr.z);
  
    // save

    rr_vec[mm] = rr;
    rr_raw_vec[mm] = rr_raw;
    vv_vec[mm] = vv;
  }

  return;
}

// ===================================================================================

__global__ void update_verlet_finish_dev(char type, float3 *vv_vec, float3 *ff_vec,
                    parametros pars)
{
  float3 vv, ff;
  float dt, hdt_omass;
  int ngrain, mm;

  dt = pars.dt;

  if (type == 'b') 
  {
    ngrain = pars.ngrain_big;
    hdt_omass = 0.5*dt/pars.mass_big;
  }
  else
  {
    ngrain = pars.ngrain_sml;
    hdt_omass = 0.5*dt/pars.mass_sml;
  }

  // another half-update velocities

  mm = threadIdx.x + blockIdx.x*blockDim.x;

  if (mm < ngrain)
  {
    // fetch

    vv = vv_vec[mm];
    ff = ff_vec[mm];

    // move velocities

    vv.x += hdt_omass*ff.x; 
    vv.y += hdt_omass*ff.y; 
    vv.z += hdt_omass*ff.z; 

    // save

    vv_vec[mm] = vv;
  }

  return;
}
