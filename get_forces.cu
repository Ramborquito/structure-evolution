#include "encabezados.h"

// ===================================================================================
// HOST
// ===================================================================================

void potential_wca_hst(float sigma, float dist_inv, float &potential, float &normal_force){
    float sor = dist_inv*sigma;
    float sor6, sor12;
    sor6 = sor*sor*sor;
    sor6 = sor6*sor6;
    sor12=sor6*sor6;

    potential += EPS * 4 * (sor12 - sor6) + 1;
    normal_force = - EPS * 4 * (12 * sor12 - 6 * sor6) * dist_inv;

}

void get_forces_same_hst(char type, float3 *rr_vec, float3 *ff_vec, float *vir_vec, 
         float *pot_vec, int *nocup_vec, int *cell_vec, parametros pars)
{
  static float3 rrm, rrn, drr, ff_pair, ffm;
  static float cutoff2, side, side_inv, cutoff, virial, potential, sigma,
                cell_side_inv, dist_inv, dist2, sor, sor6, normal_force;
  static int ii, jj, kk, i_del, j_del, k_del, iip, jjp, kkp, tag, ncell,
	     cell_index, mm, nn, tag_init, tag_end, ntags, nocup, ngrain, nrange;

  if (type == 'b')
  { 
    sigma = pars.sigma_big;
    ngrain = pars.ngrain_big;
    cell_side_inv = 1.0f/pars.cell_side_big;
    ntags = pars.ntags_big;
    ncell = pars.ncell_big;
    nrange = pars.nrange_bb;
  }
  else
  { 
    sigma = pars.sigma_sml;
    ngrain = pars.ngrain_sml;
    cell_side_inv = 1.0f/pars.cell_side_sml;
    ntags = pars.ntags_sml;
    ncell = pars.ncell_sml;
    nrange = pars.nrange_ss;
  }

    //wca
    cutoff = 1.122462048f * sigma;

    cutoff2 = cutoff*cutoff;
    side = pars.side;
    side_inv = 1.0f/side;


  for (mm = 0; mm < ngrain; mm++)
  {
    // fetch

    rrm = rr_vec[mm];
    ffm = ff_vec[mm];
    virial = vir_vec[mm];
    potential = pot_vec[mm];

    // calculate position of mm in the cells

    ii = (int) (cell_side_inv*rrm.x);
    jj = (int) (cell_side_inv*rrm.y);
    kk = (int) (cell_side_inv*rrm.z);
    if (ii == ncell) ii--;
    if (jj == ncell) jj--;
    if (kk == ncell) kk--;

    // run through neighbor cells

    for (i_del = -nrange; i_del <= nrange; i_del++) 
        for (j_del = -nrange; j_del <= nrange; j_del++)
            for (k_del = -nrange; k_del <= nrange; k_del++)
    {
      iip = (ii + i_del + ncell)%ncell;
      jjp = (jj + j_del + ncell)%ncell;
      kkp = (kk + k_del + ncell)%ncell;
      cell_index = iip + ncell*(jjp + ncell*kkp);
      nocup = nocup_vec[cell_index];
      tag_init = cell_index*ntags;
      tag_end = tag_init + nocup;
    
      // check grain in cell

      for (tag = tag_init; tag < tag_end; tag++) 
      {
        nn = cell_vec[tag];
        if (nn == mm) continue;
         
        // fetch another

        rrn = rr_vec[nn];

        // distance
  
        drr.x = rrn.x - rrm.x;
        drr.y = rrn.y - rrm.y;  
        drr.z = rrn.z - rrm.z;  

        // periodic boundary conditions

        drr.x -= side*floor(side_inv*drr.x + 0.5f);
        drr.y -= side*floor(side_inv*drr.y + 0.5f);
        drr.z -= side*floor(side_inv*drr.z + 0.5f);
  
        // distance and normal force
        
        dist2 = drr.x*drr.x + drr.y*drr.y + drr.z*drr.z;
        if (dist2 < cutoff2)
        {
            dist_inv = sqrt(1.0f/dist2);
            potential_wca_hst(sigma, dist_inv, potential, normal_force);
              
          // calculate normal force
              
          ff_pair.x = normal_force*dist_inv*drr.x;
          ff_pair.y = normal_force*dist_inv*drr.y;
          ff_pair.z = normal_force*dist_inv*drr.z;

          // sum forces and virial

          ffm.x += ff_pair.x;
          ffm.y += ff_pair.y;
          ffm.z += ff_pair.z;

          virial -= 0.5f*(drr.x*ff_pair.x + drr.y*ff_pair.y + drr.z*ff_pair.z);
        }
      }
    }
    
    // save 

    ff_vec[mm] = ffm;
    vir_vec[mm] = virial;
    pot_vec[mm] = potential;
  }
}

// ===================================================================================

void get_forces_diff_hst(char type, char type_other, float3 *rr_vec, float3 *ff_vec, 
         float3 *rr_other_vec, float *vir_vec, float *pot_vec, int *nocup_other_vec, 
         int *cell_other_vec, parametros pars)
{
  static float3 rrm, rrn, drr, ff_pair, ffm;
  static float cutoff2, side, side_inv, cutoff, sigma, virial, potential,
	       cell_side_inv, dist_inv, dist2, sor, sor6, normal_force;
  static int ii, jj, kk, i_del, j_del, k_del, iip, jjp, kkp, tag, ncell, mm, nn, 
             cell_index, tag_init, tag_end, ntags, nocup, ngrain, nrange;

  side = pars.side;
  side_inv = 1.0f/side;

  if (type == 'b')
  {
    ngrain = pars.ngrain_big;
    cell_side_inv = 1.0f/pars.cell_side_sml;
    sigma = 0.5f*(pars.sigma_big + pars.sigma_sml);
    ntags = pars.ntags_sml;
    ncell = pars.ncell_sml;
    nrange = pars.nrange_bs;
  }
  else
  {
    ngrain = pars.ngrain_sml;
    cell_side_inv = 1.0f/pars.cell_side_big;
    sigma = 0.5f*(pars.sigma_sml + pars.sigma_big);
    ntags = pars.ntags_big;
    ncell = pars.ncell_big;
    nrange = pars.nrange_sb;
  }

    //wca
    cutoff = 1.122462048f * sigma;

    cutoff2 = cutoff*cutoff;
    side = pars.side;
    side_inv = 1.0f/side;


  for (mm = 0; mm < ngrain; mm++)
  {
    // fetch

    rrm = rr_vec[mm];
    ffm = ff_vec[mm];
    virial = vir_vec[mm];
    potential = pot_vec[mm];

    // calculate cell for mm

    ii = (int) (cell_side_inv*rrm.x);
    jj = (int) (cell_side_inv*rrm.y);
    kk = (int) (cell_side_inv*rrm.z);
    if (ii == ncell) ii--;
    if (jj == ncell) jj--;
    if (kk == ncell) kk--;

    // run through neighbor cells

    for (i_del = -nrange; i_del <= nrange; i_del++) 
        for (j_del = -nrange; j_del <= nrange; j_del++)
            for (k_del = -nrange; k_del <= nrange; k_del++)
    {
      iip = (ii + i_del + ncell)%ncell;
      jjp = (jj + j_del + ncell)%ncell;
      kkp = (kk + k_del + ncell)%ncell;
      cell_index = iip + ncell*(jjp + ncell*kkp);
      nocup = nocup_other_vec[cell_index];
      tag_init = cell_index*ntags;
      tag_end = tag_init + nocup;
    
      // check others

      for (tag = tag_init; tag < tag_end; tag++) 
      {
        nn = cell_other_vec[tag];

        // fetch again

        rrn = rr_other_vec[nn];
        
        // coordinate differences
  
        drr.x = rrn.x - rrm.x;
        drr.y = rrn.y - rrm.y;  
        drr.z = rrn.z - rrm.z;  

        // periodic boundary conditions

        drr.x -= side*floor(side_inv*drr.x + 0.5f);	// per. bound. cond.
        drr.y -= side*floor(side_inv*drr.y + 0.5f);
        drr.z -= side*floor(side_inv*drr.z + 0.5f);
    
        // distance and normal force

        dist2 = drr.x*drr.x + drr.y*drr.y + drr.z*drr.z;
        if (dist2 < cutoff2)
        {
            dist_inv = sqrt(1.0f/dist2);
            potential_wca_hst(sigma, dist_inv, potential, normal_force);

          // calculate normal force
              
          ff_pair.x = normal_force*dist_inv*drr.x;
          ff_pair.y = normal_force*dist_inv*drr.y;
          ff_pair.z = normal_force*dist_inv*drr.z;

          // sum forces

          ffm.x += ff_pair.x;
          ffm.y += ff_pair.y;
          ffm.z += ff_pair.z;

          virial -= 0.5f*(drr.x*ff_pair.x + drr.y*ff_pair.y + drr.z*ff_pair.z);
        }
      }
    }

    // save

    ff_vec[mm] = ffm; 
    vir_vec[mm] = virial; 
    pot_vec[mm] = potential; 
  } 
  return;
}

// ===================================================================================
// DEVICE
// ===================================================================================

__device__ void potential_wca_dev(float sigma, float dist_inv, float &potential, float &normal_force) {
    float sor = dist_inv * sigma;
    float sor6, sor12;
    sor6 = sor * sor * sor;
    sor6 = sor6 * sor6;
    sor12 = sor6 * sor6;

    potential += EPS * 4 * (sor12 - sor6) + 1;
    normal_force = -EPS * 4 * (12 * sor12 - 6 * sor6) * dist_inv;
}

__global__ void get_forces_same_dev(char type, float3 *rr_vec, float3 *ff_vec, 
                    float *vir_vec, float *pot_vec, int *nocup_vec, 
                    int *cell_vec, parametros pars)
{
  float3 rrm, rrn, drr, ffm, ff_pair;
  float cutoff2, side, side_inv, cutoff, virial, potential, sigma, cell_side_inv,
         dist_inv, dist2, normal_force, temp;
  int ii, jj, kk, i_del, j_del, k_del, iip, jjp, kkp, tag, ncell, cell_index, mm, nn,
      tag_init, tag_end, ntags, nocup, ngrain, nrange;

  if (type == 'b')
  { 
    sigma = pars.sigma_big;
    ngrain = pars.ngrain_big;
    cell_side_inv = 1.0f/pars.cell_side_big;
    ntags = pars.ntags_big;
    ncell = pars.ncell_big;
    nrange = pars.nrange_bb;
    temp = pars.temp_set;
  }
  else
  { 
    sigma = pars.sigma_sml;
    ngrain = pars.ngrain_sml;
    cell_side_inv = 1.0f/pars.cell_side_sml;
    ntags = pars.ntags_sml;
    ncell = pars.ncell_sml;
    nrange = pars.nrange_ss;
    temp = pars.temp_set;
  }

    //wca
  cutoff = 1.122462048f * sigma;

  cutoff2 = cutoff*cutoff;
  side = pars.side;
  side_inv = 1.0f/side;


  mm = threadIdx.x + blockIdx.x*blockDim.x;

  if (mm < ngrain)
  {
    // fetch

    rrm = rr_vec[mm];
    ffm = ff_vec[mm];
    virial = vir_vec[mm];
    potential = pot_vec[mm];

    // calculate cell for mm

    ii = (int) (cell_side_inv*rrm.x);
    jj = (int) (cell_side_inv*rrm.y);
    kk = (int) (cell_side_inv*rrm.z);
    if (ii == ncell) ii--;
    if (jj == ncell) jj--;
    if (kk == ncell) kk--;

    // run through neighbor cells

    for (i_del = -nrange; i_del <= nrange; i_del++) 
        for (j_del = -nrange; j_del <= nrange; j_del++)
            for (k_del = -nrange; k_del <= nrange; k_del++)
    {
      iip = (ii + i_del + ncell)%ncell;
      jjp = (jj + j_del + ncell)%ncell;
      kkp = (kk + k_del + ncell)%ncell;
      cell_index = iip + ncell*(jjp + ncell*kkp);
      nocup = nocup_vec[cell_index];
      tag_init = cell_index*ntags;
      tag_end = tag_init + nocup;
    
      // check another grain

      for (tag = tag_init; tag < tag_end; tag++) 
      {
        nn = cell_vec[tag];
        if (nn == mm) continue;
         
        // another fetch

        rrn = rr_vec[nn];

        // coordinate differences
 
        drr.x = rrn.x - rrm.x;
        drr.y = rrn.y - rrm.y;  
        drr.z = rrn.z - rrm.z;  

        // periodic boundary conditions

        drr.x -= side*floor(side_inv*drr.x + 0.5f);	// per. bound. cond.
        drr.y -= side*floor(side_inv*drr.y + 0.5f);
        drr.z -= side*floor(side_inv*drr.z + 0.5f);
  
        // distance

        dist2 = drr.x*drr.x + drr.y*drr.y + drr.z*drr.z;
        if (dist2 < cutoff2)
        {
          dist_inv = sqrt(1.0f/dist2);
          potential_wca_dev(sigma, dist_inv, potential, normal_force);
              
          // calculate normal force
              
          ff_pair.x = normal_force*dist_inv*drr.x;
          ff_pair.y = normal_force*dist_inv*drr.y;
          ff_pair.z = normal_force*dist_inv*drr.z;

          // sum forces and virial

          ffm.x += ff_pair.x;
          ffm.y += ff_pair.y;
          ffm.z += ff_pair.z;

          virial -= 0.5f*(drr.x*ff_pair.x + drr.y*ff_pair.y + drr.z*ff_pair.z);
        }
      }
    }

    // save

    ff_vec[mm] = ffm;
    vir_vec[mm] = virial;
    pot_vec[mm] = potential;
  }
}

// ===================================================================================

__global__ void get_forces_diff_dev(char type, char type_other, float3 *rr_vec, 
                    float3 *ff_vec, float3 *rr_other_vec, float *vir_vec, 
                    float *pot_vec, int *nocup_other_vec, int *cell_other_vec,
                    parametros pars)
{
  float3 rrm, rrn, drr, ff_pair, ffm;
  float cutoff2, side, side_inv, cutoff, sigma, virial, potential, cell_side_inv,
	 dist_inv, dist2, normal_force, temp;
  int ii, jj, kk, i_del, j_del, k_del, iip, jjp, kkp, tag, ncell, mm, nn, cell_index,
      tag_init, tag_end, ntags, nocup, ngrain, nrange;

  if (type == 'b')
  {
    ngrain = pars.ngrain_big;
    cell_side_inv = 1.0f/pars.cell_side_sml;
    sigma = 0.5f*(pars.sigma_big + pars.sigma_sml);
    ntags = pars.ntags_sml;
    ncell = pars.ncell_sml;
    nrange = pars.nrange_bs;
    temp = pars.temp_set;

  }
  else
  {
    ngrain = pars.ngrain_sml;
    cell_side_inv = 1.0f/pars.cell_side_big;
    sigma = 0.5f*(pars.sigma_sml + pars.sigma_big);
    ntags = pars.ntags_big;
    ncell = pars.ncell_big;
    nrange = pars.nrange_sb;
    temp = pars.temp_set;

  }

  //wca
  cutoff = 1.122462048f * sigma;

  cutoff2 = cutoff*cutoff;
  side = pars.side;
  side_inv = 1.0f/side;

  // calcula interacciones sobre grandes

  mm = threadIdx.x + blockIdx.x*blockDim.x;

  if (mm < ngrain)
  {
    // fetch

    rrm = rr_vec[mm];
    ffm = ff_vec[mm];
    virial = vir_vec[mm];
    potential = pot_vec[mm];

    // calcula posicion de mm

    ii = (int) (cell_side_inv*rrm.x);
    jj = (int) (cell_side_inv*rrm.y);
    kk = (int) (cell_side_inv*rrm.z);
    if (ii == ncell) ii--;
    if (jj == ncell) jj--;
    if (kk == ncell) kk--;

    // recorre vecindario

    for (i_del = -nrange; i_del <= nrange; i_del++) 
        for (j_del = -nrange; j_del <= nrange; j_del++)
            for (k_del = -nrange; k_del <= nrange; k_del++)
    {
      iip = (ii + i_del + ncell)%ncell;
      jjp = (jj + j_del + ncell)%ncell;
      kkp = (kk + k_del + ncell)%ncell;
      cell_index = iip + ncell*(jjp + ncell*kkp);
      nocup = nocup_other_vec[cell_index];
      tag_init = cell_index*ntags;
      tag_end = tag_init + nocup;
    
      // checa con granos de taman~o diferente

      for (tag = tag_init; tag < tag_end; tag++) 
      {
        nn = cell_other_vec[tag];
        
        // another fetch

        rrn = rr_other_vec[nn];
        
        // corrdinate differences
  
        drr.x = rrn.x - rrm.x;
        drr.y = rrn.y - rrm.y;  
        drr.z = rrn.z - rrm.z;  

        // periodic boundary conditions

        drr.x -= side*floor(side_inv*drr.x + 0.5f);	// per. bound. cond.
        drr.y -= side*floor(side_inv*drr.y + 0.5f);
        drr.z -= side*floor(side_inv*drr.z + 0.5f);
    
        // distance

        dist2 = drr.x*drr.x + drr.y*drr.y + drr.z*drr.z ;
        if (dist2 < cutoff2)
        {
          dist_inv = sqrt(1.0f/dist2);
          potential_wca_dev(sigma, dist_inv, potential, normal_force);
              
          // calcula las fuerzas normales para este par
              
          ff_pair.x = normal_force*dist_inv*drr.x;
          ff_pair.y = normal_force*dist_inv*drr.y;
          ff_pair.z = normal_force*dist_inv*drr.z;

          // suma a fuerzas, al potencial y al virial

          ffm.x += ff_pair.x;
          ffm.y += ff_pair.y;
          ffm.z += ff_pair.z;

          virial -= 0.5f*(drr.x*ff_pair.x + drr.y*ff_pair.y + drr.z*ff_pair.z);
        }
      }
    }

    // save

    ff_vec[mm] = ffm;
    vir_vec[mm] = virial;
    pot_vec[mm] = potential;
  }
}
