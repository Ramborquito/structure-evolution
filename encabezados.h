# include <cuda.h>
# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# define cudaCheckErrors(msg) \
    do { \
        cudaError_t __err = cudaGetLastError(); \
        if (__err != cudaSuccess) { \
            fprintf(stderr, "Fatal error: %s (%s at %s:%d)\n", \
                msg, cudaGetErrorString(__err), \
                __FILE__, __LINE__); \
            fprintf(stderr, "*** FAILED - ABORTING\n"); \
            exit(1); \
        } \
    } while (0)
# define PI 3.14159266f
# define CUTOFF_FAC 1.12246205f
# define PROB 0.2f
# define RAN rannew64
# define EPS_MIN 1.0e-6f
# define EPS 1.0f

typedef struct 
{
  float side;
  float dt;
  float mass_big;
  float mass_sml;
  float sigma_big;
  float sigma_sml;
  float cell_side_big;
  float cell_side_sml;
  float temp_set;
  float bin_size_gder;
  float range_gder;
  int ngrain_big;
  int ngrain_sml;
  int ngrain_tot;
  int ntags_big;
  int ntags_sml;
  int ncell_big;
  int ncell_med;
  int ncell_sml;
  int nrange_bb;
  int nrange_bs;
  int nrange_sb;
  int nrange_ss;
  int nbins_gder;
  int not_finished;
  int NH;
  int WCA_flag;
  long idum;
} 
parametros;

//====================================================================================

float RAN(long *);

void set_vec_int_hst(int *, int, int);
void set_vec_float_hst(float *, int, float);
void set_vec_float3_hst(float3 *, int, float3);

void cell_locate_hst(char, float3 *, int *, int *, parametros);

void update_verlet_init_hst(char, float3 *, float3 *, float3 *, float3 *, parametros);
void update_verlet_finish_hst(char, float3 *, float3 *, parametros);

void get_forces_same_hst(char, float3 *, float3 *, float *, float *, int *, int *, 
         parametros);
void get_forces_diff_hst(char, char, float3 *, float3 *, float3 *, float *, float *, 
         int *, int *, parametros);

void potential_wca_hst(float sigma, float dist_inv, float &potential, float &normal_force);

void set_temp(float3 *, float3 *, float [], parametros);

float calculate_temp_linear(float T0, float Tf, float t0, float tf, float t);

void get_gder_hst(char, char, float3 *, float3 *, float *, parametros);

// ===================================================================================

__device__ void potential_wca_dev(float sigma, float dist_inv, float &potential, float &normal_force);

__global__ void set_vec_int_dev(int *, int, int); 
__global__ void set_vec_float_dev(float *, int, float); 
__global__ void set_vec_float3_dev(float3 *, int, float3); 

__global__ void cell_locate_dev(char, float3 *, int *, int *, parametros);

__global__ void update_verlet_init_dev(char, float3 *, float3 *, float3 *, float3 *, 
                    parametros);
__global__ void update_verlet_finish_dev(char, float3 *, float3 *, parametros);

__global__ void get_forces_same_dev(char, float3 *, float3 *, float *, float *, 
                    int *, int *, parametros);
__global__ void get_forces_diff_dev(char, char, float3 *, float3 *, float3 *, 
                    float *, float *, int *, int *, parametros);

__global__ void get_gder_dev(char, char, float3 *, float3 *, float *, parametros);
