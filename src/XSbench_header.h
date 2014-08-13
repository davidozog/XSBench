#ifndef __XSBENCH_HEADER_H__
#define __XSBENCH_HEADER_H__

#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<string.h>
#include<strings.h>
#include<math.h>
#include<omp.h>
#include<unistd.h>
#include<sys/time.h>

// Papi Header
#ifdef PAPI
#include "papi.h"
#endif

// I/O Specifiers
#define INFO 1
#define DEBUG 1
#define SAVE 1

// Structures
typedef struct{
	double energy;
	double total_xs;
	double elastic_xs;
	double absorbtion_xs;
	double fission_xs;
	double nu_fission_xs;
} NuclideGridPoint;

typedef struct{
	double *energy;
	double *total_xs;
	double *elastic_xs;
	double *absorbtion_xs;
	double *fission_xs;
	double *nu_fission_xs;
} NuclideGridPoint_SOA;

typedef struct{
	double energy;
	int * xs_ptrs;
} GridPoint;

typedef struct{
	double *energy;
	int * xs_ptrs;
} GridPoint_SOA;

typedef struct{
	int nthreads;
	long n_isotopes;
	long n_gridpoints;
	int lookups;
	char * HM;
} Inputs;

// Function Prototypes
void logo(int version);
void center_print(const char *s, int width);
void border_print(void);
void fancy_int(long a);

NuclideGridPoint ** gpmatrix(size_t m, size_t n);
NuclideGridPoint_SOA * gpmatrix_SOA(size_t m, size_t n);

void gpmatrix_free( NuclideGridPoint ** M );
void gpmatrix_free_SOA( NuclideGridPoint_SOA * M );

int NGP_compare( const void * a, const void * b );
int NGP_compare_SOA( const void * a, const void * b );
int compare_double( const void *a, const void * b);

void generate_grids( NuclideGridPoint ** nuclide_grids,
                     long n_isotopes, long n_gridpoints );
void generate_grids_v( NuclideGridPoint ** nuclide_grids,
                     long n_isotopes, long n_gridpoints );
void generate_grids_v_SOA( NuclideGridPoint_SOA * nuclide_grids,
                     long n_isotopes, long n_gridpoints );

void sort_nuclide_grids( NuclideGridPoint ** nuclide_grids, long n_isotopes,
                         long n_gridpoints );
void sort_nuclide_grids_SOA( NuclideGridPoint_SOA * nuclide_grids, long n_isotopes,
                         long n_gridpoints );

GridPoint * generate_energy_grid( long n_isotopes, long n_gridpoints,
                                  NuclideGridPoint ** nuclide_grids);
GridPoint * generate_energy_grid_SOA( long n_isotopes, long n_gridpoints,
                                  NuclideGridPoint_SOA * nuclide_grids);

void set_grid_ptrs( GridPoint * energy_grid, NuclideGridPoint ** nuclide_grids,
                    long n_isotopes, long n_gridpoints );
void set_grid_ptrs_SOA( GridPoint * energy_grid, NuclideGridPoint_SOA * nuclide_grids,
                    long n_isotopes, long n_gridpoints );

int binary_search( NuclideGridPoint * A, double quarry, int n );
int binary_search_SOA( double * A, double quarry, int n );

void calculate_macro_xs(   double p_energy, int mat, long n_isotopes,
                           long n_gridpoints, int * restrict num_nucs,
                           double ** restrict concs,
						   GridPoint * restrict energy_grid,
                           NuclideGridPoint ** restrict nuclide_grids,
						   int ** restrict mats,
                           double * restrict macro_xs_vector );
void calculate_macro_xs_SOA(   double p_energy, int mat, long n_isotopes,
                           long n_gridpoints, int * restrict num_nucs,
                           double ** restrict concs,
						   GridPoint * restrict energy_grid,
						   GridPoint_SOA * restrict energy_grid_SOA,
                           NuclideGridPoint ** restrict nuclide_grids,
                           NuclideGridPoint_SOA * restrict nuclide_grids_SOA,
						   int ** restrict mats,
                           double * restrict macro_xs_vector );

void calculate_micro_xs(   double p_energy, int nuc, long n_isotopes,
                           long n_gridpoints,
                           GridPoint * restrict energy_grid,
                           NuclideGridPoint ** restrict nuclide_grids, int idx,
                           double * restrict xs_vector );
void calculate_micro_xs_SOA(   double p_energy, int nuc, long n_isotopes,
                           long n_gridpoints,
                           GridPoint * restrict energy_grid,
                           GridPoint_SOA * restrict energy_grid_SOA,
                           NuclideGridPoint ** restrict nuclide_grids, 
                           NuclideGridPoint_SOA * restrict nuclide_grids_SOA, 
                           int idx,
                           double * restrict xs_vector );

long grid_search( long n, double quarry, GridPoint * A);
long grid_search_SOA( long n, double quarry, GridPoint_SOA * A);

int * load_num_nucs(long n_isotopes);
int ** load_mats( int * num_nucs, long n_isotopes );
double ** load_concs( int * num_nucs );
double ** load_concs_v( int * num_nucs );
int pick_mat(unsigned long * seed);
double rn(unsigned long * seed);
int rn_int(unsigned long * seed);
void counter_stop( int * eventset, int num_papi_events );
void counter_init( int * eventset, int * num_papi_events );
void do_flops(void);
void do_loads( int nuc,
               NuclideGridPoint ** restrict nuclide_grids,
		       long n_gridpoints );	
Inputs read_CLI( int argc, char * argv[] );
void print_CLI_error(void);
double rn_v(void);
double round_double( double input );
unsigned int hash(unsigned char *str, int nbins);
size_t estimate_mem_usage( Inputs in );
void print_inputs(Inputs in, int nprocs, int version);
void print_results( Inputs in, int mype, double runtime, int nprocs, unsigned long long vhash );
void binary_dump(long n_isotopes, long n_gridpoints, NuclideGridPoint ** nuclide_grids, GridPoint * energy_grid);
void binary_read(long n_isotopes, long n_gridpoints, NuclideGridPoint ** nuclide_grids, GridPoint * energy_grid);
void qsort_SOA(double *base, 
               double *base1,
               double *base2,
               double *base3,
               double *base4,
               double *base5,
    unsigned num, unsigned width, int (*comp)(const void *, const void *));

void copy_AOS_to_SOA(GridPoint *energy_grid, GridPoint_SOA **energy_grid_SOA, 
           NuclideGridPoint **nuclide_grids, NuclideGridPoint_SOA *nuclide_grids_SOA, 
           int n_isotopes, int n_gridpoints);

#endif
