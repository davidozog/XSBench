#include "XSbench_header.h"
#include <immintrin.h>


// Calculates the microscopic cross section for a given nuclide & energy
void calculate_micro_xs(   double p_energy, int nuc, int prefetch_nuc, int prefetch_nuc0,
                           long n_isotopes, long n_gridpoints,
                           GridPoint * restrict energy_grid,
                           NuclideGridPoint ** restrict nuclide_grids,
                           int idx, double * restrict xs_vector ){
	
	// Variables
	double f;
	NuclideGridPoint * low, * high;

//#pragma noprefetch nuclide_grids
  //_mm_prefetch((const char *)&nuclide_grids[prefetch_nuc][energy_grid[idx].xs_ptrs[prefetch_nuc]-1], _MM_HINT_T1); // vprefetch1
  //_mm_prefetch((const char *)&nuclide_grids[prefetch_nuc0][energy_grid[idx].xs_ptrs[prefetch_nuc0]-1], _MM_HINT_T0); // vprefetch0
	// pull ptr from energy grid and check to ensure that
	// we're not reading off the end of the nuclide's grid
	if( energy_grid[idx].xs_ptrs[nuc] == n_gridpoints - 1 )
		low = &nuclide_grids[nuc][energy_grid[idx].xs_ptrs[nuc] - 1];
	else
		low = &nuclide_grids[nuc][energy_grid[idx].xs_ptrs[nuc]];
	
	high = low + 1;
	
//#pragma noprefetch high, low
//Not sure why this seems to help:
//#pragma prefetch high:1:30
//#pragma prefetch high:0:5
//#pragma prefetch low:1:30
//#pragma prefetch low:0:5
	// calculate the re-useable interpolation factor
	f = (high->energy - p_energy) / (high->energy - low->energy);

	// Total XS
	xs_vector[0] = high->total_xs - f * (high->total_xs - low->total_xs);
	
	// Elastic XS
	xs_vector[1] = high->elastic_xs - f * (high->elastic_xs - low->elastic_xs);
	
	// Absorbtion XS
	xs_vector[2] = high->absorbtion_xs - f * (high->absorbtion_xs - low->absorbtion_xs);
	
	// Fission XS
	xs_vector[3] = high->fission_xs - f * (high->fission_xs - low->fission_xs);
	
	// Nu Fission XS
	xs_vector[4] = high->nu_fission_xs - f * (high->nu_fission_xs - low->nu_fission_xs);
	
	//test
	/*	
	if( omp_get_thread_num() == 0 )
	{
		printf("Lookup: Energy = %lf, nuc = %d\n", p_energy, nuc);
		printf("e_h = %lf e_l = %lf\n", high->energy , low->energy);
		printf("xs_h = %lf xs_l = %lf\n", high->elastic_xs, low->elastic_xs);
		printf("total_xs = %lf\n\n", xs_vector[1]);
	}
	*/
	
}

// Calculates macroscopic cross section based on a given material & energy 
void calculate_macro_xs( double p_energy, int mat, long n_isotopes,
                         long n_gridpoints, int * restrict num_nucs,
                         double ** restrict concs,
                         GridPoint * restrict energy_grid,
                         NuclideGridPoint ** restrict nuclide_grids,
                         int ** restrict mats,
                         double * restrict macro_xs_vector ){
	__attribute__((align(64))) double xs_vector[5];
	int nucs, p_nuc, prefetch_nuc, prefetch_nuc0; // the nuclide we are looking up
	long idx = 0;	
	double conc; // the concentration of the nuclide in the material

	// cleans out macro_xs_vector
	for( int k = 0; k < 5; k++ ) {
		macro_xs_vector[k] = 0;
  }

	// binary search for energy on unionized energy grid (UEG)
	idx = grid_search( n_isotopes * n_gridpoints, p_energy,
	                   energy_grid);	

  nucs = num_nucs[mat];
  _mm_prefetch((const char *)&energy_grid[idx], _MM_HINT_T1); // vprefetch1
	
	// Once we find the pointer array on the UEG, we can pull the data
	// from the respective nuclide grids, as well as the nuclide
	// concentration data for the material
	// Each nuclide from the material needs to have its micro-XS array
	// looked up & interpolatied (via calculate_micro_xs). Then, the
	// micro XS is multiplied by the concentration of that nuclide
	// in the material, and added to the total macro XS array.
#pragma novector
	for( int j = 0; j < nucs; j++ )
	{
    //printf("%d, ", j);
    _mm_prefetch((const char *)&mats[mat][j+1], _MM_HINT_T1); // vprefetch1
    _mm_prefetch((const char *)&mats[mat][j+1], _MM_HINT_T0); // vprefetch0
    _mm_prefetch((const char *)&concs[mat][j+1], _MM_HINT_T1); // vprefetch1
    _mm_prefetch((const char *)&concs[mat][j+1], _MM_HINT_T0); // vprefetch0
		p_nuc = mats[mat][j];
		conc = concs[mat][j];
    //printf("%d, ", p_nuc);
    //int *myp = &p_nuc+1;
    //prefetch_nuc = *myp;
//    _mm_prefetch((const char *)&energy_grid[idx].xs_ptrs[prefetch_nuc]-1, _MM_HINT_T1); // vprefetch1
//    _mm_prefetch((const char *)&energy_grid[idx].xs_ptrs[prefetch_nuc]-1, _MM_HINT_T0); // vprefetch0
    //if (j+100 < nucs)	
    //prefetch_nuc = mats[mat][j+1];
    //else
    //  prefetch_nuc = mats[mat][j+1];
    //if (j+5 < nucs)	
    prefetch_nuc0 = j+1;
    //else
    //  prefetch_nuc0 = mats[mat][j+1];
		calculate_micro_xs( p_energy, p_nuc, prefetch_nuc, prefetch_nuc0, n_isotopes,
		                    n_gridpoints, energy_grid,
		                    nuclide_grids, idx, xs_vector );

			macro_xs_vector[0] += xs_vector[0] * conc;
			macro_xs_vector[1] += xs_vector[1] * conc;
			macro_xs_vector[2] += xs_vector[2] * conc;
			macro_xs_vector[3] += xs_vector[3] * conc;
			macro_xs_vector[4] += xs_vector[4] * conc;
	}

//	for( int j = 0; j < num_nucs[mat]; j++ )
//  {
//		conc = concs[mat][j];
//		macro_xs_vector[0] = xs_vector[0] * conc;
//		macro_xs_vector[1] = xs_vector[1] * conc;
//		macro_xs_vector[2] = xs_vector[2] * conc;
//		macro_xs_vector[3] = xs_vector[3] * conc;
//		macro_xs_vector[4] = xs_vector[4] * conc;
////		for( int k = 0; k < 5; k++ ) {
////     // int slot = k*64;
////			macro_xs_vector[k] += xs_vector[k] * conc;
////    }
//	}

//  #pragma simd
//	for( int j = 0; j < num_nucs[mat]; j++ )
//	{
//		p_nuc = mats[mat][j];
//		conc = concs[mat][j];
//		calculate_micro_xs( p_energy, p_nuc, n_isotopes,
//		                    n_gridpoints, energy_grid,
//		                    nuclide_grids, idx, xs_vector );
//		//for( int k = 0; k < 5; k++ )
//			macro_xs_vector[0] += xs_vector[0] * conc;
//			macro_xs_vector[1] += xs_vector[1] * conc;
//			macro_xs_vector[2] += xs_vector[2] * conc;
//			macro_xs_vector[3] += xs_vector[3] * conc;
//			macro_xs_vector[4] += xs_vector[4] * conc;
//	}
	
	//test
	/*
	for( int k = 0; k < 5; k++ )
		printf("Energy: %lf, Material: %d, XSVector[%d]: %lf\n",
		       p_energy, mat, k, macro_xs_vector[k]);
	*/
}


// (fixed) binary search for energy on unionized energy grid
// returns lower index
long grid_search( long n, double quarry, GridPoint * A)
{
	long lowerLimit = 0;
	long upperLimit = n-1;
	long examinationPoint;
	long length = upperLimit - lowerLimit;

	while( length > 1 )
	{
		examinationPoint = lowerLimit + ( length / 2 );
		
		if( A[examinationPoint].energy > quarry )
			upperLimit = examinationPoint;
		else
			lowerLimit = examinationPoint;
		
		length = upperLimit - lowerLimit;
	}
	
	return lowerLimit;
}
