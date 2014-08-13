#include "XSbench_header.h"

// Calculates the microscopic cross section for a given nuclide & energy
void calculate_micro_xs(   double p_energy, int nuc, long n_isotopes,
                           long n_gridpoints,
                           GridPoint * restrict energy_grid,
                           NuclideGridPoint ** restrict nuclide_grids,
                           int idx, double * restrict xs_vector ){
	
	// Variables
	double f;
	NuclideGridPoint * low, * high;

	// pull ptr from energy grid and check to ensure that
	// we're not reading off the end of the nuclide's grid
	if( energy_grid[idx].xs_ptrs[nuc] == n_gridpoints - 1 )
		low = &nuclide_grids[nuc][energy_grid[idx].xs_ptrs[nuc] - 1];
	else
		low = &nuclide_grids[nuc][energy_grid[idx].xs_ptrs[nuc]];
	
	high = low + 1;
	
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

// Struct-of-Arrays version:
void calculate_micro_xs_SOA(   double p_energy, int nuc, long n_isotopes,
                           long n_gridpoints,
                           GridPoint * restrict energy_grid,
                           GridPoint_SOA * restrict energy_grid_SOA,
                           NuclideGridPoint ** restrict nuclide_grids,
                           NuclideGridPoint_SOA * restrict nuc_grids_SOA,
                           int idx, double * restrict xs_vector ){
	
	// Variables
	double f;
	NuclideGridPoint * low, * high;

  // TODO: get rid of this by duplicating last grid point?
	// pull ptr from energy grid and check to ensure that
	// we're not reading off the end of the nuclide's grid
	if( energy_grid[idx].xs_ptrs[nuc] == n_gridpoints - 1 ) {
		low = &nuclide_grids[nuc][energy_grid[idx].xs_ptrs[nuc] - 1];
		//low = &nuclide_grids[nuc][energy_grid_SOA->xs_ptrs[idx*n_isotopes + nuc - 1]];
		//double low1 = nuc_grids_SOA->energy[nuc*n_gridpoints + energy_grid[idx].xs_ptrs[nuc]];
  }
	else {
		low = &nuclide_grids[nuc][energy_grid[idx].xs_ptrs[nuc]];
		//low = &nuclide_grids[nuc][energy_grid_SOA->xs_ptrs[idx*n_isotopes + nuc]];
		//double low1 = nuc_grids_SOA->energy[nuc*n_gridpoints + energy_grid[idx].xs_ptrs[nuc]];
    //printf("myegrid[idx] = %d\n", energy_grid_SOA->xs_ptrs[idx*n_isotopes + nuc]);
  }
	
  int j = energy_grid_SOA->xs_ptrs[idx*n_isotopes + nuc];
  int nuc_idx = nuc*n_gridpoints;

	//high = low + 1;
  int hi = nuc_idx+j+1;
  int lo = nuc_idx+j;
	
	// calculate the re-useable interpolation factor
	//f = (high->energy - p_energy) / (high->energy - low->energy);
	f = (nuc_grids_SOA->energy[hi] - p_energy) / (nuc_grids_SOA->energy[hi] - nuc_grids_SOA->energy[lo]);

	// Total XS
	//xs_vector[0] = high->total_xs - f * (high->total_xs - low->total_xs);
	xs_vector[0] = nuc_grids_SOA->total_xs[hi] - f * (nuc_grids_SOA->total_xs[hi] - nuc_grids_SOA->total_xs[lo]);
	
	// Elastic XS
	//xs_vector[1] = high->elastic_xs - f * (high->elastic_xs - low->elastic_xs);
	xs_vector[1] = nuc_grids_SOA->elastic_xs[hi] - f * (nuc_grids_SOA->elastic_xs[hi] - nuc_grids_SOA->elastic_xs[lo]);
	
	// Absorbtion XS
	//xs_vector[2] = high->absorbtion_xs - f * (high->absorbtion_xs - low->absorbtion_xs);
	xs_vector[2] = nuc_grids_SOA->absorbtion_xs[hi] - f * (nuc_grids_SOA->absorbtion_xs[hi] - nuc_grids_SOA->absorbtion_xs[lo]);
	
	// Fission XS
	//xs_vector[3] = high->fission_xs - f * (high->fission_xs - low->fission_xs);
	xs_vector[3] = nuc_grids_SOA->fission_xs[hi] - f * (nuc_grids_SOA->fission_xs[hi] - nuc_grids_SOA->fission_xs[lo]);
	
	// Nu Fission XS
	//xs_vector[4] = high->nu_fission_xs - f * (high->nu_fission_xs - low->nu_fission_xs);
	xs_vector[4] = nuc_grids_SOA->nu_fission_xs[hi] - f * (nuc_grids_SOA->nu_fission_xs[hi] - nuc_grids_SOA->nu_fission_xs[lo]);
	
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
	double xs_vector[5];
	int p_nuc; // the nuclide we are looking up
	long idx = 0;	
	double conc; // the concentration of the nuclide in the material

	// cleans out macro_xs_vector
	for( int k = 0; k < 5; k++ )
		macro_xs_vector[k] = 0;

	// binary search for energy on unionized energy grid (UEG)
	idx = grid_search( n_isotopes * n_gridpoints, p_energy,
	                   energy_grid);	
	
	// Once we find the pointer array on the UEG, we can pull the data
	// from the respective nuclide grids, as well as the nuclide
	// concentration data for the material
	// Each nuclide from the material needs to have its micro-XS array
	// looked up & interpolatied (via calculate_micro_xs). Then, the
	// micro XS is multiplied by the concentration of that nuclide
	// in the material, and added to the total macro XS array.
	for( int j = 0; j < num_nucs[mat]; j++ )
	{
		p_nuc = mats[mat][j];
		conc = concs[mat][j];
		calculate_micro_xs( p_energy, p_nuc, n_isotopes,
		                    n_gridpoints, energy_grid,
		                    nuclide_grids, idx, xs_vector );
		for( int k = 0; k < 5; k++ )
			macro_xs_vector[k] += xs_vector[k] * conc;
	}
	
	//test
	/*
	for( int k = 0; k < 5; k++ )
		printf("Energy: %lf, Material: %d, XSVector[%d]: %lf\n",
		       p_energy, mat, k, macro_xs_vector[k]);
	*/
}

// Struct-of-Array version
void calculate_macro_xs_SOA( double p_energy, int mat, long n_isotopes,
                         long n_gridpoints, int * restrict num_nucs,
                         double ** restrict concs,
                         GridPoint * restrict energy_grid,
                         GridPoint_SOA * restrict energy_grid_SOA,
                         NuclideGridPoint ** restrict nuclide_grids,
                         NuclideGridPoint_SOA * restrict nuclide_grids_SOA,
                         int ** restrict mats,
                         double * restrict macro_xs_vector ){
	double xs_vector[5];
	int p_nuc; // the nuclide we are looking up
	long idx = 0;	
	double conc; // the concentration of the nuclide in the material

	// cleans out macro_xs_vector
	for( int k = 0; k < 5; k++ )
		macro_xs_vector[k] = 0;

	// binary search for energy on unionized energy grid (UEG)
	idx = grid_search_SOA( n_isotopes * n_gridpoints, p_energy, energy_grid_SOA );	
	
	// Once we find the pointer array on the UEG, we can pull the data
	// from the respective nuclide grids, as well as the nuclide
	// concentration data for the material
	// Each nuclide from the material needs to have its micro-XS array
	// looked up & interpolatied (via calculate_micro_xs). Then, the
	// micro XS is multiplied by the concentration of that nuclide
	// in the material, and added to the total macro XS array.
#pragma simd
//#pragma novector
	for( int j = 0; j < num_nucs[mat]; j++ )
	{
		p_nuc = mats[mat][j];
		conc = concs[mat][j];
		//calculate_micro_xs( p_energy, p_nuc, n_isotopes,
		//                    n_gridpoints, energy_grid,
		//                    nuclide_grids, idx, xs_vector );
		calculate_micro_xs_SOA( p_energy, p_nuc, n_isotopes,
		                    n_gridpoints, energy_grid, energy_grid_SOA,
		                    nuclide_grids, nuclide_grids_SOA, idx, xs_vector );
		for( int k = 0; k < 5; k++ )
			macro_xs_vector[k] += xs_vector[k] * conc;
	}
	
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

long grid_search_SOA( long n, double quarry, GridPoint_SOA * A)
{
	long lowerLimit = 0;
	long upperLimit = n-1;
	long examinationPoint;
	long length = upperLimit - lowerLimit;

	while( length > 1 )
	{
		examinationPoint = lowerLimit + ( length / 2 );
		
		if( A->energy[examinationPoint] > quarry )
			upperLimit = examinationPoint;
		else
			lowerLimit = examinationPoint;
		
		length = upperLimit - lowerLimit;
	}
	
	return lowerLimit;
}
