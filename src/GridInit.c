#include "XSbench_header.h"

#ifdef MPI
#include<mpi.h>
#endif

// Generates randomized energy grid for each nuclide
// Note that this is done as part of initialization (serial), so
// rand() is used.
void generate_grids( NuclideGridPoint ** nuclide_grids,
                     long n_isotopes, long n_gridpoints ) {
	for( long i = 0; i < n_isotopes; i++ )
		for( long j = 0; j < n_gridpoints; j++ )
		{
			nuclide_grids[i][j].energy       =((double)rand()/(double)RAND_MAX);
			nuclide_grids[i][j].total_xs     =((double)rand()/(double)RAND_MAX);
			nuclide_grids[i][j].elastic_xs   =((double)rand()/(double)RAND_MAX);
			nuclide_grids[i][j].absorbtion_xs=((double)rand()/(double)RAND_MAX);
			nuclide_grids[i][j].fission_xs   =((double)rand()/(double)RAND_MAX);
			nuclide_grids[i][j].nu_fission_xs=((double)rand()/(double)RAND_MAX);
		}
}

// Verification version of this function (tighter control over RNG)
void generate_grids_v( NuclideGridPoint ** nuclide_grids,
                     long n_isotopes, long n_gridpoints ) {
	for( long i = 0; i < n_isotopes; i++ )
		for( long j = 0; j < n_gridpoints; j++ )
		{
			nuclide_grids[i][j].energy       = rn_v();
			nuclide_grids[i][j].total_xs     = rn_v();
			nuclide_grids[i][j].elastic_xs   = rn_v();
			nuclide_grids[i][j].absorbtion_xs= rn_v();
			nuclide_grids[i][j].fission_xs   = rn_v();
			nuclide_grids[i][j].nu_fission_xs= rn_v();
		}
}

// Verification version with MIC-friendly struct-of-array data structure
void generate_grids_v_SOA( NuclideGridPoint_SOA * nuclide_grids,
                     long n_isotopes, long n_gridpoints ) {

  	double * energy = (double *) _mm_malloc( n_isotopes * n_gridpoints * sizeof(double), 64 );
  	double * total_xs = (double *) _mm_malloc( n_isotopes * n_gridpoints * sizeof(double), 64 );
  	double * elastic_xs = (double *) _mm_malloc( n_isotopes * n_gridpoints * sizeof(double), 64 );
  	double * absorbtion_xs = (double *) _mm_malloc( n_isotopes * n_gridpoints * sizeof(double), 64 );
  	double * fission_xs = (double *) _mm_malloc( n_isotopes * n_gridpoints * sizeof(double), 64 );
  	double * nu_fission_xs = (double *) _mm_malloc( n_isotopes * n_gridpoints * sizeof(double), 64 );

    nuclide_grids->energy = energy;
    nuclide_grids->total_xs = total_xs;
    nuclide_grids->elastic_xs = elastic_xs;
    nuclide_grids->absorbtion_xs = absorbtion_xs;
    nuclide_grids->fission_xs = fission_xs;
    nuclide_grids->nu_fission_xs = nu_fission_xs;

	  for( long i = 0; i < n_isotopes; i++ ) 
    {
	  	for( long j = 0; j < n_gridpoints; j++ )
	  	{
	  		nuclide_grids->energy[i*n_isotopes + j]        = rn_v();
	  		nuclide_grids->total_xs[i*n_isotopes + j]      = rn_v();
	  		nuclide_grids->elastic_xs[i*n_isotopes + j]    = rn_v();
	  		nuclide_grids->absorbtion_xs[i*n_isotopes + j] = rn_v();
	  		nuclide_grids->fission_xs[i*n_isotopes + j]    = rn_v();
	  		nuclide_grids->nu_fission_xs[i*n_isotopes + j] = rn_v();
	  	}
    }
}

// Sorts the nuclide grids by energy (lowest -> highest)
void sort_nuclide_grids( NuclideGridPoint ** nuclide_grids, long n_isotopes,
                         long n_gridpoints )
{
	int (*cmp) (const void *, const void *);
	cmp = NGP_compare;
	
	for( long i = 0; i < n_isotopes; i++ )
		qsort( nuclide_grids[i], n_gridpoints, sizeof(NuclideGridPoint),
		       cmp );
	
	// error debug check
	//for( int i = 0; i < n_isotopes; i++ )
	//{
	//	printf("NUCLIDE %d ==============================\n", i);
	//	for( int j = 0; j < n_gridpoints; j++ )
	//		printf("e%d = %lf\n", j, nuclide_grids[i][j].energy);
	//}
}

// Struct-of-Array version
void sort_nuclide_grids_SOA( NuclideGridPoint_SOA * nuclide_grids, long n_isotopes,
                         long n_gridpoints )
{
	int (*cmp) (const void *, const void *);
	cmp = compare_double;

	for( long i = 0; i < n_isotopes; i++ )
  {

    int idx = i*n_gridpoints;
		qsort_SOA( &nuclide_grids->energy[idx], 
               &nuclide_grids->total_xs[idx],
               &nuclide_grids->elastic_xs[idx],
               &nuclide_grids->absorbtion_xs[idx],
               &nuclide_grids->fission_xs[idx],
               &nuclide_grids->nu_fission_xs[idx],
               n_gridpoints, sizeof(double), cmp );
  }
	// error debug check
	//for( int i = 0; i < n_isotopes; i++ )
	//{
	//	printf("NUCLIDE %d ==============================\n", i);
	//	for( int j = 0; j < n_gridpoints; j++ )
	//		printf("E%d = %lf\n", j, nuclide_grids->energy[i*n_gridpoints+j]);
	//}
}

// Allocates unionized energy grid, and assigns union of energy levels
// from nuclide grids to it.
GridPoint * generate_energy_grid( long n_isotopes, long n_gridpoints,
                                  NuclideGridPoint ** nuclide_grids) {
	int mype = 0;

	#ifdef MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &mype);
	#endif
	
	if( mype == 0 ) printf("Generating Unionized Energy Grid...\n");
	
	long n_unionized_grid_points = n_isotopes*n_gridpoints;
	int (*cmp) (const void *, const void *);
	cmp = NGP_compare;
	
	GridPoint * energy_grid = (GridPoint *)malloc( n_unionized_grid_points
	                                               * sizeof( GridPoint ) );
	if( mype == 0 ) printf("Copying and Sorting all nuclide grids...\n");
	
	NuclideGridPoint ** n_grid_sorted = gpmatrix( n_isotopes, n_gridpoints );
	
	  	
	memcpy( n_grid_sorted[0], nuclide_grids[0], n_isotopes*n_gridpoints*
	                                      sizeof( NuclideGridPoint ) );
	
	qsort( &n_grid_sorted[0][0], n_unionized_grid_points,
	       sizeof(NuclideGridPoint), cmp);

	if( mype == 0 ) printf("Assigning energies to unionized grid...\n");
	
	for( long i = 0; i < n_unionized_grid_points; i++ )
		energy_grid[i].energy = n_grid_sorted[0][i].energy;
	

	gpmatrix_free(n_grid_sorted);
	
	int * full = (int *) malloc( n_isotopes * n_unionized_grid_points
	                             * sizeof(int) );
	
	for( long i = 0; i < n_unionized_grid_points; i++ )
		energy_grid[i].xs_ptrs = &full[n_isotopes * i];
	
	// debug error checking
	/*
	for( int i = 0; i < n_unionized_grid_points; i++ )
		printf("E%d = %lf\n", i, energy_grid[i].energy);
	*/

	return energy_grid;
}

// Struct-of-Array version
GridPoint * generate_energy_grid_SOA( long n_isotopes, long n_gridpoints,
                                  NuclideGridPoint_SOA * nuclide_grids) {
	int mype = 0;

	#ifdef MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &mype);
	#endif
	
	if( mype == 0 ) printf("Generating Unionized Energy Grid...\n");
	
	long n_unionized_grid_points = n_isotopes*n_gridpoints;
	int (*cmp) (const void *, const void *);
	cmp = compare_double;
	
	NuclideGridPoint_SOA * n_grid_sorted = gpmatrix_SOA( n_isotopes, n_gridpoints );
	GridPoint * energy_grid = (GridPoint *)malloc( n_unionized_grid_points
	                                               * sizeof( GridPoint ) );
	if( mype == 0 ) printf("Copying and Sorting all nuclide grids...\n");
	
  double * energy = (double *) malloc( n_isotopes * n_gridpoints * sizeof(double) );
  double * total_xs = (double *) malloc( n_isotopes * n_gridpoints * sizeof(double) );
  double * elastic_xs = (double *) malloc( n_isotopes * n_gridpoints * sizeof(double) );
  double * absorbtion_xs = (double *) malloc( n_isotopes * n_gridpoints * sizeof(double) );
  double * fission_xs = (double *) malloc( n_isotopes * n_gridpoints * sizeof(double) );
  double * nu_fission_xs = (double *) malloc( n_isotopes * n_gridpoints * sizeof(double) );

  n_grid_sorted->total_xs = total_xs;
  n_grid_sorted->elastic_xs = elastic_xs;
  n_grid_sorted->absorbtion_xs = absorbtion_xs;
  n_grid_sorted->fission_xs = fission_xs;
  n_grid_sorted->nu_fission_xs = nu_fission_xs;
	  	
	memcpy( energy, nuclide_grids->energy, 
                                        n_isotopes*n_gridpoints*sizeof( double ) );
	memcpy( n_grid_sorted->total_xs, nuclide_grids->total_xs, 
                                        n_isotopes*n_gridpoints*sizeof( double ) );
	memcpy( n_grid_sorted->elastic_xs, nuclide_grids->elastic_xs, 
                                        n_isotopes*n_gridpoints*sizeof( double ) );
	memcpy( n_grid_sorted->absorbtion_xs, nuclide_grids->absorbtion_xs, 
                                        n_isotopes*n_gridpoints*sizeof( double ) );
	memcpy( n_grid_sorted->fission_xs, nuclide_grids->fission_xs, 
                                        n_isotopes*n_gridpoints*sizeof( double ) );
	memcpy( n_grid_sorted->nu_fission_xs, nuclide_grids->nu_fission_xs, 
                                        n_isotopes*n_gridpoints*sizeof( double ) );
	
	qsort_SOA( &energy[0], 
      &n_grid_sorted->total_xs[0],
      &n_grid_sorted->elastic_xs[0],
      &n_grid_sorted->absorbtion_xs[0],
      &n_grid_sorted->fission_xs[0],
      &n_grid_sorted->nu_fission_xs[0],
      n_unionized_grid_points, sizeof(double), cmp);

	if( mype == 0 ) printf("Assigning energies to unionized grid...\n");
	
	for( long i = 0; i < n_unionized_grid_points; i++ )
		energy_grid[i].energy = energy[i];

	int * full = (int *) malloc( n_isotopes * n_unionized_grid_points * sizeof(int) );
	
	for( long i = 0; i < n_unionized_grid_points; i++ )
		energy_grid[i].xs_ptrs = &full[n_isotopes * i];
	
	// debug error checking
	/*
	for( int i = 0; i < n_unionized_grid_points; i++ )
		printf("E%d = %lf\n", i, energy_grid[i].energy);
	*/

	return energy_grid;
}

void copy_AOS_to_SOA(GridPoint *energy_grid, GridPoint_SOA **energy_grid_SOA, 
           NuclideGridPoint **nuclide_grids, NuclideGridPoint_SOA *nuclide_grids_SOA, 
           int n_isotopes, int n_gridpoints) 
{

  double * energy = (double *) _mm_malloc( n_isotopes * n_gridpoints * sizeof(double), 64 );
  double * total_xs = (double *) _mm_malloc( n_isotopes * n_gridpoints * sizeof(double), 64 );
  double * elastic_xs = (double *) _mm_malloc( n_isotopes * n_gridpoints * sizeof(double), 64 );
  double * absorbtion_xs = (double *) _mm_malloc( n_isotopes * n_gridpoints * sizeof(double), 64 );
  double * fission_xs = (double *) _mm_malloc( n_isotopes * n_gridpoints * sizeof(double), 64 );
  double * nu_fission_xs = (double *) _mm_malloc( n_isotopes * n_gridpoints * sizeof(double), 64 );

  nuclide_grids_SOA->energy = energy;
  nuclide_grids_SOA->total_xs = total_xs;
  nuclide_grids_SOA->elastic_xs = elastic_xs;
  nuclide_grids_SOA->absorbtion_xs = absorbtion_xs;
  nuclide_grids_SOA->fission_xs = fission_xs;
  nuclide_grids_SOA->nu_fission_xs = nu_fission_xs;

  for (int i=0; i<n_isotopes; i++)
  {
    for (int j=0; j<n_gridpoints; j++)
    {
      nuclide_grids_SOA->energy[i*n_gridpoints + j] = nuclide_grids[i][j].energy;
      nuclide_grids_SOA->total_xs[i*n_gridpoints + j] = nuclide_grids[i][j].total_xs;
      nuclide_grids_SOA->elastic_xs[i*n_gridpoints + j] = nuclide_grids[i][j].elastic_xs;
      nuclide_grids_SOA->absorbtion_xs[i*n_gridpoints + j] = nuclide_grids[i][j].absorbtion_xs;
      nuclide_grids_SOA->fission_xs[i*n_gridpoints + j] = nuclide_grids[i][j].fission_xs;
      nuclide_grids_SOA->nu_fission_xs[i*n_gridpoints + j] = nuclide_grids[i][j].nu_fission_xs;
    }
  }

  long n_unionized_grid_points = n_isotopes*n_gridpoints;
  double * e_grid = (double *) _mm_malloc( n_unionized_grid_points * sizeof(double), 64 );
  int * xs_ptrs = (int *) _mm_malloc( n_unionized_grid_points * n_isotopes * sizeof(int), 64 );

  (*energy_grid_SOA) = (GridPoint_SOA *) _mm_malloc( sizeof(GridPoint_SOA), 64 );
  (*energy_grid_SOA)->energy = e_grid;
  (*energy_grid_SOA)->xs_ptrs = xs_ptrs;

	for( long i = 0; i < n_isotopes * n_gridpoints ; i++ )
  {
    (*energy_grid_SOA)->energy[i] = energy_grid[i].energy; 
    memcpy(&(*energy_grid_SOA)->xs_ptrs[i*n_isotopes], energy_grid[i].xs_ptrs, n_isotopes*sizeof(int));
  }
}

// Searches each nuclide grid for the closest energy level and assigns
// pointer from unionized grid to the correct spot in the nuclide grid.
// This process is time consuming, as the number of binary searches
// required is:  binary searches = n_gridpoints * n_isotopes^2
void set_grid_ptrs( GridPoint * energy_grid, NuclideGridPoint ** nuclide_grids,
                    long n_isotopes, long n_gridpoints )
{
	int mype = 0;

	#ifdef MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &mype);
	#endif
	
	if( mype == 0 ) printf("Assigning pointers to Unionized Energy Grid...\n");
	#pragma omp parallel for default(none) \
	shared( energy_grid, nuclide_grids, n_isotopes, n_gridpoints, mype )
	for( long i = 0; i < n_isotopes * n_gridpoints ; i++ )
	{
		double quarry = energy_grid[i].energy;
		if( INFO && mype == 0 && omp_get_thread_num() == 0 && i % 200 == 0 )
			printf("\rAligning Unionized Grid...(%.0lf%% complete)",
			       100.0 * (double) i / (n_isotopes*n_gridpoints /
				                         omp_get_num_threads())     );
		for( long j = 0; j < n_isotopes; j++ )
		{
			// j is the nuclide i.d.
			// log n binary search
			energy_grid[i].xs_ptrs[j] = 
				binary_search( nuclide_grids[j], quarry, n_gridpoints);
		}
	}
	if( mype == 0 ) printf("\n");

	//test
	/*
	for( int i=0; i < n_isotopes * n_gridpoints; i++ )
		for( int j = 0; j < n_isotopes; j++ )
			printf("E = %.4lf\tNuclide %d->%p->%.4lf\n",
			       energy_grid[i].energy,
                   j,
				   energy_grid[i].xs_ptrs[j],
				   (energy_grid[i].xs_ptrs[j])->energy
				   );
	*/
}

// Struct-of-Array version
void set_grid_ptrs_SOA( GridPoint * energy_grid, NuclideGridPoint_SOA * nuclide_grids,
                    long n_isotopes, long n_gridpoints )
{
	int mype = 0;

	#ifdef MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &mype);
	#endif
	
	if( mype == 0 ) printf("Assigning pointers to Unionized Energy Grid...\n");
	#pragma omp parallel for default(none) \
	shared( energy_grid, nuclide_grids, n_isotopes, n_gridpoints, mype )
	for( long i = 0; i < n_isotopes * n_gridpoints ; i++ )
	{
		double quarry = energy_grid[i].energy;
		if( INFO && mype == 0 && omp_get_thread_num() == 0 && i % 200 == 0 )
			printf("\rAligning Unionized Grid...(%.0lf%% complete)",
			       100.0 * (double) i / (n_isotopes*n_gridpoints /
				                         omp_get_num_threads())     );
		for( long j = 0; j < n_isotopes; j++ )
		{
			// j is the nuclide i.d.
			// log n binary search
			energy_grid[i].xs_ptrs[j] = 
				binary_search_SOA( &nuclide_grids->energy[j*n_gridpoints], quarry, n_gridpoints);
		}
	}
	if( mype == 0 ) printf("\n");

	////test
	///*
	//for( int i=0; i < n_isotopes * n_gridpoints; i++ )
	//	for( int j = 0; j < n_isotopes; j++ )
	//		printf("E = %.4lf\tNuclide %d->%p->%.4lf\n",
	//		       energy_grid[i].energy,
  //                 j,
	//			   energy_grid[i].xs_ptrs[j],
	//			   (energy_grid[i].xs_ptrs[j])->energy
	//			   );
	//*/
}
