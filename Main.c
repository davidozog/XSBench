#include "XSbench_header.h"

int main( int argc, char* argv[] )
{
	int n_isotopes = 68;
	int n_gridpoints = 10000;
	int lookups = 100000000;
	int i, thread, nthreads, mat;
	double omp_start, omp_end, p_energy;
	int max_procs = omp_get_num_procs();
	
	srand(time(NULL));

	if( argc == 2 )
		nthreads = atoi(argv[1]);
	else
		nthreads = max_procs;
	
	omp_set_num_threads(nthreads); 

	logo();
	
	// Allocate & fill energy grids
	if( DEBUG ) printf("Generating Nuclide Energy Grids...\n");
	
	NuclideGridPoint ** nuclide_grids = gpmatrix( n_isotopes, n_gridpoints );
	generate_grids( nuclide_grids, n_isotopes, n_gridpoints );	
	
	// Sort grids by energy
	sort_nuclide_grids( nuclide_grids, n_isotopes );

	// Prepare Unionized Energy Grid Framework
	GridPoint * energy_grid = generate_energy_grid( n_isotopes, n_gridpoints,
	                                                nuclide_grids ); 	

	// Double Indexing. Filling in energy_grid with pointers to the
	// nuclide_energy_grids.
	set_grid_ptrs( energy_grid, nuclide_grids, n_isotopes, n_gridpoints );
	
	// Get material data
	if( INFO ) printf("Loading Mats...\n");
	int *num_nucs = load_num_nucs();
	int **mats = load_mats(num_nucs);
	double **concs = load_concs(num_nucs);

	if( INFO ) printf("Using %d threads.\n", nthreads);

	omp_start = omp_get_wtime();
	
	// Energy grid built. Now to enter parallel region
	#pragma omp parallel default(none) \
	private(i, thread, p_energy, mat) \
	shared( max_procs, n_isotopes, n_gridpoints, \
	energy_grid, nuclide_grids, lookups, nthreads, \
	mats, concs, num_nucs)
	{	
		thread = omp_get_thread_num();

		#pragma omp for
		for( i = 0; i < lookups; i++ )
		{
			// Status text
			if( DEBUG && thread == 0 && i % 10000 == 0 )
				printf("\rCalculating XS's... (%.1lf%% completed)",
						i / ( lookups / (double) nthreads ) * 100.0);

			// Randomly pick an energy and material for the particle
			//p_energy = (double) rand() / (double) RAND_MAX;
			p_energy = rn();
			mat = pick_mat(); 
		
			// This returns the macro_xs, but we're not going to do anything
			// with it in this program, so return value is not stored.
			calculate_macro_xs( p_energy, mat, n_isotopes,
			                    n_gridpoints, num_nucs, concs,
			                    energy_grid, nuclide_grids, mats );
		}	
	}
	if( DEBUG ) printf("\n" );

	omp_end = omp_get_wtime();

	// Print the results
	if( INFO ) printf("Runtime:   %.3lf seconds\n", omp_end-omp_start);
	if( INFO ) printf("Lookups:   %d\n", lookups);
	if( INFO ) printf("Lookups/s: %.0lf\n",
		               (double) lookups / (omp_end-omp_start));

	return 0;
}
