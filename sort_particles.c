#include "pd.h"

int 
sort_particles()  //cell[icell].begin_of_list (last particle in the cell) -> particle_list[last] (second last particle) -> particle_list[second last] (third last)  this is the link series  Particle_list always link to one prior particle in the same cell until eventually go back to -1 
{
	int             i, ix, iy, iz, icell;

	for (i = 0; i < number_of_cells; i++) {
		cell[i].begin_of_list = -1;
		cell[i].active = cell[i].toggled = 0;
	}

	for (i = 0; i < number_of_particles; i++) {
        
		ix = (int) ((particle_x[i]) / size_of_cells);
		iy = (int) ((particle_y[i]) / size_of_cells);
		iz = (int) ((particle_z[i]) / size_of_cells);
    //printf("WRONG 2 sort, number %i\n", i);
		icell = ix + iz * number_x + iy * number_x * number_z;
        
		particle_list[i] = cell[icell].begin_of_list;
		particle_cell_id[i] = icell;
		cell[icell].begin_of_list = i;
		cell[icell].active = 1;
		if (i>=number_of_wall_particles && !cell[icell].toggled) {
			int counter, jcell;
			
			cell[icell].toggled = 1;
			for (counter = 0; counter < CELL_COORD_NUMBER; counter++) {
				jcell = cell[icell].anti_neighbors[counter];  //jcell is icell's neighbor
				if (jcell>=0) cell[jcell].active = 1;  // is this necessary? seems active=1 for each particle
			}
		}

	}

	return (0);
}
