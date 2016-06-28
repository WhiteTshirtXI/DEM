#include "pd.h"

void
contact_detect()
{
    Nx_index=0;
    W_Nx_index=0;
    
	int             i, j, icell, jcell, counter;
	for (icell = 0; icell < number_of_cells; icell++) {
		if (cell[icell].active) {
			i = cell[icell].begin_of_list;

			while (i >= 0) {

				j = particle_list[i];

				while (j >= 0) {

                    if (i > j && i>=number_of_wall_particles) //this prevents double calculation. (i,j) & (j,i)
                    {
						calculate_contact_force(i, j); 
                    }
					else if (j>=number_of_wall_particles)
						calculate_contact_force(j, i);
					j = particle_list[j];  // go back to one particle in the same cell. Do the contact force calculation again   --- loop

				}

				for (counter = 0; counter < CELL_COORD_NUMBER; counter++) {

					jcell = cell[icell].neighbors[counter];

					if (jcell >= 0 && cell[jcell].active) {  //active means the cell has particles inside it

						j = cell[jcell].begin_of_list;

						while (j >= 0) {

							if (i > j && i>=number_of_wall_particles)
								calculate_contact_force(i, j);  //calculate contact force between one particle in a cell and anotehr particle in its neighbor's cell
							else if (j>=number_of_wall_particles)
								calculate_contact_force(j, i);
							j = particle_list[j];

						}

					}
				}

				i = particle_list[i];  //go to one privious particle in the cell

			}

		}
	}
    if (cohesive==0) {
        fprintf(fpout_Nx_index, "Time(s) %e Contact_number %d number_of_time_steps %lli\n", elapsed_time*TIME, Nx_index, number_of_time_steps);
    }
    else if (cohesive==1){
                fprintf(fpout_Nx_index, "Time(s) %e Contact_number %d Wet_Contact_number %d number_of_time_steps %lli\n", elapsed_time*TIME, Nx_index, W_Nx_index, number_of_time_steps);
    }
}
