#include "pd.h"

void
contact_detect()
{
	int             i, j, icell, jcell, counter;
	for (icell = 0; icell < number_of_cells; icell++) {
		if (cell[icell].active) {

			i = cell[icell].begin_of_list;

			while (i >= 0) {

				j = particle_list[i];

				while (j >= 0) {

                    if (i > j && i>=number_of_wall_particles)
                    {
						calculate_contact_force(i, j); //something wrong here.
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

}
