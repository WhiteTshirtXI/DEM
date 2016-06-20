#include "pd.h"
void
get_cell_neighbors()
{
	int             i, l, jx, jy, jz, kx, ky, kz;
	int             cellx, celly, cellz;
	
	/** figure cell neighbors **/
	/** zero cell number counter **/

	i = 0;

	/** loop over all cell coordinates **/

	for (celly = 0; celly < number_y; celly++) {
		for (cellz = 0; cellz < number_z; cellz++) {
			for (cellx = 0; cellx < number_x; cellx++) {

				/** zero neighbor number counter and remember the cell coordinates **/
				
				l = 0;
				cell[i].x = cellx;
				cell[i].y = celly;
				cell[i].z = cellz;

				for (kz = -1; kz <= 1; kz++) {
					for (ky = -1; ky <= 1; ky++) {
						for (kx = -1; kx <= 1; kx++) {
							if (kx != 0 || ky != 0 || kz != 0) {
								jx = cellx + kx;
								jy = celly + ky;
								jz = cellz + kz;

								if (XPERIODIC) {
									if (jx >= number_x)
										jx -= number_x;
									else if (jx < 0)
										jx += number_x;
								}
								if (YPERIODIC) {
									if (jy >= number_y)
										jy -= number_y;
									else if (jy < 0)
										jy += number_y;
								}
								if (ZPERIODIC) {
									if (jz >= number_z)
										jz -= number_z;
									else if (jz < 0)
										jz += number_z;
								}
								if (jx >= 0 && jx < number_x && jy >= 0 && jy < number_y && jz >= 0 && jz < number_z) {
									if (l < 13)
										cell[i].neighbors[l] = jx + jz * number_x + jy * number_x * number_z;
									else
										cell[i].anti_neighbors[l - 13] = jx + jz * number_x + jy * number_x * number_z;
								} else {  //when will this happen?
									if (l < 13)
										cell[i].neighbors[l] = -1;
									else
										cell[i].anti_neighbors[l - 13] = -1;
								}
								l++;
							}
						}
					}
				}
				i++;
			}
		}
	}
}
