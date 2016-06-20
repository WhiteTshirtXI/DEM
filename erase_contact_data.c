#include "pd.h"

void 
erase_contact_data(int q)
{
	particle_touching[q] = 0;
	particle_contact_radius[q] = 0.0;
	particle_force_tan_old_x[q] =
		particle_force_tan_old_y[q] =
		particle_force_tan_old_z[q] = 0.0;
	particle_force_normal_old[q] =
		particle_force_normal_max[q] =
		particle_force_normal_status[q] =
		particle_Rp[q] = 0.0;
}
