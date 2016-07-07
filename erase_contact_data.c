#include "pd.h"

void 
erase_contact_data(int q)
{
	particle_initial_time_step[q] = number_of_time_steps;
	particle_sign_of_disp[q] = 1.0;
	particle_norm_old_x[q] =
		particle_norm_old_y[q] =
		particle_norm_old_z[q] = 0.0;
	particle_force_tan_old_mag[q] =
		particle_force_tan_old_x[q] =
		particle_force_tan_old_y[q] =
		particle_force_tan_old_z[q] = 0.0;
	particle_disp_tan_old_mag[q] =
		particle_disp_tan_old_x[q] =
		particle_disp_tan_old_y[q] =
		particle_disp_tan_old_z[q] = 0.0;
	particle_force_normal_max[q] =
		particle_disp_normal_max[q] =
		particle_Rp[q] = 0.0;
}
