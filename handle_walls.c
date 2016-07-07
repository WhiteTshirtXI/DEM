#include "pd.h"

void 
handle_walls()
{
	int             i;
	double			rx1, rx2, ry1, ry2, rot_sin, rot_cos, freq;
	double			XOFFSET, YOFFSET, dangle;

	freq = (FREQ * 2.0 * M_PI) * TIME;
	
	if (fabs(accumulated_angle+freq*angle_sign*dt)>MAX_ANGLE) angle_sign*=-1.0;
	
	accumulated_angle += freq*angle_sign*dt;
	
	dangle = freq*angle_sign*dt;
	
	rot_sin = sin(dangle);
	rot_cos = cos(dangle);

	for (i = 0; i < number_of_wall_particles; i++) {
		
		if ((double)i<(0.5*(double)number_of_wall_particles)) XOFFSET = XOFFSET_LEFT;
		else XOFFSET = XOFFSET_RIGHT;
	
		YOFFSET = wall_y(number_of_particles+bottom, particle_x[i], particle_y[i], particle_z[i]);
		
		rx1 = (particle_x[i] - XOFFSET) * rot_cos;
		rx2 = -(particle_y[i] - YOFFSET) * rot_sin;
		ry1 = (particle_x[i] - XOFFSET) * rot_sin;
		ry2 = (particle_y[i] - YOFFSET) * rot_cos;

		particle_x[i] = rx1 + rx2 + XOFFSET;
		particle_y[i] = ry1 + ry2 + YOFFSET;
		
		particle_velocity_x[i] = (particle_y[i] - YOFFSET) * freq;
		particle_velocity_y[i] = -(particle_x[i] - XOFFSET) * freq;
		particle_velocity_z[i] = 0.0;

		particle_angular_velocity_x[i] = 0.0;
		particle_angular_velocity_y[i] = 0.0;
		particle_angular_velocity_z[i] = 0.0;

		particle_force_x[i] = particle_force_y[i] = particle_force_z[i] = 0.0;
		particle_torque_x[i] = particle_torque_y[i] = particle_torque_z[i] = 0.0;
	}
	
}

void 
handle_moving_walls()
{

}
