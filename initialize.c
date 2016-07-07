#include "pd.h"

void
initialize(char *infile)
{
	int		m, q;
	double junk;
	double radius=0.0;
	double gravity = 9.81;	/* m/s^2 */
	double youngs_modulus, yield_stress, poisson_ratio, density, shear_modulus, lambda, particle_diameter_min;
	FILE *fpin;
	
	/** default values for particle properties (assumes glass beads) **/

	density = 2700.0;	/* kg/m^3 */
	youngs_modulus = 68.95e9;	/* Pa=kg/m s^2		 */
	friction_coefficient = 0.30;
	poisson_ratio = 0.33;
	yield_stress = 68.95e7;	/* Pa = kg/m s^2	 */
	
	LENGTH = particle_diameter;
	MASS = (4.0 / 3.0) * M_PI * density * pow((LENGTH / 2.0), 3.0);
	TIME = sqrt(LENGTH / gravity);
	particle_diameter_min = LENGTH;
	
	/* read in wall particle details */

	if ((fpin = fopen(infile, "r")) == NULL) {
		printf("Cannot Open File\n");
		exit(1);
	}
	
	for (m = 0; m < number_of_particles; m++) {
		/* read solid particle positions in from file with format of: x y z vx vy vz angvx angvy angvz r temp color_int */
		fscanf(fpin, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %i\n", &(particle_x[m]), &(particle_y[m]), &(particle_z[m]), &(particle_velocity_x[m]), &(particle_velocity_y[m]), &(particle_velocity_z[m]), &(particle_angular_velocity_x[m]), &(particle_angular_velocity_y[m]), &(particle_angular_velocity_z[m]), &(particle_radius[m]), &(junk), &(particle_color[m]));
		
		radius = particle_radius[m] * LENGTH;	/* make radius dimensional */

		if (2.0 * radius < particle_diameter_min)
			particle_diameter_min = 2.0 * radius;

		particle_mass[m] = (4.0 / 3.0) * M_PI * pow(radius, 3.0) * density / MASS;

		if (particle_color[m]==6) particle_mass[m] *= DENSITY_RATIO;
		
		particle_moment[m]=(2.0/5.0)*particle_mass[m]*particle_radius[m]*particle_radius[m];
		particle_force_x[m] = particle_force_y[m] = particle_force_z[m] = 0.0;
		particle_torque_x[m] = particle_torque_y[m] = particle_torque_z[m] = 0.0;

		
		for (q = m+1; q < number_of_units; q++) {
			// Set the "last" touching time step to be so small that we will erase all contact data at the start
			particle_number_of_time_steps[pair(q,m)] = -2;
			
			// Set the "last" initial time step to be very small so that we catch contacts that are too short
			particle_initial_time_step[pair(q,m)] = -200;
		}
		
	}
	
	/** calculate time-step using dimensional parameters (from Kafui, 2002) **/

	shear_modulus = 0.5 * youngs_modulus / (1.0 + poisson_ratio);
	lambda = 0.1631 * poisson_ratio + 0.8766;

	/** 0.25* to make it insensitive **/
	dt = 0.125 * ((M_PI * particle_diameter_min / (2.0 * lambda)) * sqrt(density / shear_modulus))/TIME;
	
	printf("Time step is %e s. Scaling variables are MASS=%e TIME=%e LENGTH=%e\n", dt*TIME, MASS, TIME, LENGTH);
	/** make collision parameters dimensionless for each particle **/

	for (m = 0; m < number_of_particles; m++) {

		particle_G_star[m] = (shear_modulus * (LENGTH * TIME * TIME / MASS)) /
			(2.0 - poisson_ratio);
		particle_E_star[m] = (youngs_modulus * (LENGTH * TIME * TIME / MASS)) /
			(1.0 - poisson_ratio * poisson_ratio);
		particle_YS[m] = yield_stress * (LENGTH * TIME * TIME / MASS);

	}
	
	damping = 0.0004 * (youngs_modulus * (LENGTH * TIME * TIME / MASS)) /
		(1.0 - poisson_ratio * poisson_ratio);
	
	fclose(fpin);
	
}	
