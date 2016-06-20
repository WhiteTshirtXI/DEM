#include "pd.h"

void
initialize(char *infile)
{
	int		m, q;
	double junk;
	double radius=0.0;
	double gravity = 9.81;	/* m/s^2 */
	double youngs_modulus, yield_stress, poisson_ratio, density, shear_modulus, lambda, particle_diameter_min, conductivity, heat_capacity, thermal_expansion_coef;
    double L_density, H_density;
    double dt_old;
	FILE *fpin;
	
    get_material_properties(material, &thermal_expansion_coef, &youngs_modulus, &density, &poisson_ratio, &yield_stress, &conductivity, &heat_capacity);
    
    /** calculate scaling variables **/
	
	LENGTH = particle_diameter;
	MASS = (4.0 / 3.0) * M_PI * density * pow((LENGTH / 2.0), 3.0);
	TIME = sqrt(LENGTH / gravity);
	particle_diameter_min = LENGTH;
	
    /** 0.25* to make it insensitive **/
    dt_old = 1.0;
    
	/* read in wall particle details */

	if ((fpin = fopen(infile, "r")) == NULL) {
		printf("Cannot Open File\n");
		exit(1);
	}
	
    sim_material LIGHT_PARTICLE = ACETATE;
    //sim_material LIGHT_PARTICLE = ACETATE;
    sim_material HEAVY_PARTICLE = GLASS;
    
    get_material_properties(LIGHT_PARTICLE, &thermal_expansion_coef, &youngs_modulus, &density, &poisson_ratio, &yield_stress, &conductivity, &heat_capacity);
    L_density=density;
    
    get_material_properties(HEAVY_PARTICLE, &thermal_expansion_coef, &youngs_modulus, &density, &poisson_ratio, &yield_stress, &conductivity, &heat_capacity);
    H_density=density;
    
	for (m = 0; m < number_of_particles; m++) {
		/* read solid particle positions in from file with format of: x y z vx vy vz angvx angvy angvz r temp color_int */
		fscanf(fpin, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %i\n", &(particle_x[m]), &(particle_y[m]), &(particle_z[m]), &(particle_velocity_x[m]), &(particle_velocity_y[m]), &(particle_velocity_z[m]), &(particle_angular_velocity_x[m]), &(particle_angular_velocity_y[m]), &(particle_angular_velocity_z[m]), &(particle_radius[m]), &(junk), &(particle_color[m]));
		
        sim_material this_material;
        
        if (particle_color[m]>5)
        {
            this_material=LIGHT_PARTICLE;
        }
        else if (particle_color[m]<4)
        {
            this_material=HEAVY_PARTICLE;
        }
        else this_material=SOFT;
        
        get_material_properties(this_material, &thermal_expansion_coef, &youngs_modulus, &density, &poisson_ratio, &yield_stress, &conductivity, &heat_capacity);
        
		radius = particle_radius[m] * LENGTH;	/* make radius dimensional */

		if (2.0 * radius < particle_diameter_min)
			particle_diameter_min = 2.0 * radius;

		particle_mass[m] = (4.0 / 3.0) * M_PI * pow(radius, 3.0) * density / MASS;

		particle_moment[m]=(2.0/5.0)*particle_mass[m]*particle_radius[m]*particle_radius[m];
		particle_force_x[m] = particle_force_y[m] = particle_force_z[m] = 0.0;
        particle_force_y_ng[m]=(4.0 / 3.0) * M_PI * pow(radius, 3.0) * (H_density-L_density) / MASS;
		particle_torque_x[m] = particle_torque_y[m] = particle_torque_z[m] = 0.0;
        //printf("particle_force_y_ng is %lf and m %i \n", particle_force_y_ng[m], m);
        shear_modulus = 0.5 * youngs_modulus / (1.0 + poisson_ratio);
        lambda = 0.1631 * poisson_ratio + 0.8766;

        /** calculate time-step using dimensional parameters (from Kafui, 2002) **/
        dt = 0.02 * ((M_PI * particle_diameter_min / (2.0 * lambda)) * sqrt(density / shear_modulus)) / TIME;  // I changed the position of this function to get the minimum dt value
        
        if (dt < dt_old)
            dt_old = dt;
        
        /** make collision parameters dimensionless for each particle **/
        particle_G_star[m] = (shear_modulus * (LENGTH * TIME * TIME / MASS)) /
        (2.0 - poisson_ratio);
        particle_E_star[m] = (youngs_modulus * (LENGTH * TIME * TIME / MASS)) /
        (1.0 - poisson_ratio * poisson_ratio);
        particle_YS[m] = yield_stress * (LENGTH * TIME * TIME / MASS);

        
		// Set the "last" touching time step to be so small that we will erase contact data even for a particle pair that is initially touching
		for (q = m+1; q < number_of_particles; q++) particle_number_of_time_steps[pair(q,m)] = -2;
	}
	

    printf("Time step dt is %e s. Write_time is %e s.  Scaling variables are dt=%e MASS=%e TIME=%e LENGTH=%e\n", dt, (1.0 / 30.0) * (1.0 / TIME), dt*TIME, MASS, TIME, LENGTH);

    dt = dt_old;
    //dt=0.000029;
    printf("value of dt with unit of second %e\n", dt * TIME);
    
    printf("L_density_value %lf // H_density_value %lf\n", L_density, H_density);
	

	
	damping = 0.001875 * (youngs_modulus * (LENGTH * TIME * TIME / MASS)) /
		(1.0 - poisson_ratio * poisson_ratio);  //this damping is slightly different than the old version. Is there a reason?
    
	
	fclose(fpin);
	
}	
