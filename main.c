#include "pd.h"

List           *cell;

sim_contact_mechanics contact_mechanics = PLASTIC;

long long	number_of_time_steps;

double particle_norm_old_x[number_of_interaction_pairs], particle_norm_old_y[number_of_interaction_pairs], particle_norm_old_z[number_of_interaction_pairs];
double particle_Rp[number_of_interaction_pairs], particle_disp_normal_max[number_of_interaction_pairs], particle_force_normal_max[number_of_interaction_pairs];
double particle_force_tan_old_mag[number_of_interaction_pairs], particle_force_tan_old_x[number_of_interaction_pairs], particle_force_tan_old_y[number_of_interaction_pairs], particle_force_tan_old_z[number_of_interaction_pairs];
double particle_disp_tan_old_mag[number_of_interaction_pairs], particle_disp_tan_old_x[number_of_interaction_pairs], particle_disp_tan_old_y[number_of_interaction_pairs], particle_disp_tan_old_z[number_of_interaction_pairs];
double particle_sign_of_disp[number_of_interaction_pairs];
long long	particle_number_of_time_steps[number_of_interaction_pairs], particle_initial_time_step[number_of_interaction_pairs];

int particle_nextx[number_of_interaction_pairs], particle_nexty[number_of_interaction_pairs], particle_nextz[number_of_interaction_pairs];

int particle_color[number_of_particles];
 
double particle_x[number_of_particles], particle_y[number_of_particles], particle_z[number_of_particles];
double particle_radius[number_of_particles];
double particle_YS[number_of_particles], particle_E_star[number_of_particles], particle_G_star[number_of_particles];	
double particle_mass[number_of_particles], particle_d_mass[number_of_particles], particle_moment[number_of_particles];
double particle_velocity_x[number_of_particles], particle_velocity_y[number_of_particles], particle_velocity_z[number_of_particles];
double particle_force_x[number_of_particles], particle_force_y[number_of_particles], particle_force_z[number_of_particles];
double particle_torque_x[number_of_particles], particle_torque_y[number_of_particles], particle_torque_z[number_of_particles];
double particle_angular_velocity_x[number_of_particles], particle_angular_velocity_y[number_of_particles], particle_angular_velocity_z[number_of_particles];

int     number_of_cells;
int     number_x, number_y, number_z;

double	accumulated_angle, angle_sign;
double	bottom_wall_force, damping, friction_coefficient;
double 	dt, TIME, MASS, LENGTH;
double 	elapsed_time, particle_diameter;
double wall_force = 0.0, wall_speed = WALL_SPEED, wall_mass = WALL_MASS;

double          min_x, max_x = 0.0;
double          min_y, max_y = 0.0;
double          min_z, max_z = 0.0;

double
Min(double x, double y)
{
	if (x < y)
		return (x);
	else
		return (y);
}

double
Max(double x, double y)
{
	if (x > y)
		return (x);
	else
		return (y);
}

void
usage()
{
	printf("\n");
	fprintf(stderr, "usage:  PD -i (name) [h] \n\
	lattice boltzmann code\n\
\n\
-i (file)\t input file (required)\n(format: x y z solid_node_boolean)\n\
\n OPTIONAL arguments\n\
-h   \t show this message  \n\n");

}

int
main(int argc, char **argv)
{
	int m, number_of_files = 0;
	double running_time=2000.0, write_time = 0.0;
	char *infile = NULL;
	char outfile[256];
	FILE           *fpout;
	
	/** I/O handling **/

	extern char    *optarg;
	int             option_letter = 0, usage_flag = 0;

	/** parse commmand line options **/

	while ((option_letter = getopt(argc, argv, "i:h")) != -1) {

		if (option_letter == 'i') {
			infile = optarg;
			usage_flag++;
		} else if (option_letter == 'h') {
			usage();
			exit(1);
		}
	}

	if (usage_flag < 0) {
		usage();
		exit(1);
	}

	angle_sign = 1.0;
	accumulated_angle = 0.0;
	particle_diameter = 0.005;
	
	number_of_time_steps=0;
	
	initialize(infile);
	
	/* info for moving wall simulations */

	wall_speed = wall_speed * TIME / LENGTH;
	wall_mass = wall_mass / MASS;		
	
	while (elapsed_time <= running_time * (1.0 / TIME)) {
		
		min_x = (double) Nx;
		max_x = 0.0;
		min_y = (double) Ny;
		max_y = 0.0;
		min_z = (double) Nz;
		max_z = 0.0;

		for (m=number_of_wall_particles;m<number_of_particles;m++) {
			double px, py, pz;
				
			particle_velocity_x[m] += update_velocity(particle_mass[m], particle_force_x[m]);
			particle_velocity_y[m] += update_velocity(particle_mass[m], particle_force_y[m]);
			particle_velocity_z[m] += update_velocity(particle_mass[m], particle_force_z[m]);
				
			px = particle_x[m] += move_particle(particle_mass[m], particle_velocity_x[m], particle_force_x[m]);
			if (px < 0.0)
				particle_x[m] = px + (double)Nx;
			else if (px >= Nx)
				particle_x[m] = px - (double)Nx;
				
			py = particle_y[m] += move_particle(particle_mass[m], particle_velocity_y[m], particle_force_y[m]);
			if (py < 0.0)
				particle_y[m] = py + (double)Ny;
			else if (py >= Ny)
				particle_y[m] = py - (double)Ny;
				
			pz = particle_z[m] += move_particle(particle_mass[m], particle_velocity_z[m], particle_force_z[m]);
			if (pz < 0.0)
				particle_z[m] = pz + (double)Nz;
			else if (pz >= Nz)
				particle_z[m] = pz - (double)Nz;
				
			particle_angular_velocity_x[m] += update_velocity(particle_moment[m], particle_torque_x[m]);
			particle_angular_velocity_y[m] += update_velocity(particle_moment[m], particle_torque_y[m]);
			particle_angular_velocity_z[m] += update_velocity(particle_moment[m], particle_torque_z[m]);
				
			particle_force_x[m] = 0.0;
			particle_force_z[m] = 0.0;
			particle_force_y[m] = -particle_mass[m];
			
			particle_torque_x[m] = particle_torque_y[m] = particle_torque_z[m] = 0.0;
			
			if (px > max_x)
				max_x = px;
			if (py > max_y)
				max_y = py;
			if (pz > max_z)
				max_z = pz;
			if (px < min_x)
				min_x = px;
			if (py < min_y)
				min_y = py;
			if (pz < min_z)
				min_z = pz;
			
		}
		
		max_x = Min(max_x + 1.0, (double) Nx);
		max_y = Min(max_y + 1.0, (double) Ny);
		max_z = Min(max_z + 1.0, (double) Nz);
		min_x = Max(min_x - 1.0, 0.0);
		min_y = Max(min_y - 1.0, 0.0);
		min_z = Max(min_z - 1.0, 0.0);
		
		/** move walls and calculate wall stresses/temps, if necessary **/

		handle_walls();
		
		/** find potential collisions/contacts between particles **/
		contact_detect();	/** then loop over cells and add particles to potential collision list **/
	
		for (m=number_of_wall_particles;m<number_of_particles;m++) {
		/** if there are solid walls **/
			if (XPERIODIC != 1 ) {
				calculate_contact_force(number_of_particles+left, m); /* use one of these for each wall */
				calculate_contact_force(number_of_particles+right, m); /* use one of these for each wall */
			}
				
			if (YPERIODIC != 1 ) {
				calculate_contact_force(number_of_particles+top, m); /* use one of these for each wall */
				calculate_contact_force(number_of_particles+bottom, m); /* use one of these for each wall */
			}
				
			if (ZPERIODIC != 1 ) {
				calculate_contact_force(number_of_particles+front, m); /* use one of these for each wall */
				calculate_contact_force(number_of_particles+back, m); /* use one of these for each wall */
			}
				
		}
		
		/** move walls and calculate wall stresses/temps, if necesary **/
		handle_moving_walls();
	
		if (elapsed_time >= write_time) {
			sprintf(outfile, "file%04d", number_of_files);
			if ((fpout = fopen(outfile, "w")) == NULL) {
				printf("Cannot Open File\n");
				exit(1);
			}
			for (m = 0; m < number_of_particles; m++)
				fprintf(fpout, "%e %e %e %e %e %e %e %e %e %e 0.0 %d\n", particle_x[m], particle_y[m], particle_z[m], particle_velocity_x[m], particle_velocity_y[m], particle_velocity_z[m], particle_angular_velocity_x[m], particle_angular_velocity_y[m],	particle_angular_velocity_z[m], particle_radius[m], particle_color[m]);
			
			fclose(fpout);
						
			number_of_files++;
			write_time += (1.0/30.0) * (1.0 / TIME);
			
		}
			
		elapsed_time += dt;
		number_of_time_steps++;

	}		
	
	return(0);
} 
