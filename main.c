#include "pd.h"

List           *cell;

int particle_touching[number_of_interaction_pairs], particle_neighbor_id[number_of_interaction_pairs];	

sim_normal_contact_status particle_force_normal_status[number_of_interaction_pairs];
sim_contact_mechanics contact_mechanics = PLASTIC;

long long	number_of_time_steps;

double particle_Rp[number_of_interaction_pairs], particle_force_normal_old[number_of_interaction_pairs], particle_force_normal_max[number_of_interaction_pairs];
double particle_force_tan_old_x[number_of_interaction_pairs], particle_force_tan_old_y[number_of_interaction_pairs], particle_force_tan_old_z[number_of_interaction_pairs];
double particle_contact_radius[number_of_interaction_pairs];
long long	particle_number_of_time_steps[number_of_interaction_pairs];

int particle_color[number_of_particles];
int particle_nextx[number_of_particles], particle_nexty[number_of_particles], particle_nextz[number_of_particles];	
int particle_list[number_of_particles], particle_cell_id[number_of_particles];
int particle_coordination_number_old[number_of_particles], particle_coordination_number_new[number_of_particles];	
double particle_force_y_ng[number_of_particles];
double particle_x[number_of_particles], particle_y[number_of_particles], particle_z[number_of_particles];
double particle_radius[number_of_particles];
double particle_YS[number_of_particles], particle_E_star[number_of_particles], particle_G_star[number_of_particles];	
double particle_mass[number_of_particles], particle_d_mass[number_of_particles], particle_moment[number_of_particles];
double particle_velocity_x[number_of_particles], particle_velocity_y[number_of_particles], particle_velocity_z[number_of_particles];
double particle_force_x[number_of_particles], particle_force_y[number_of_particles], particle_force_z[number_of_particles];
double particle_torque_x[number_of_particles], particle_torque_y[number_of_particles], particle_torque_z[number_of_particles];
double particle_angular_velocity_x[number_of_particles], particle_angular_velocity_y[number_of_particles], particle_angular_velocity_z[number_of_particles];

double  wall_force = 0.0;
double  wall_mass = 100; //dimensionless, default MASS is glass so the wall particle mass in this case is 0.1 glass mass
double  flip_time=5.0;  //unit is second. 3 second for each trial

int     number_of_cells;
int     number_x, number_y, number_z;

double	bottom_wall_force, damping, friction_coefficient;
double 	TIME, MASS, LENGTH;
double  wall_speed;
double 	dt;
double 	elapsed_time, particle_diameter;

sim_material    material = GLASS;

int		incontact=-1;  //what does this value do?
int     cohesive = 0;  //dry particles for now
double   vinit = -8675309.0;	/* this is just a dummy number. it
					 * has no significance ... a number
					 * that can never be true */
double   S_crit, A, B, C;
double   surface_tension, viscosity;

void
usage()
{
	printf("\n");
	fprintf(stderr, "usage:  PD -i (name) [h] \n\
	lattice boltzmann code\n\
\n\
-i (file)\t input file (required)\n(format: x y z solid_node_boolean)\n\
\n OPTIONAL arguments\n\
-c 	 \t use wet (cohesive) particles    \n\
-h   \t show this message  \n\n");

}

int
main(int argc, char **argv)
{
	int m, number_of_files = 0;
	double running_time=30.0, write_time = 0.0;  //running_time=2 secobd. Maybe a test?
	char *infile = NULL;
	char outfile[256];
    
    double   stored_time=0.0;
	FILE           *fpout;
	
	/** I/O handling **/

	extern char    *optarg;
	int             option_letter = 0, usage_flag = 0;
    

	/** parse commmand line options **/

	while ((option_letter = getopt(argc, argv, "i:ch")) != -1) {

		if (option_letter == 'i')
        {
			infile = optarg;
			usage_flag++;
		}
        else if (option_letter == 'c')
        {
            cohesive = 1;
        }
        else if (option_letter == 'h')
        {
			usage();
			exit(1);
		}
	}

	if (usage_flag < 0) {
		usage();
		exit(1);
	}

	particle_diameter = 0.005;
	
	number_of_time_steps=0;
	
	initialize(infile);
	number_x = (int) (ceil) ((double) Nx / size_of_cells);
	number_y = (int) (ceil) ((double) Ny / size_of_cells);
	number_z = (int) (ceil) ((double) Nz / size_of_cells);

	number_of_cells = number_x * number_y * number_z;

	cell = malloc(number_of_cells * sizeof(List));

	get_cell_neighbors();			
	
    /* info for piston wall simulations */
    wall_speed = WALL_SPEED * TIME / LENGTH;
    printf("wall_speed is %f \n", wall_speed);
    printf("cohesive_value %d, wall_mass %f, flip_time %e \n", cohesive, wall_mass, flip_time);
   
	while (elapsed_time <= running_time * (1.0 / TIME)) {
        
		for (m=number_of_wall_particles;m<number_of_particles;m++) { //this part is initialize
            
			double px, py, pz;
            
            Gdirection=-1.0;  //switch gravity
            
            S_crit = pow(VOLUME, (1.0 / 3.0));	/* VOLUME is defined above -
                                                 * this is from Lian/Thornton */
            A = -1.1 * pow(VOLUME, -0.53);	/* we use a correlation for Fc from
                                             * Horio */
            B = -0.019 * log(VOLUME) + 0.48;	/* no need to correct for
                                                 * different */
            C = 0.0042 * log(VOLUME) + 0.078;	/* length scale since our S =
                                                 * 2h in their notation */
            
			px = particle_x[m] += move_particle(particle_mass[m], particle_velocity_x[m], particle_force_x[m]);  //px is the x coordinate of one particle [m]
			if (px < 0.0)
				particle_x[m] = px + (double)Nx;// Nx=number of cells in the x-direction, size of cell =1
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
            
			particle_velocity_x[m] += update_velocity(particle_mass[m], particle_force_x[m]);
			particle_velocity_y[m] += update_velocity(particle_mass[m], particle_force_y[m]);
			particle_velocity_z[m] += update_velocity(particle_mass[m], particle_force_z[m]);
            
			particle_angular_velocity_x[m] += update_velocity(particle_moment[m], particle_torque_x[m]);
			particle_angular_velocity_y[m] += update_velocity(particle_moment[m], particle_torque_y[m]);
			particle_angular_velocity_z[m] += update_velocity(particle_moment[m], particle_torque_z[m]);
            
			particle_force_x[m] = particle_force_z[m] = 0.0;
            
            if (particle_color[m]==7)
            {
                particle_force_y[m]= 0.0;  //ignore gravity for lighter particles  , in our case, acetate
            }
            else if (particle_color[m]==3)
            {
                particle_force_y[m] = particle_force_y_ng[m]*Gdirection;  //this is a new way to define gravity free condition. Not sure working yet
            }
        
			particle_torque_x[m] = particle_torque_y[m] = particle_torque_z[m] = 0.0;

            
		}

        
        handle_walls();
        /* info for wet simulations */
        
        surface_tension = SURFACE_TENSION / (MASS / (TIME * TIME));
        viscosity = VISCOSITY / (MASS / (LENGTH * TIME));

		/** find potential collisions/contacts between particles **/
		sort_particles();	/** first sort the particles into cells **/
        
		contact_detect();	/** then loop over cells and add particles to potential collision list **/  //This function also include contact_force calculation
        

		for (m=number_of_wall_particles ;m<number_of_particles;m++) {

			particle_velocity_x[m] += update_velocity(particle_mass[m], particle_force_x[m]);
			particle_velocity_y[m] += update_velocity(particle_mass[m], particle_force_y[m]);

            particle_velocity_z[m] += update_velocity(particle_mass[m], particle_force_z[m]);
            
			particle_angular_velocity_x[m] += update_velocity(particle_moment[m], particle_torque_x[m]);
			particle_angular_velocity_y[m] += update_velocity(particle_moment[m], particle_torque_y[m]);
			particle_angular_velocity_z[m] += update_velocity(particle_moment[m], particle_torque_z[m]);
			
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
				fprintf(fpout, "%e %e %e %e %e %e %e %e %e %e 0.0 %i\n", particle_x[m], particle_y[m], particle_z[m], particle_velocity_x[m], particle_velocity_y[m], particle_velocity_z[m], particle_angular_velocity_x[m], particle_angular_velocity_y[m],	particle_angular_velocity_z[m], particle_radius[m], particle_color[m]);
			
			fclose(fpout);
						
			number_of_files++;
            
            printf("gdirection %f  stored_time(s) %f wall_force %f \n ", Gdirection,stored_time*TIME, wall_force);
			write_time += (1.0 / 30.0) * (1.0 / TIME);
            //write_time += dt;
			
		}
			
		elapsed_time += dt;
		number_of_time_steps++;
        stored_time  += dt;

        if (stored_time>=flip_time/TIME)
        {
            stored_time=0.0;
            Gdirection=-Gdirection;
            printf("flip gravity with Gdirection %f (negative is down)\n", Gdirection);
        }
	}
	
	return(0);
} 
