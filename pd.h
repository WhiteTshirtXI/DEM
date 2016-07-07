#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <math.h>
#include <malloc/malloc.h>
#include <sys/types.h>

#define	number_of_particles 7040
#define	number_of_wall_particles 1160
#define	number_of_moving_wall_particles 0

#define	number_of_walls 6
#define	number_of_units (number_of_particles+number_of_walls)

#define	number_of_interaction_pairs (number_of_units*(number_of_units+1)/2)

#define ASP		0.01	// asperity height for lubrication forces

#define XPERIODIC	1
#define YPERIODIC	0
#define	ZPERIODIC	1

#define  Nx 32    // number of cells in the x-direction
#define  Ny 50    // number of cells in the y-direction
#define  Nz 10     // number of cells in the z-direction

#define WALL_SPEED	0.0
#define WALL_MASS	0.0

#define	BETA			0.01
#define DENSITY_RATIO	6.0
#define FREQ			0.1			// in Hz
#define MAX_H			40.0
#define MAX_ANGLE		(9.0*M_PI/180.0)
#define	WIDTH			15.0

#define XOFFSET_LEFT	8.661	// initial position of the left wall 
								// (0.5*MAX_H*tan(MAX_ANLGLE))
#define XOFFSET_RIGHT	23.661	// initial position of the right wall 
								// (0.5*MAX_H*tan(MAX_ANLGLE)+WIDTH)

#define DEBUG	0
#define CELL_COORD_NUMBER 	(13)	/* (4+3*N_x/size_of_cells) */
#define size_of_cells	1.0

#define VW		0.0	// speed of wall (dimensionless?)
#define	PERIOD	6000.0
#define VWALL(elapsed_time)	(VW*sin(2.0*M_PI*elapsed_time/PERIOD))

#define right	0
#define left	1
#define back	2
#define front	3
#define top		4
#define bottom	5

typedef enum {
	ELASTIC,
	PLASTIC
} sim_contact_mechanics;

typedef enum
{
  	NONE,
	LOADING,
	UNLOADING
} sim_normal_contact_status;

typedef struct {
 int done, head;
} List;

#define pair(i, j) (j-i+(i*(i+1)/2))

extern double	accumulated_angle, angle_sign;
extern double 	dt, TIME, MASS, LENGTH, elapsed_time, particle_diameter;
extern double	bottom_wall_force, damping, friction_coefficient;
extern double   max_x, min_x, max_y, min_y, max_z, min_z;
extern double	vinit;
extern int		incontact;
extern long long	number_of_time_steps;

extern int      number_of_cells;
extern int      number_x, number_y, number_z;

extern sim_contact_mechanics contact_mechanics;

extern double particle_norm_old_x[number_of_interaction_pairs], particle_norm_old_y[number_of_interaction_pairs], particle_norm_old_z[number_of_interaction_pairs];
extern double particle_Rp[number_of_interaction_pairs], particle_disp_normal_max[number_of_interaction_pairs], particle_force_normal_max[number_of_interaction_pairs];
extern double particle_force_tan_old_mag[number_of_interaction_pairs], particle_force_tan_old_x[number_of_interaction_pairs], particle_force_tan_old_y[number_of_interaction_pairs], particle_force_tan_old_z[number_of_interaction_pairs];
extern double particle_disp_tan_old_mag[number_of_interaction_pairs], particle_disp_tan_old_x[number_of_interaction_pairs], particle_disp_tan_old_y[number_of_interaction_pairs], particle_disp_tan_old_z[number_of_interaction_pairs];
extern double particle_sign_of_disp[number_of_interaction_pairs];
extern long long	particle_number_of_time_steps[number_of_interaction_pairs], particle_initial_time_step[number_of_interaction_pairs];
extern int particle_nextx[number_of_interaction_pairs], particle_nexty[number_of_interaction_pairs], particle_nextz[number_of_interaction_pairs];

extern int particle_color[number_of_particles];
 
extern double particle_x[number_of_particles], particle_y[number_of_particles], particle_z[number_of_particles];
extern double particle_radius[number_of_particles];
extern double particle_YS[number_of_particles], particle_E_star[number_of_particles], particle_G_star[number_of_particles];	
extern double particle_mass[number_of_particles], particle_d_mass[number_of_particles], particle_moment[number_of_particles];
extern double particle_velocity_x[number_of_particles], particle_velocity_y[number_of_particles], particle_velocity_z[number_of_particles];
extern double particle_force_x[number_of_particles], particle_force_y[number_of_particles], particle_force_z[number_of_particles];
extern double particle_torque_x[number_of_particles], particle_torque_y[number_of_particles], particle_torque_z[number_of_particles];
extern double particle_angular_velocity_x[number_of_particles], particle_angular_velocity_y[number_of_particles], particle_angular_velocity_z[number_of_particles];

extern void		initialize(char *infile);
extern double move_particle(double mass, double velocity, double force);
extern double update_velocity(double mass, double force);
extern double wall_x(int n, double x, double y, double z);
extern double wall_y(int n, double x, double y, double z);
extern double wall_z(int n, double x, double y, double z);
extern double wall_vx(int n);
extern double wall_vy(int n);
extern double wall_vz(int n);
extern void contact_detect();
extern void calculate_contact_force(int i, int j);
extern void add_to_list(int i, int j);
extern void erase_contact_data(int q);
extern void handle_walls();
extern void handle_moving_walls();
extern void cross(double a_x, double a_y, double a_z, double b_x, double b_y, double b_z, double *c_x, double *c_y, double *c_z);


