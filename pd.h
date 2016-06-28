#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <math.h>
#include <malloc/malloc.h>
#include <sys/types.h>

#define	number_of_particles 3552
#define	number_of_wall_particles  352
#define number_of_moving_wall_particles 176
#define	number_of_units (number_of_particles)

#define SURFACE_TENSION 0.072
#define VISCOSITY 0.001


#define	number_of_interaction_pairs (number_of_units*(number_of_units+1)/2)

#define XPERIODIC	1
#define YPERIODIC	1     //change to periodic, no concept solid wall
#define	ZPERIODIC	1

#define  Nx 20    // number of cells in the x-direction
#define  Ny 40    // number of cells in the y-direction
#define  Nz 8     // number of cells in the z-direction

#define DEBUG	0
#define CELL_COORD_NUMBER 	(13)	/* (4+3*simulation_size_x/size_of_cells) */
#define size_of_cells	1.0

#define ASP 0.001				/* 1.0x10^-6/dp */  //asperity height for lubrication forces
#define VOLUME 0.001		/* dimensionless volume of a *single* bridge */
/* vol_percent_liq*[(4/3)*PI*rp^3]/[solid_frac*5 bridges per particle] */
/* since we use dp to make our equations dimensionless */
/* vol_percent_liq*[(1/3)*PI]/[solid_frac*10] because of the factor of 8 */

#define WALL_SPEED	1  //m/s

#define VW		0.0	// speed of wall (dimensionless?)
#define	PERIOD	6000.0
#define VWALL(elapsed_time)	(VW*sin(2.0*M_PI*elapsed_time/PERIOD))  //what is this?????? Something related to a circle?

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
 int done, head, active, toggled;
 int x, y, z, begin_of_list, neighbors[CELL_COORD_NUMBER], anti_neighbors[CELL_COORD_NUMBER];
} List;

typedef enum
{
    GLASS,
    STEEL,
    ALUMINUM,
    ACRYLIC,
    ACETATE,
    SOFT,
    DUKE
} sim_material;

#define pair(i, j) (j-i+(i*(i+1)/2))

double Gdirection;

extern double S_crit, A, B, C;
extern double   dt;
extern double 	TIME, MASS, LENGTH, elapsed_time, particle_diameter;
extern double	bottom_wall_force, damping, friction_coefficient;
extern double	vinit;
extern int		incontact;
extern long long	number_of_time_steps;
extern sim_material material;
extern int      cohesive;

extern int      number_of_cells;
extern int      number_x, number_y, number_z;
extern List    *cell;

extern int particle_touching[number_of_interaction_pairs];

extern sim_normal_contact_status particle_force_normal_status[number_of_interaction_pairs];
extern sim_contact_mechanics contact_mechanics;

extern double particle_Rp[number_of_interaction_pairs], particle_force_normal_old[number_of_interaction_pairs], particle_force_normal_max[number_of_interaction_pairs];
extern double particle_force_tan_old_x[number_of_interaction_pairs], particle_force_tan_old_y[number_of_interaction_pairs], particle_force_tan_old_z[number_of_interaction_pairs];
extern double particle_contact_radius[number_of_interaction_pairs];
extern long long	particle_number_of_time_steps[number_of_interaction_pairs];

extern int particle_color[number_of_particles];
extern int particle_nextx[number_of_particles], particle_nexty[number_of_particles], particle_nextz[number_of_particles];	
extern int particle_list[number_of_particles], particle_cell_id[number_of_particles];
extern int particle_coordination_number_old[number_of_particles], particle_coordination_number_new[number_of_particles];	
 
extern double particle_x[number_of_particles], particle_y[number_of_particles], particle_z[number_of_particles];
extern double particle_radius[number_of_particles];
extern double particle_YS[number_of_particles], particle_E_star[number_of_particles], particle_G_star[number_of_particles];	
extern double particle_mass[number_of_particles], particle_d_mass[number_of_particles], particle_moment[number_of_particles];
extern double particle_velocity_x[number_of_particles], particle_velocity_y[number_of_particles], particle_velocity_z[number_of_particles];
extern double particle_force_x[number_of_particles], particle_force_y[number_of_particles], particle_force_z[number_of_particles];
extern double particle_force_y_ng[number_of_particles];  //extra line for no gravity condition
extern double particle_torque_x[number_of_particles], particle_torque_y[number_of_particles], particle_torque_z[number_of_particles];
extern double particle_angular_velocity_x[number_of_particles], particle_angular_velocity_y[number_of_particles], particle_angular_velocity_z[number_of_particles];

extern void		initialize(char *infile);
extern double move_particle(double mass, double velocity, double force);
extern double update_velocity(double mass, double force);
extern int sort_particles();
extern void get_cell_neighbors();
extern double wall_x(int n, double x, double y, double z);
extern double wall_y(int n, double x, double y, double z);
extern double wall_z(int n, double x, double y, double z);
extern double wall_vx(int n);
extern double wall_vy(int n);
extern double wall_vz(int n);
extern void handle_walls();
extern void handle_moving_walls();
extern void contact_detect();
extern void calculate_contact_force(int i, int j);
extern void add_to_list(int i, int j);
extern void erase_contact_data(int q);
extern void  get_material_properties(sim_material this_material, double *thermal_expansion_coef, double *youngs_modulus, double *density, double *poisson_ratio, double *yield_stress, double *conductivity, double *heat_capacity);

extern double surface_tension, viscosity;
extern double wall_force;
extern double wall_mass;
extern double wall_speed;
extern int calculate_wet_force;


//Rest is for Contact Index files
extern FILE           *fpout_Nx;
extern char outfile_Nx[256];
extern FILE           *fpout_Nx_index;
extern char outfile_Nx_index[256];
extern int  Nx_index;
extern int  W_Nx_index;