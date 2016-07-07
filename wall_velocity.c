#include "pd.h"

double 
wall_vx(int n)
{
	if (n==(number_of_particles+top)) return(VWALL(elapsed_time));
	else return (0.0);
}

double 
wall_vy(int n)
{
	return (0.0);
}

double 
wall_vz(int n)
{
	return (0.0);
}

