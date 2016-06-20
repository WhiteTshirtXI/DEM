#include "pd.h"

double 
wall_x(int n, double x, double y, double z) /* returns the x position of each smooth wall */
{
	/* x position of the left wall */
        if (n==(number_of_particles+left))
			return(1.5); 
	/* x position of the right wall */
        else if (n==(number_of_particles+right))
			return(Nx-1.5);
		else return(x);
}

double 
wall_y(int n, double x, double y, double z) /* returns the x position of each smooth wall */
{
	/* y position of the bottom wall */
        if (n==(number_of_particles+bottom))
			return(1.5); 
	/* y position of the top wall */
        else if (n==(number_of_particles+top))
			return(Ny-1.5);
		else return(y);
}

double 
wall_z(int n, double x, double y, double z) /* returns the z position of each smooth wall */
{
	/* z position of the front wall */
        if (n==(number_of_particles+back))
			return(1.5); 
	/* z position of the back wall */
        else if (n==(number_of_particles+front))
			return(Nz-1.5);
		else return(z);
}
