#include "pd.h"

void
cross(double a_x, double a_y, double a_z, double b_x, double b_y, double b_z, double *c_x, double *c_y, double *c_z)
{
	*c_x = a_y * b_z - a_z * b_y;
	*c_y = -(a_x * b_z - a_z * b_x);
	*c_z = a_x * b_y - a_y * b_x;
}
