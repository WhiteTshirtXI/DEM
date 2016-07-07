#include "pd.h"

double 
update_velocity(double mass, double force)
{
	return (dt * force / mass);
}
