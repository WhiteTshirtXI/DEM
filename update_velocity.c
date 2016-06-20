#include "pd.h"

double 
update_velocity(double mass, double force)
{
	return (0.5 * dt * force / mass);
}
