#include "pd.h"

double 
move_particle(double mass, double velocity, double force)
{
	return (dt * velocity + 0.5 * dt * dt * force / mass);
}
