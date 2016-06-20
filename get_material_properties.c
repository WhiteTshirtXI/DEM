/*
 * This file is part of freeDEM.
 * 
 * Copyright (C) 1998-2004 Joseph McCarthy <jjmcc@pitt.edu>
 * 
 * freeDEM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 2 of the License, or (at your option) any later
 * version.
 * 
 * freeDEM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 * 
 * You should have received a copy of the GNU General Public License along with
 * freeDEM; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 * 
 */

#include "pd.h"

void 
get_material_properties(sim_material this_material, double *thermal_expansion_coef, double *youngs_modulus, double *density, double *poisson_ratio, double *yield_stress, double *conductivity, double *heat_capacity)
{

	/** default values **/

	*youngs_modulus = 68.95e9;	/* Pa=kg/m s^2		 */
	*density = 2700.0;	/* kg/m^3					 */
	friction_coefficient = 0.30;
	*poisson_ratio = 0.33;
	*yield_stress = 68.95e7;	/* Pa = kg/m s^2	 */
	*conductivity = 1.7;	/* W/mK=kg/K s^3	 */
	*heat_capacity = 800.0;	/* W/kgK=m/K s^3	 */
	*thermal_expansion_coef = 17.5e-6;	/* [ 1/K] */

	/** several materials **/
	if (this_material == GLASS) {
		*youngs_modulus = 68.95e9;	/* Pa=kg/m s^2		 */
		*density = 2700.0;	/* kg/m^3					 */
		friction_coefficient = 0.30;
		*poisson_ratio = 0.33;
		*yield_stress = 68.95e7;	/* Pa = kg/m s^2	 */
		*conductivity = 1.7;	/* W/mK=kg/K s^3	 */
		*heat_capacity = 800.0;	/* W/kgK=m/K s^3	 */
	} else if (this_material == STEEL) {
		*youngs_modulus = 193.0e9;	/* Pa=kg/m s^2		 */
		*density = 7900.0;	/* kg/m^3					 */
		friction_coefficient = 0.30;
		*poisson_ratio = 0.29;
		*yield_stress = 265.0e7;	/* Pa = kg/m s^2	 */
		*conductivity = 15.0;	/* W/mK=kg/K s^3	 */
		*heat_capacity = 477.0;	/* W/kgK=m/K s^3	 */
	} else if (this_material == ALUMINUM) {
		*youngs_modulus = 110.0e9;	/* Pa=kg/m s^2		 */
		*density = 2700.0;	/* kg/m^3					 */
		friction_coefficient = 0.30;
		*poisson_ratio = 0.33;
		*yield_stress = 320.0e6;	/* Pa = kg/m s^2	 */
		*conductivity = 180.0;	/* W/mK=kg/K s^3	 */
		*heat_capacity = 900.0;	/* W/kgK=m/K s^3	 */
	} else if (this_material == ACRYLIC) {
		*youngs_modulus = 2.9e9;	/* Pa=kg/m s^2		 */
		*density = 1350.0;	/* kg/m^3					 */
		friction_coefficient = 0.30;
		*poisson_ratio = 0.43;
		*yield_stress = 45.0e6;	/* Pa = kg/m s^2	- mostly a
					 * guess */
		*conductivity = 0.14;	/* W/mK=kg/K s^3	 */
		*heat_capacity = 1200.0;	/* W/kgK=m/K s^3	 */
	} else if (this_material == ACETATE) {
		*youngs_modulus = 1.5e9;	/* Pa=kg/m s^2		 */
		*density = 1300.0;	/* kg/m^3					 */
		friction_coefficient = 0.30;
		*poisson_ratio = 0.43;
		*yield_stress = 30.0e6;	/* Pa = kg/m s^2	- mostly a
					 * guess */
		*conductivity = 0.20;	/* W/mK=kg/K s^3	 */
		*heat_capacity = 1400.0;	/* W/kgK=m/K s^3	 */
	} else if (this_material == SOFT) {
		*youngs_modulus = 3.0e7;	/* Pa=kg/m s^2		 */
		*density = 1000.0;	/* kg/m^3					 */
		friction_coefficient = 0.30;
		*poisson_ratio = 0.33;
		*yield_stress = 3.0e5;	/* Pa = kg/m s^2	- mostly a
					 * guess */
		*conductivity = 30.0;	/* W/mK=kg/K s^3	 */
		*heat_capacity = 800.0;	/* W/kgK=m/K s^3	 */
	}
	
}
