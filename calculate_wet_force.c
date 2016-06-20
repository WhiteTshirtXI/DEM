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

int 
calculate_wet_force(int i, int j)
{

	double          norm_x, norm_y, norm_z, xij, yij, zij, vxij, vyij,
	                vzij, angvi_x, angvi_y, angvi_z, angvj_x, angvj_y,
	                angvj_z;
	double          separation, separation_sq, radii_sum;
	double          wet_radii_sum, wet_radii_sum_sq, S;
	double          fstn, fvn, fluid_fn, fluid_kt;
	double          R_eff, fs_x, fs_y, fs_z;
	double          vdisp_x = 0.0, vdisp_y = 0.0, vdisp_z = 0.0;
	double          tx = 0.0, ty = 0.0, tz = 0.0, fx = 0.0, fy = 0.0,
	                fz = 0.0, vn, rparti, rpartj;

	rparti = particle[i].radius;
	
	angvi_x = particle[i].angular_velocity_x;
	angvi_y = particle[i].angular_velocity_y;
	angvi_z = particle[i].angular_velocity_z;

	if  (j>=number_of_particles) {
	
		rpartj = 0.0;
		radii_sum = rparti + rpartj;

		xij = particle[i].x - wall_x(j, elapsed_time, particle[i].x, particle[i].y, particle[i].z, particle[i].radius);
		yij = particle[i].y - wall_y(j, elapsed_time, particle[i].x, particle[i].y, particle[i].z, particle[i].radius);
		zij = particle[i].z - wall_z(j, elapsed_time, particle[i].x, particle[i].y, particle[i].z, particle[i].radius);
		
		angvj_x = angvj_y = angvj_z = 0.0;

		vxij = particle[i].velocity_x - wall_vx(j, elapsed_time);
		vyij = particle[i].velocity_y - wall_vy(j, elapsed_time);
		vzij = particle[i].velocity_z - wall_vz(j, elapsed_time);
		
		R_eff = rparti;
		
	}
	
	else { /* for all "regular" particles */

		rpartj = particle[j].radius;

		radii_sum = rparti + rpartj;

		xij = particle[i].x - particle[j].x;
		if (xij > 0.5 * simulation_size_x)
			xij -= simulation_size_x;
		else if (xij < -0.5 * simulation_size_x)
			xij += simulation_size_x;

		yij = particle[i].y - particle[j].y;
		if (yij > 0.5 * simulation_size_y)
			yij -= simulation_size_y;
		else if (yij < -0.5 * simulation_size_y)
			yij += simulation_size_y;

		zij = particle[i].z - particle[j].z;
		if (zij > 0.5 * simulation_size_z)
			zij -= simulation_size_z;
		else if (zij < -0.5 * simulation_size_z)
			zij += simulation_size_z;

		angvj_x = particle[j].angular_velocity_x;
		angvj_y = particle[j].angular_velocity_y;
		angvj_z = particle[j].angular_velocity_z;

		vxij = particle[i].velocity_x - particle[j].velocity_x;
		vyij = particle[i].velocity_y - particle[j].velocity_y;
		vzij = particle[i].velocity_z - particle[j].velocity_z;
		
		R_eff = rparti * rpartj / (rparti + rpartj);

	}

	separation_sq = xij * xij + yij * yij + zij * zij;

	wet_radii_sum = radii_sum + S_crit;
	wet_radii_sum_sq = wet_radii_sum * wet_radii_sum;

	if (wet_radii_sum_sq > separation_sq) {

		separation = sqrt(separation_sq);

		S = separation - radii_sum;
		if (S < ASP)
			S = ASP;

		norm_x = xij / separation;
		norm_y = yij / separation;
		norm_z = zij / separation;

		vn = norm_x * vxij + norm_y * vyij + norm_z * vzij;

		fstn = -M_PI * R_eff * surface_tension * (exp(A * S + B) + C);
		fvn = -6.0 * M_PI * R_eff * viscosity * vn * (R_eff / S);
		fluid_fn = fstn + fvn;

		vdisp_x = (norm_z * vxij - norm_x * vzij) * norm_z - (norm_x * vyij - norm_y * vxij) * norm_y - rparti * (angvi_y * norm_z - angvi_z * norm_y) - rpartj * (angvj_y * norm_z - angvj_z * norm_y);

		vdisp_y = (norm_x * vyij - norm_y * vxij) * norm_x - (norm_y * vzij - norm_z * vyij) * norm_z - rparti * (angvi_z * norm_x - angvi_x * norm_z) - rpartj * (angvj_z * norm_x - angvj_x * norm_z);

		vdisp_z = (norm_y * vzij - norm_z * vyij) * norm_y - (norm_z * vxij - norm_x * vzij) * norm_x - rparti * (angvi_x * norm_y - angvi_y * norm_x) - rpartj * (angvj_x * norm_y - angvj_y * norm_x);

		fluid_kt = 6.0 * M_PI * R_eff * viscosity * (0.53333 * log(R_eff / S) + 0.9588);

		fs_x = -fluid_kt * vdisp_x;
		fs_y = -fluid_kt * vdisp_y;
		fs_z = -fluid_kt * vdisp_z;

		tx = (norm_y * fs_z - norm_z * fs_y);
		ty = (norm_z * fs_x - norm_x * fs_z);
		tz = (norm_x * fs_y - norm_y * fs_x);

		fx = norm_x * fluid_fn + fs_x;
		fy = norm_y * fluid_fn + fs_y;
		fz = norm_z * fluid_fn + fs_z;

		particle[i].torque_x -= rparti * tx;
		particle[i].torque_y -= rparti * ty;
		particle[i].torque_z -= rparti * tz;

		particle[i].force_x += fx;
		particle[i].force_y += fy;
		particle[i].force_z += fz;

		if (j < number_of_particles) {
			particle[j].torque_x -= rpartj * tx;
			particle[j].torque_y -= rpartj * ty;
			particle[j].torque_z -= rpartj * tz;

			particle[j].force_x -= fx;
			particle[j].force_y -= fy;
			particle[j].force_z -= fz;
		}
		
		return (1);

	} else
		return (0);

}
