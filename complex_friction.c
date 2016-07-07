/*
 * This file is part of freeDEM. Copyright (C) 1998-2004 Joseph McCarthy
 * <jjmcc@pitt.edu>
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
 */

#include "pd.h"

void
calculate_contact_force(int touching, int l, int i, int j)
{
	sim_normal_contact_status fn_status;
	int             the1, correction_index;

	double          kt, CDF = 1.0, CTD, disp_crit, UFL = 1.0;

	double          norm_x, norm_y, norm_z, xij, yij, zij, vxij, vyij, vzij,
	                angvi_x, angvi_y, angvi_z, angvj_x, angvj_y, angvj_z;
	double          alpha, separation, separation_sq, radii_sum;
	double          E_eff, YS, R_eff, F_yield, Kn, Kne, G_eff;
	double          theta_base = 0.0, fs_x = 0.0, fs_y = 0.0, fs_z = 0.0,
	                fs_mag = 0.0, fs_old_mag = 0.0, tx = 0.0, ty = 0.0,
	                tz = 0.0;
	double          Rp, fn, fn_old, fn_max = 0.0, fn_test;
	double          disp_x = 0.0, disp_y = 0.0, disp_z = 0.0, fx = 0.0, fy = 0.0,
	                fz = 0.0, vn, rparti, rpartj;
	double          disp_mag, disp_abs, d_disp_x, d_disp_y, d_disp_z,
	                d_disp_mag, theta = 1.0;

	double          T_star, T_star_star, disp_x_old, disp_y_old, disp_z_old, disp_mag_old;
	double          norm_x_old, norm_y_old, norm_z_old, v_disp_x, v_disp_y,
	                v_disp_z;
	double          td_x, td_y, td_z, fz_x, fz_y, fz_z, f_x, f_y, f_z,
	                thdr_x, thdr_y, thdr_z, cs_x, cs_y, cs_z, trs_x,
	                trs_y, trs_z;

	double          contact_radius, dalpha;
	Particle        particle_i = particle[i];
	Particle        particle_j = particle[j];

	rparti = particle_i.radius;
	rpartj = particle_j.radius;

	radii_sum = rparti + rpartj;

	xij = particle_i.x - particle_j.x;
	if (xij > 0.5 * N_x)
		xij -= N_x;
	else if (xij < -0.5 * N_x)
		xij += N_x;

	yij = particle_i.y - particle_j.y;
	if (yij > 0.5 * N_y)
		yij -= N_y;
	else if (yij < -0.5 * N_y)
		yij += N_y;

	zij = particle_i.z - particle_j.z;
	if (zij > 0.5 * N_z)
		zij -= N_z;
	else if (zij < -0.5 * N_z)
		zij += N_z;

	separation_sq = xij * xij + yij * yij + zij * zij;
	radii_sum = rparti + rpartj;

	separation = sqrt(xij * xij + yij * yij + zij * zij);
	alpha = radii_sum - separation;

	norm_x = xij / separation;
	norm_y = yij / separation;
	norm_z = zij / separation;

	Rp = particle_i.contact_data_new[l].Rp;
	fn_old = particle_i.contact_data_new[l].force_normal_old;
	fn_max = particle_i.contact_data_new[l].force_normal_max;
	fn_status = particle_i.contact_data_new[l].force_normal_status;

	CTD = particle_i.contact_data_new[l].CTD;

	T_star = particle_i.contact_data_new[l].T_star;
	T_star_star = particle_i.contact_data_new[l].T_star_star;

	norm_x_old = particle_i.contact_data_new[l].norm_x_old;
	norm_y_old = particle_i.contact_data_new[l].norm_y_old;
	norm_z_old = particle_i.contact_data_new[l].norm_z_old;

	disp_x_old = particle_i.contact_data_new[l].disp_x_old;
	disp_y_old = particle_i.contact_data_new[l].disp_y_old;
	disp_z_old = particle_i.contact_data_new[l].disp_z_old;

	disp_mag_old = particle_i.contact_data_new[l].disp_mag_old;

	disp_abs = particle_i.contact_data_new[l].disp_abs;

	fs_old_mag = particle_i.contact_data_new[l].fs_old_mag;

	angvi_x = particle_i.angular_velocity_x;
	angvi_y = particle_i.angular_velocity_y;
	angvi_z = particle_i.angular_velocity_z;

	angvj_x = particle_j.angular_velocity_x;
	angvj_y = particle_j.angular_velocity_y;
	angvj_z = particle_j.angular_velocity_z;

	vxij = particle_i.velocity_x - particle_j.velocity_x;
	vyij = particle_i.velocity_y - particle_j.velocity_y;
	vzij = particle_i.velocity_z - particle_j.velocity_z;

	vn = norm_x * vxij + norm_y * vyij + norm_z * vzij;

	dalpha = -vn * dt;

	YS = particle_i.YS;
	if (YS > particle_j.YS)
		YS = particle_j.YS;

	E_eff = particle_i.E_star * particle_j.E_star / (particle_i.E_star + particle_j.E_star);
	G_eff = particle_i.G_star * particle_j.G_star / (particle_i.G_star + particle_j.G_star);
	R_eff = 2.0 * rparti * rpartj / (rparti + rpartj);

	F_yield = YS * YS * YS * M_PI * M_PI * M_PI * R_eff * R_eff / (6.0 * E_eff * E_eff);

	/* for loading */

	if (vn <= 0.0) {

		if (fn_status == NONE || fn_old >= fn_max)
			Rp = R_eff;
		fn_status = LOADING;

		if (fn_old < fn_max && fn_old != 0.0)
			contact_radius = particle_i.contact_data_new[l].contact_radius = pow((3.0 * Rp * fn_old / (4.0 * E_eff)), 0.333);
		else
			contact_radius = particle_i.contact_data_new[l].contact_radius = sqrt(Rp * alpha);

		Kn = Kne = 2.0 * E_eff * contact_radius;

		/* elastic loading force */

		fn = fn_old + Kn * dalpha;

		/* if we have visco-damping */

		if (contact_mechanics == ELASTIC)
			fn -= damping * vn * alpha;

		/* if we have plastic deformation */

		else if (contact_mechanics == PLASTIC)
			if (fn > F_yield) {
				contact_radius = particle_i.contact_data_new[l].contact_radius = sqrt((2.0 * fn_old + F_yield) / (2.0 * M_PI * YS));
				Kn = M_PI * R_eff * YS;
				fn = fn_old + Kn * dalpha;
			}
		/* remember this for some contact models */

		if (fn > fn_max)
			fn_max = fn;

	}
	/* for unloading */

	else {

		/* if we have plastic deformation */

		if (contact_mechanics == PLASTIC) {
			if (fn_status == LOADING) {
				fn_status = UNLOADING;
				contact_radius = particle_i.contact_data_new[l].contact_radius = sqrt((2.0 * fn_old + F_yield) / (2.0 * M_PI * YS));
				Rp = 4.0 * E_eff * contact_radius * contact_radius * contact_radius / (3.0 * fn_max);
			}
		}
		contact_radius = particle_i.contact_data_new[l].contact_radius = pow((3.0 * Rp * fn_old / (4.0 * E_eff)), 0.333);
		Kn = 2.0 * E_eff * contact_radius;
		fn = fn_old + Kn * dalpha;

		/* if we have visco-damping */

		if (contact_mechanics == ELASTIC)
			fn -= damping * vn * alpha;

	}

	if (fn < 0.0)
		fn = 0.0;

	/* start tangential force */

	/*
	 * calculate the translation portion of the tangential displacement
	 * as n X v X n
	 */

	td_x = ((norm_z * vxij - norm_x * vzij) * norm_z - (norm_x * vyij - norm_y * vxij) * norm_y) * dt;
	td_y = ((norm_x * vyij - norm_y * vxij) * norm_x - (norm_y * vzij - norm_z * vyij) * norm_z) * dt;
	td_z = ((norm_y * vzij - norm_z * vyij) * norm_y - (norm_z * vxij - norm_x * vzij) * norm_x) * dt;

	/*
	 * since the particle contact has potentially rotated, the new
	 * tangential displacement will not necessarily align properly with
	 * the old, even if they are coincident in the particle frame of
	 * reference since we "remember" the old displacement in the system
	 * frame; therefore we need to use the old and new normal vectors to
	 * "rotate" the old tangential displacement (and force) direction so
	 * that it is still correct in the particle-contact frame
	 */

	/*
	 * first, we check to make sure there IS an old normal vector (if
	 * not, skip this)
	 */

	if (norm_x_old != 0.0 || norm_y_old != 0.0 || norm_z_old != 0.0) {
		double          n_dot_n, n_X_n_x, n_X_n_y, n_X_n_z;
		double          A1, A2, A3, theta;	/* trig functions in
							 * rotation matrix */

		/*
		 * In order to calculate the angle of rotation, we need the
		 * dot product between n's
		 */

		n_dot_n = norm_x * norm_x_old + norm_y * norm_y_old + norm_z * norm_z_old;

		/*
		 * The axis of the rotation comes from the cross product
		 * between n's
		 */

		cross(norm_x, norm_y, norm_z, norm_x_old, norm_y_old, norm_z_old, &n_X_n_x, &n_X_n_y, &n_X_n_z);

		/* We now calculate some trig functions of the rotation angle */

		if (n_dot_n > 1.0)
			theta = 0.0;
		else
			theta = acos(n_dot_n);

		A3 = cos(theta);
		A1 = 1.0 - A3;
		A2 = sin(theta);

		/*
		 * We need to save the old components of disp_old so that
		 * they all change "at once"
		 */

		fz_x = disp_x_old;
		fz_y = disp_y_old;
		fz_z = disp_z_old;

		/*
		 * now we use the rotation matrix to move the old
		 * displacement
		 */

		disp_x_old = (n_X_n_x * n_X_n_x * A1 + A3) * fz_x + (n_X_n_x * n_X_n_y * A1 - n_X_n_z * A2) * fz_y + (n_X_n_x * n_X_n_z * A1 + n_X_n_y * A2) * fz_z;
		disp_y_old = (n_X_n_x * n_X_n_y * A1 + n_X_n_z * A2) * fz_x + (n_X_n_y * n_X_n_y * A1 + A3) * fz_y + (n_X_n_y * n_X_n_z * A1 - n_X_n_x * A2) * fz_z;
		disp_z_old = (n_X_n_x * n_X_n_z * A1 + n_X_n_y * A2) * fz_x + (n_X_n_z * n_X_n_y * A1 + n_X_n_x * A2) * fz_y + (n_X_n_z * n_X_n_z * A1 + A3) * fz_z;

	}
	/* calculate the total relative (rotational) velocity/speed */

	trs_x = angvi_x * rparti + angvj_x * rpartj;
	trs_y = angvi_y * rparti + angvj_y * rpartj;
	trs_z = angvi_z * rparti + angvj_z * rpartj;

	/*
	 * calculate the components of the tangential velocity by taking the
	 * cross product of rotational velocity/speed with the normal vector
	 */

	cross(trs_x, trs_y, trs_z, norm_x, norm_y, norm_z, &v_disp_x, &v_disp_y, &v_disp_z);

	/*
	 * calculate the total tangential displacement by combining the
	 * translational part with the rotational part
	 */

	d_disp_x = td_x - v_disp_x * dt;
	d_disp_y = td_y - v_disp_y * dt;
	d_disp_z = td_z - v_disp_z * dt;

	/*
	 * calculate the magnitude of the corrected tangential displacement
	 * in each direction
	 */

	disp_x = disp_x_old + d_disp_x;
	disp_y = disp_y_old + d_disp_y;
	disp_z = disp_z_old + d_disp_z;

	/*
	 * determine whether the new increment of displacement is in the same
	 * or opposite direction of the total tangential displacement, and
	 * correct the old magnitude
	 */

	if ((disp_x * disp_x_old + disp_y * disp_y_old + disp_z * disp_z_old) < 0.0)
		CTD = -1.0 * CTD;

	disp_mag = CTD * sqrt(disp_x * disp_x + disp_y * disp_y + disp_z * disp_z);
	d_disp_mag = disp_mag - disp_mag_old;

	/* determine if the tangential force is increasing or decreasing */

	if (fs_old_mag < 0.0)
		CDF = -1.0;
	else
		CDF = 1.0;

	disp_crit = friction_coefficient * (fn - fn_old) / (8.0 * G_eff * contact_radius);

	if (fn >= fn_old) {
		if ((d_disp_mag * CDF > 0.0) && (d_disp_mag * CDF >= disp_crit + disp_abs) && fabs(fs_old_mag) < fabs(T_star))
			the1 = 3;
		else if ((d_disp_mag * CDF > 0.0) && (d_disp_mag * CDF >= disp_crit + disp_abs))
			the1 = 1;
		else if (d_disp_mag * CDF > 0.0)
			the1 = 4;
		else if ((fabs(d_disp_mag) >= disp_crit + disp_abs) && (fs_old_mag * CDF > T_star_star * CDF) && (T_star_star != 0.0))
			the1 = 3;
		else if (fabs(d_disp_mag) >= disp_crit + disp_abs) {
			the1 = 2;
			UFL = -1.0;
		} else {
			the1 = 4;
			if (T_star_star == 0.0)
				UFL = -1.0;
		}
	}
	if (fn < fn_old) {
		if ((d_disp_mag * CDF > 0.0) && fabs(fs_old_mag) < fabs(T_star))
			the1 = 3;
		else if (d_disp_mag * CDF > 0.0)
			the1 = 1;
		else if ((fs_old_mag * CDF > T_star_star * CDF) && (T_star_star != 0.0))
			the1 = 3;
		else {
			the1 = 2;
			UFL = -1.0;
		}
	}
	fn_test = friction_coefficient * fn;

	if (fn_test != 0.0) {
		if (the1 == 1)
			theta_base = (1.0 - (fs_old_mag * CDF + friction_coefficient * (fn - fn_old)) / fn_test);
		if (the1 == 2)
			theta_base = (1.0 - (T_star * CDF - fs_old_mag * CDF + 2.0 * friction_coefficient * (fn - fn_old)) / (2.0 * fn_test));
		if (the1 == 3)
			theta_base = (1.0 - (fs_old_mag * CDF - T_star_star * CDF + 2.0 * friction_coefficient * (fn - fn_old)) / (2.0 * fn_test));
		if (the1 == 4)
			theta_base = 1.0;

		if (theta_base < 0.0)
			theta_base = 0.0;

		theta = pow(theta_base, 1.0 / 3.0);
		fs_mag = fs_old_mag + CDF * (8.0 * G_eff * contact_radius * d_disp_mag * CDF * theta + friction_coefficient * (fn - fn_old) * UFL * (1.0 - theta));
		kt = CDF * (8.0 * G_eff * contact_radius * CDF * theta + friction_coefficient * (fn - fn_old) * UFL * (1.0 - theta) / d_disp_mag);
	} else {
		fs_mag = 0.0;
		kt = 0.0;
	}


	if (disp_mag != 0.0) {
		fs_x = -fs_mag * disp_x / disp_mag;
		fs_y = -fs_mag * disp_y / disp_mag;
		fs_z = -fs_mag * disp_z / disp_mag;
	} else
		fs_mag = fs_x = fs_y = fs_z = 0.0;

	/* update disp_abs, corresponds to DD in trubal */
	if ((fn - fn_old) >= 0.0) {
		if (fabs(fs_mag - fs_old_mag) != friction_coefficient * (fn - fn_old))
			disp_abs -= fabs(d_disp_mag) - disp_crit;
	} else
		disp_abs += disp_crit;
	if (disp_abs < 0.0)
		disp_abs = 0.0;

	/* find T* and T** ... */
	if ((d_disp_mag * CDF > 0.0) && (T_star != 0.0) && (T_star_star == 0.0))
		T_star_star = fs_old_mag;
	if ((d_disp_mag * CDF < 0.0) && (T_star == 0.0) && (T_star_star == 0.0))
		T_star = fs_old_mag;

	if (fs_mag != 0.0) {
		if (fn_test < fabs(fs_mag)) {
			fs_x = fs_x * fn_test / fabs(fs_mag);
			fs_y = fs_y * fn_test / fabs(fs_mag);
			fs_z = fs_z * fn_test / fabs(fs_mag);
			fs_mag = fn_test * fs_mag / fabs(fs_mag);
		}
	} else
		fs_mag = fs_x = fs_y = fs_z = 0.0;

	if ((T_star != 0.0) && (T_star != fs_old_mag))
		T_star += CDF * friction_coefficient * (fn - fn_old);
	if ((T_star_star != 0.0) && (T_star_star != fs_old_mag))
		T_star_star -= CDF * friction_coefficient * (fn - fn_old);

	if ((T_star != 0.0) && fabs(fs_mag) >= fabs(T_star)) {
		T_star = 0.0;
		T_star_star = 0.0;
	}
	if ((T_star_star != 0.0) && (fs_mag * CDF < T_star_star * CDF))
		T_star_star = 0.0;

	tx = (norm_y * fs_z - norm_z * fs_y);
	ty = (norm_z * fs_x - norm_x * fs_z);
	tz = (norm_x * fs_y - norm_y * fs_x);

	fx = norm_x * fn + fs_x;
	fy = norm_y * fn + fs_y;
	fz = norm_z * fn + fs_z;

	particle_i.torque_x -= rparti * tx;
	particle_i.torque_y -= rparti * ty;
	particle_i.torque_z -= rparti * tz;

	particle_i.force_x += fx;
	particle_i.force_y += fy;
	particle_i.force_z += fz;

	particle_i.contact_data_new[l].touching = 1;
	particle_i.contact_data_new[l].neighbor_id = j;
	particle_i.contact_data_new[l].Rp = Rp;
	particle_i.contact_data_new[l].force_normal_old = fn;
	particle_i.contact_data_new[l].force_normal_max = fn_max;
	particle_i.contact_data_new[l].force_normal_status = fn_status;

	particle_i.contact_data_new[l].CTD = CTD;
	particle_i.contact_data_new[l].T_star = T_star;
	particle_i.contact_data_new[l].T_star_star = T_star_star;
	particle_i.contact_data_new[l].disp_x_old = disp_x;
	particle_i.contact_data_new[l].disp_y_old = disp_y;
	particle_i.contact_data_new[l].disp_z_old = disp_z;
	particle_i.contact_data_new[l].norm_x_old = norm_x;
	particle_i.contact_data_new[l].norm_y_old = norm_y;
	particle_i.contact_data_new[l].norm_z_old = norm_z;
	particle_i.contact_data_new[l].fs_old_mag = fs_mag;
	particle_i.contact_data_new[l].disp_mag_old = disp_mag;
	particle_i.contact_data_new[l].disp_abs = disp_abs;

	particle_j.torque_x -= rpartj * tx;
	particle_j.torque_y -= rpartj * ty;
	particle_j.torque_z -= rpartj * tz;

	particle_j.force_x -= fx;
	particle_j.force_y -= fy;
	particle_j.force_z -= fz;

	particle[i] = particle_i;
	particle[j] = particle_j;

	if ((type == COUETTE || type == PISTON) && j < number_of_moving_wall_particles)
		wall_force -= fy;

}
