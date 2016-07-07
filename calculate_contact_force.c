#include "pd.h"

void
calculate_contact_force(int i, int j)
{
	int 			q;
	double          fs_old_mag, fs_old_x, fs_old_y, fs_old_z;
	double          disp_mag, ddisp_mag, disp_old_mag, disp_old_x, disp_old_y, disp_old_z;
	double          norm_old_x, norm_old_y, norm_old_z;
	double          norm_x, norm_y, norm_z, xij, yij, zij, vxij, vyij,
	                vzij, angvi_x, angvi_y, angvi_z, angvj_x, angvj_y,
	                angvj_z;
	double			sign_of_disp, sign_of_fs, alpha, alpha_max, alpha_y, separation, separation_sq, radii_sum, radii_sum_sq;
	double          E_eff, YS, R_eff, Kn, Kt, G_eff, M_eff;
	double          Rp, fn, fn_max, fn_test, fs_x, fs_y,
	                fs_z, fs_mag;
	double          ddisp_x, ddisp_y, ddisp_z, disp_x, disp_y, disp_z, tx, ty, tz,
	                fx, fy, fz, vn, rparti, rpartj;
	double 			odisp_x, odisp_y, odisp_z, ofs_x, ofs_y, ofs_z;
	double          contact_radius, dalpha, F_yield;
	
	/* find the proper index for this potential pair of particles, note that the index i must be greater than j */
	
	q = pair(i, j);
	
	/* since j will always be a "real" particle, we set its radius and angular velocity values */
	
	rpartj = particle_radius[j];
	
	angvj_x = particle_angular_velocity_x[j];
	angvj_y = particle_angular_velocity_y[j];
	angvj_z = particle_angular_velocity_z[j];

	/* if i is *not* a regular particle, we use the wall routines to get the wall positions and velocities */
	/* as well as to calculate the differences in position and velocity */
	
	if  (i>=number_of_particles) {
	
		rparti = 0.0;
		
		xij = wall_x(i, particle_x[j], particle_y[j], particle_z[j]) - particle_x[j];
		if (xij > 0.5 * Nx && XPERIODIC)
			xij -= Nx;
		else if (xij < -0.5 * Nx && XPERIODIC)
			xij += Nx;
			
		yij = wall_y(i, particle_x[j], particle_y[j], particle_z[j]) - particle_y[j];
		if (yij > 0.5 * Ny && YPERIODIC)
			yij -= Ny;
		else if (yij < -0.5 * Ny && YPERIODIC)
			yij += Ny;
			
		zij = wall_z(i, particle_x[j], particle_y[j], particle_z[j]) - particle_z[j];
		if (zij > 0.5 * Nz && ZPERIODIC)
			zij -= Nz;
		else if (zij < -0.5 * Nz && ZPERIODIC)
			zij += Nz;
		
		angvi_x = angvi_y = angvi_z = 0.0;

		vxij = wall_vx(i) - particle_velocity_x[j];
		vyij = wall_vy(i) - particle_velocity_y[j];
		vzij = wall_vz(i) - particle_velocity_z[j];
		
		/* we always assume that the wall is the same material as the particle, for now */
		
		YS = particle_YS[j];
		
		E_eff = particle_E_star[j];
		G_eff = particle_G_star[j];
		R_eff = rpartj;
		M_eff = particle_mass[j];
		
	}
	
	else { /* for all "regular" particles, we can simply take the differences in position and velocity directly */

		rparti = particle_radius[i];

		xij = particle_x[i] - particle_x[j];
		if (xij > 0.5 * Nx && XPERIODIC)
			xij -= Nx;
		else if (xij < -0.5 * Nx && XPERIODIC)
			xij += Nx;

		yij = particle_y[i] - particle_y[j];
		if (yij > 0.5 * Ny && YPERIODIC)
			yij -= Ny;
		else if (yij < -0.5 * Ny && YPERIODIC)
			yij += Ny;

		zij = particle_z[i] - particle_z[j];
		if (zij > 0.5 * Nz && ZPERIODIC)
			zij -= Nz;
		else if (zij < -0.5 * Nz && ZPERIODIC)
			zij += Nz;

		angvi_x = particle_angular_velocity_x[i];
		angvi_y = particle_angular_velocity_y[i];
		angvi_z = particle_angular_velocity_z[i];

		vxij = particle_velocity_x[i] - particle_velocity_x[j];
		vyij = particle_velocity_y[i] - particle_velocity_y[j];
		vzij = particle_velocity_z[i] - particle_velocity_z[j];
		
		YS = particle_YS[j];
		if (YS > particle_YS[i])
			YS = particle_YS[i];

		E_eff = particle_E_star[i] * particle_E_star[j] / (particle_E_star[i] + particle_E_star[j]);
		G_eff = particle_G_star[i] * particle_G_star[j] / (particle_G_star[i] + particle_G_star[j]);
		R_eff = rparti * rpartj / (rparti + rpartj);
		M_eff = particle_mass[i] * particle_mass[j] / (particle_mass[i] + particle_mass[j]);
		
	}

	radii_sum = rparti + rpartj;
	radii_sum_sq = radii_sum * radii_sum;

	separation_sq = xij * xij + yij * yij + zij * zij;
	
	/* check for particle overlap before taking the sqrt */
	
	if (radii_sum_sq > separation_sq) {
	
		/* we check if the current contact is a continuation (i.e., only 1 time step away from the pervious contact) */
		/* otherwise, we erase all of the previous contact data */
		
		if ((number_of_time_steps-particle_number_of_time_steps[q])>1) erase_contact_data(q);
	
		/* recover the old contact data that is necessary to continue the contact */
		
		Rp = particle_Rp[q];
		/* if Rp is zero, we know that this is a new contact, so the average curvature is the true value */
		if (Rp == 0.0) Rp = R_eff;
		
		norm_old_x = particle_norm_old_x[q];
		norm_old_y = particle_norm_old_y[q];
		norm_old_z = particle_norm_old_z[q];
		
		fn_max = particle_force_normal_max[q];
		alpha_max = particle_disp_normal_max[q];
		
		fs_old_mag = particle_force_tan_old_mag[q];
		fs_old_x = particle_force_tan_old_x[q];
		fs_old_y = particle_force_tan_old_y[q];
		fs_old_z = particle_force_tan_old_z[q];

		sign_of_disp = particle_sign_of_disp[q];
	
		disp_old_mag = particle_disp_tan_old_mag[q];
		disp_old_x = particle_disp_tan_old_x[q];
		disp_old_y = particle_disp_tan_old_y[q];
		disp_old_z = particle_disp_tan_old_z[q];
		
		/* calculate the degree of overlap, setting its value to be positive when there is overlap */
		
		separation = sqrt(xij * xij + yij * yij + zij * zij);
		alpha = radii_sum - separation;

		/* our normal vector points from particle j to particle i */
		
		norm_x = xij / separation;
		norm_y = yij / separation;
		norm_z = zij / separation;

		/* calculate the velocity parallel to the normal vector */
		
		vn = norm_x * vxij + norm_y * vyij + norm_z * vzij;

		/* calculate the relative normal displacement increment */
		
		dalpha = vn * dt;

		/* determine normal force */
		
		if (contact_mechanics == ELASTIC) {
		
			contact_radius = sqrt(R_eff * alpha);

			Kn = (4.0/3.0) * E_eff/R_eff;
			fn = Kn * contact_radius*contact_radius*contact_radius - damping * vn * alpha;

		}
		
		else if (contact_mechanics == PLASTIC) {
		
			/* determine the yield force using py = 2.5 YS (as suggest by Thornton (1997)) */
			
			F_yield = 15.625 * YS * YS * YS * M_PI * M_PI * M_PI * R_eff * R_eff / (6.0 * E_eff * E_eff);

			/* for loading */

			if (vn <= 0.0) {

				/* if we are below the yield force, we have Hertzian loading */
				/* so let's assume that we are still elastic */
				
				/* we use Rp in case this is reloading */
				contact_radius = sqrt(Rp * alpha);
				Kn = (4.0/3.0) * E_eff/Rp;
				fn = Kn * contact_radius*contact_radius*contact_radius;
				
				/* if our guess was wrong, we need to fix it */
				if (fn > F_yield) {
				
					/* we are now in the linear portion, so we need the alpha intercept */
					alpha_y = (M_PI * M_PI * 6.25 * YS * YS * Rp)/(4.0 * E_eff * E_eff);
					fn = F_yield + M_PI * Rp * 2.5 * YS * (alpha - alpha_y);
					
				}
				
			/* remember this for some contact models */

				if (fn > fn_max) {
					fn_max = fn;
					alpha_max = alpha;
				}

			}
			
			/* for unloading */

			else {
				double alpha_p;
				
				/* we need to calculate Rp if we have not already done so */
				/* Rp will be R_eff unless we have exceeded the yield point */
				if (fn_max > F_yield && Rp == R_eff)
					Rp = R_eff*(F_yield/fn_max)*pow((2.0*fn_max+F_yield)/(3.0*F_yield), 1.5);
					
				Kn = (4.0/3.0) * E_eff/Rp;
			
				contact_radius = sqrt(Rp * alpha);
				alpha_p = alpha_max - pow(((3.0/4.0)*(fn_max/(E_eff*sqrt(Rp)))), (2.0/3.0));
				fn = (4.0/3.0)*E_eff*sqrt(Rp*(alpha-alpha_p))*(alpha-alpha_p);

			}
			
			fn -= 2.0*BETA*sqrt(Kn*M_eff)*vn;
			
		}
		
		else {
			printf("Undefined force model!\n");
			exit(1);
		}
		
		if (fn < 0.0)
			fn = 0.0;

		/* start tangential force */

		fn_test = friction_coefficient * fn;

		/* calculate the relative tangential displacement increment */
		/* as n X v X n plus rotation */

		ddisp_x = (((norm_z * vxij - norm_x * vzij) * norm_z - (norm_x * vyij - norm_y * vxij) * norm_y) - (rparti * (angvi_y * norm_z - angvi_z * norm_y) + rpartj * (angvj_y * norm_z - angvj_z * norm_y)))*dt;
		ddisp_y = (((norm_x * vyij - norm_y * vxij) * norm_x - (norm_y * vzij - norm_z * vyij) * norm_z) - (rparti * (angvi_z * norm_x - angvi_x * norm_z) + rpartj * (angvj_z * norm_x - angvj_x * norm_z)))*dt;
		ddisp_z = (((norm_y * vzij - norm_z * vyij) * norm_y - (norm_z * vxij - norm_x * vzij) * norm_x) - (rparti * (angvi_x * norm_y - angvi_y * norm_x) + rpartj * (angvj_x * norm_y - angvj_y * norm_x)))*dt;

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
		 
		if (norm_old_x != 0.0 || norm_old_y != 0.0 || norm_old_z != 0.0) {
			double          n_dot_n, n_X_n_x, n_X_n_y, n_X_n_z;
			double          A1, A2, A3, theta;	/* trig functions in rotation matrix */

			/*
			 * In order to calculate the angle of rotation, we need the
			 * dot product between the old and new normal vectors
			 */

			n_dot_n = norm_x * norm_old_x + norm_y * norm_old_y + norm_z * norm_old_z;

			/*
			 * The axis of the rotation comes from the cross product
			 * between the old and new normal vectors
			 */

			cross(norm_x, norm_y, norm_z, norm_old_x, norm_old_y, norm_old_z, &n_X_n_x, &n_X_n_y, &n_X_n_z);

			/* We now calculate the trig functions needed for the rotation matrix */

			if (n_dot_n > 1.0) // if the normal vectors stayed the same then we rotate nothing
				theta = 0.0;
			else	// otherwise the angle comes from the arc cosine
				theta = acos(n_dot_n);

			A3 = cos(theta);
			A1 = 1.0 - A3;
			A2 = sin(theta);

			/* we adjust the old tangential displacement to account for the rotation */
			/* first we need to save the old displacement values in temporary odisp vector */
		
			odisp_x = disp_old_x;
			odisp_y = disp_old_y;
			odisp_z = disp_old_z;
		
			/* use the rotation matrix to rotate about the axis n_X_n by angle theta */
			
			disp_old_x = (n_X_n_x * n_X_n_x * A1 + A3) * odisp_x + (n_X_n_x * n_X_n_y * A1 - n_X_n_z * A2) * odisp_y + (n_X_n_x * n_X_n_z * A1 + n_X_n_y * A2) * odisp_z;
			disp_old_y = (n_X_n_x * n_X_n_y * A1 + n_X_n_z * A2) * odisp_x + (n_X_n_y * n_X_n_y * A1 + A3) * odisp_y + (n_X_n_y * n_X_n_z * A1 - n_X_n_x * A2) * odisp_z;
			disp_old_z = (n_X_n_x * n_X_n_z * A1 + n_X_n_y * A2) * odisp_x + (n_X_n_z * n_X_n_y * A1 + n_X_n_x * A2) * odisp_y + (n_X_n_z * n_X_n_z * A1 + A3) * odisp_z;
		
			/* similarly we adjust the old tangential force to account for the rotation */
			/* first we need to save the old force values in temporary ofs vector */
		
			ofs_x = fs_old_x;
			ofs_y = fs_old_y;
			ofs_z = fs_old_z;
		
			/* use the rotation matrix to rotate about the axis n_X_n by angle theta */
			
			fs_old_x = (n_X_n_x * n_X_n_x * A1 + A3) * ofs_x + (n_X_n_x * n_X_n_y * A1 - n_X_n_z * A2) * ofs_y + (n_X_n_x * n_X_n_z * A1 + n_X_n_y * A2) * ofs_z;
			fs_old_y = (n_X_n_x * n_X_n_y * A1 + n_X_n_z * A2) * ofs_x + (n_X_n_y * n_X_n_y * A1 + A3) * ofs_y + (n_X_n_y * n_X_n_z * A1 - n_X_n_x * A2) * ofs_z;
			fs_old_z = (n_X_n_x * n_X_n_z * A1 + n_X_n_y * A2) * ofs_x + (n_X_n_z * n_X_n_y * A1 + n_X_n_x * A2) * ofs_y + (n_X_n_z * n_X_n_z * A1 + A3) * ofs_z;
			
		}
		
		/* calculate the magnitude of the total tangential displacement thus far in each direction */
		
		disp_x = disp_old_x + ddisp_x;
		disp_y = disp_old_y + ddisp_y;
		disp_z = disp_old_z + ddisp_z;
		
		/* determine whether the new increment of displacement is in the same */
	 	/* or opposite direction of the total tangential displacement, and */
	 	/* correct the old magnitude */

		if ((disp_x * disp_old_x + disp_y * disp_old_y + disp_z * disp_old_z) < 0.0)
			sign_of_disp *= -1.0;
		
		disp_mag = sign_of_disp*sqrt(disp_x*disp_x + disp_y*disp_y + disp_z*disp_z);
		
		/* determine the current net tangential displacement increment */

		if (disp_old_mag<0.0)
			disp_old_mag = -sqrt(disp_old_x*disp_old_x + disp_old_y*disp_old_y + disp_old_z*disp_old_z);
		else
			disp_old_mag = sqrt(disp_old_x*disp_old_x + disp_old_y*disp_old_y + disp_old_z*disp_old_z);
			
		ddisp_mag = disp_mag - disp_old_mag;
		
		/* determine if the tangential force is increasing or decreasing */

		if (fs_old_mag < 0.0)
			sign_of_fs = -1.0;
		else
			sign_of_fs = 1.0;
		
		fs_old_mag = sign_of_fs*sqrt(fs_old_x*fs_old_x + fs_old_y*fs_old_y + fs_old_z*fs_old_z);
		
		/* determine the tangential force */

		Kt = 8.0 * G_eff * contact_radius;
		
		fs_mag = fs_old_mag + Kt * ddisp_mag;
		
		//if (fs_mag*ddisp_mag>0.0) fs_mag -= 2.0*BETA*sqrt(Kt*M_eff)*ddisp_mag/dt;
		//else fs_mag += 2.0*BETA*sqrt(Kt*M_eff)*ddisp_mag/dt;
		
		/* decompose the tangential force into its primary directions (in system coordinates) */

		if (disp_mag != 0.0) {
			fs_x = -fs_mag * disp_x / disp_mag;
			fs_y = -fs_mag * disp_y / disp_mag;
			fs_z = -fs_mag * disp_z / disp_mag;
		} 
		
		else
			fs_mag = fs_x = fs_y = fs_z = 0.0;
		
		/* limit the tangential force to the Coulomb limit, if necessary */

		if (fs_mag != 0.0) {
			if (fn_test < fabs(fs_mag)) {
				fs_x = fs_x * fn_test / fabs(fs_mag);
				fs_y = fs_y * fn_test / fabs(fs_mag);
				fs_z = fs_z * fn_test / fabs(fs_mag);
				fs_mag = fn_test * fs_mag / fabs(fs_mag);
			}
		} 
		
		else
			fs_mag = fs_x = fs_y = fs_z = 0.0;
		
		fx = norm_x * fn + fs_x;
		fy = norm_y * fn + fs_y;
		fz = norm_z * fn + fs_z;
 
 		tx = (norm_y * fs_z - norm_z * fs_y);
		ty = (norm_z * fs_x - norm_x * fs_z);
		tz = (norm_x * fs_y - norm_y * fs_x);
	
 		if (DEBUG) {
			if (vinit==-8675309.0) vinit=vn;
			printf("%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %d %d %e %lld\n", particle_x[j], particle_y[j], fn*MASS*LENGTH/(TIME*TIME), alpha*LENGTH, fs_x*MASS*LENGTH/(TIME*TIME), -disp_x*LENGTH, fs_y*MASS*LENGTH/(TIME*TIME), -disp_y*LENGTH, fs_z*MASS*LENGTH/(TIME*TIME), -disp_z*LENGTH, disp_mag, disp_old_mag, ddisp_mag, ty, angvi_z, tz, i, j, -vn/vinit, number_of_time_steps);
    	}
    	
		particle_torque_x[j] -= rpartj * tx;
		particle_torque_y[j] -= rpartj * ty;
		particle_torque_z[j] -= rpartj * tz;

		particle_force_x[j] -= fx;
		particle_force_y[j] -= fy;
		particle_force_z[j] -= fz;
			
		particle_Rp[q] = Rp;
		
		particle_norm_old_x[q] = norm_x;
		particle_norm_old_y[q] = norm_y;
		particle_norm_old_z[q] = norm_z;
		
		particle_force_normal_max[q] = fn_max;
		particle_disp_normal_max[q] = alpha_max;
		
		particle_force_tan_old_mag[q] = fs_mag;
		particle_force_tan_old_x[q] = fs_x;
		particle_force_tan_old_y[q] = fs_y;
		particle_force_tan_old_z[q] = fs_z;
		
		particle_sign_of_disp[q] = sign_of_disp;
	
		particle_disp_tan_old_mag[q] = disp_mag;
		particle_disp_tan_old_x[q] = disp_x;
		particle_disp_tan_old_y[q] = disp_y;
		particle_disp_tan_old_z[q] = disp_z;

		if (i < number_of_particles) {
			particle_torque_x[i] -= rparti * tx;
			particle_torque_y[i] -= rparti * ty;
			particle_torque_z[i] -= rparti * tz;

			particle_force_x[i] += fx;
			particle_force_y[i] += fy;
			particle_force_z[i] += fz;

		}
		
		else if (i == (number_of_particles+bottom)) bottom_wall_force+=fx;
		
		particle_number_of_time_steps[q] = number_of_time_steps;
		
	}
	
	else if ((number_of_time_steps-particle_initial_time_step[q])<30) {
		//printf("Problem. The contact for %d and %d only took %lld time steps at time = %lld time steps! %e %e\n", i, j, (number_of_time_steps-particle_initial_time_step[q]), number_of_time_steps, radii_sum_sq, separation_sq);
		particle_initial_time_step[q] = -200;
	}
	
}
