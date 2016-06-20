#include "pd.h"

void
calculate_contact_force(int i, int j)
{
	sim_normal_contact_status fn_status;
	int 			q;
	double          fs_old_x = 0.0, fs_old_y = 0.0, fs_old_z = 0.0;
	double          norm_x, norm_y, norm_z, xij, yij, zij, vxij, vyij,
	                vzij, angvi_x, angvi_y, angvi_z, angvj_x, angvj_y,
	                angvj_z;
	double          alpha, separation, separation_sq, radii_sum, radii_sum_sq;
	double          E_eff, YS, R_eff, Kn, G_eff;
	double          Rp, fn, fn_old, fn_max = 0.0, fn_test, fs_x, fs_y,
	                fs_z, fs_mag;
	double          disp_x, disp_y, disp_z, tx = 0.0, ty = 0.0, tz = 0.0,
	                fx = 0.0, fy = 0.0, fz = 0.0, vn, rparti, rpartj;
	double          contact_radius, dalpha, F_yield;
	
	q = pair(i, j);
	rpartj = particle_radius[j];
	
	angvj_x = particle_angular_velocity_x[j];
	angvj_y = particle_angular_velocity_y[j];
	angvj_z = particle_angular_velocity_z[j];

        /* for all "regular" particles */
    
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

	

	radii_sum = rparti + rpartj;
    
	separation_sq = xij * xij + yij * yij + zij * zij;
	radii_sum_sq = radii_sum * radii_sum;

    
	if (radii_sum_sq > separation_sq) {
    
		if ((number_of_time_steps-particle_number_of_time_steps[q])>1) erase_contact_data(q);
	
		separation = sqrt(xij * xij + yij * yij + zij * zij);
        
		alpha = radii_sum - separation;  //alpha is way too high

		norm_x = xij / separation;
		norm_y = yij / separation;
		norm_z = zij / separation;

		Rp = particle_Rp[q];
		fn_old = particle_force_normal_old[q]; //why is this zero?

		fn_max = particle_force_normal_max[q];
		fn_status = particle_force_normal_status[q];
		fs_old_x = particle_force_tan_old_x[q];
		fs_old_y = particle_force_tan_old_y[q];
		fs_old_z = particle_force_tan_old_z[q];

		vn = norm_x * vxij + norm_y * vyij + norm_z * vzij;
        
        /*
        if (i==672)
        {
            printf("vn %e, vxij %e, vyij %e, vzij %e alpha %e, fn_old %e, particle_velocity_x[i] %e\n", vn, vxij, vyij, vzij, alpha, fn_old, particle_velocity_x[i]);
            
            printf("fn_old %e ******* \n", fn_old);
        }
         */
		dalpha = -vn * dt;

	F_yield = YS * YS * YS * M_PI * M_PI * M_PI * R_eff * R_eff / (6.0 * E_eff * E_eff);

	/* for loading */

	if (vn <= 0.0) {  // particles are loading for no reason

		if (fn_status == NONE || fn_old >= fn_max)
			Rp = R_eff;
        
		fn_status = LOADING;
    
		if (fn_old < fn_max && fn_old != 0.0)
			contact_radius = pow((3.0 * Rp * fn_old / (4.0 * E_eff)), 0.333);
		else
			contact_radius = sqrt(Rp * alpha);
        
		Kn = 2.0 * E_eff * contact_radius;

		/* elastic loading force */

		fn = fn_old + Kn * dalpha;

		/* if we have visco-damping */

		if (contact_mechanics == ELASTIC)
			fn -= damping * vn * alpha;

		/* if we have plastic deformation */

		else if (contact_mechanics == PLASTIC)
			if (fn > F_yield) {
				contact_radius = sqrt((2.0 * fn_old + F_yield) / (2.0 * M_PI * YS));
				Kn = M_PI * R_eff * YS;
				fn = fn_old + Kn * dalpha;
			}
		/* remember this for some contact models */
        
		if (fn > fn_max)
			fn_max = fn;
        
        /*
        if (i==2997 || i==3017)
        {
            printf("START_loaidng\n");
            printf("fn %e, Loading Kn %e, dalpha %e, Kn * dalpha %e, i %i, j %i\n",fn, Kn, dalpha, Kn * dalpha, i, j);
            printf("fn_old %e, fn_max  %e\n", fn_old, fn_max);
        }
        */
	}
        
	/* for unloading */
	else {
		/* if we have plastic deformation */

		if (contact_mechanics == PLASTIC) {//It is plastic
			if (fn_status == LOADING) {
				fn_status = UNLOADING;
				contact_radius = sqrt((2.0 * fn_old + F_yield) / (2.0 * M_PI * YS));
                
				Rp = 4.0 * E_eff * contact_radius * contact_radius * contact_radius / (3.0 * fn_max);
                
                /*
                if (i==2997 || i==3017)
                {
                    printf("START_UN_loaidng\n");
                    printf("Rp %e, E_eff %e, contact_radius %e, fn_max %e\n",Rp, E_eff, contact_radius, fn_max);
                }
                 */
			}
		}
		contact_radius = pow((3.0 * Rp * fn_old / (4.0 * E_eff)), 0.333);
		Kn = 2.0 * E_eff * contact_radius;
        
		fn = fn_old + Kn * dalpha;
        
		if (contact_mechanics == ELASTIC)
			fn -= damping * vn * alpha;

        
	}

	if (fn < 0.0)
		fn = 0.0;

		/* start tangential force */

		fn_test = friction_coefficient * fn;

		disp_x = ((norm_z * vxij - norm_x * vzij) * norm_z - (norm_x * vyij - norm_y * vxij) * norm_y - rparti * (angvi_y * norm_z - angvi_z * norm_y) - rpartj * (angvj_y * norm_z - angvj_z * norm_y)) * dt;
		disp_y = ((norm_x * vyij - norm_y * vxij) * norm_x - (norm_y * vzij - norm_z * vyij) * norm_z - rparti * (angvi_z * norm_x - angvi_x * norm_z) - rpartj * (angvj_z * norm_x - angvj_x * norm_z)) * dt;
		disp_z = ((norm_y * vzij - norm_z * vyij) * norm_y - (norm_z * vxij - norm_x * vzij) * norm_x - rparti * (angvi_x * norm_y - angvi_y * norm_x) - rpartj * (angvj_x * norm_y - angvj_y * norm_x)) * dt;

		fs_x = fs_old_x - 8.0 * G_eff * contact_radius * disp_x;
		fs_y = fs_old_y - 8.0 * G_eff * contact_radius * disp_y;
		fs_z = fs_old_z - 8.0 * G_eff * contact_radius * disp_z;
        
		if ((fs_mag = sqrt(fs_x * fs_x + fs_y * fs_y + fs_z * fs_z)) > 0.0) {
			if (fn_test < fs_mag) {

				fs_x = fs_x * fn_test / fs_mag;
				fs_y = fs_y * fn_test / fs_mag;
				fs_z = fs_z * fn_test / fs_mag;
				fs_mag = fn_test;

			}
		} else
			fs_mag = fs_x = fs_y = fs_z = 0.0;
			
		tx = (norm_y * fs_z - norm_z * fs_y);
		ty = (norm_z * fs_x - norm_x * fs_z);
		tz = (norm_x * fs_y - norm_y * fs_x);
        
		fx = norm_x * fn + fs_x;
		fy = norm_y * fn + fs_y;
		fz = norm_z * fn + fs_z;
 
 		if (DEBUG) {
			if (vinit==-8675309.0) vinit=vn;
			printf("%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %d %d\n", particle_x[j], particle_y[j], fn*MASS*LENGTH/(TIME*TIME), alpha*LENGTH, fs_x*MASS*LENGTH/(TIME*TIME), disp_x*LENGTH, fs_y*MASS*LENGTH/(TIME*TIME), disp_y*LENGTH, fs_z*MASS*LENGTH/(TIME*TIME), disp_z*LENGTH, angvi_x, tx, angvi_y, ty, angvi_z, tz, i, j);
    	}
    	
    	incontact=1;
    
 		particle_torque_x[j] -= rpartj * tx;
		particle_torque_y[j] -= rpartj * ty;
		particle_torque_z[j] -= rpartj * tz;

		particle_force_x[j] -= fx;
		particle_force_y[j] -= fy;
		particle_force_z[j] -= fz;
			
		particle_touching[q] = 1;
		particle_Rp[q] = Rp;
		particle_force_normal_old[q] = fn;
		particle_force_normal_max[q] = fn_max;
		particle_force_normal_status[q] = fn_status;
		particle_force_tan_old_x[q] = fs_x;
		particle_force_tan_old_y[q] = fs_y;
		particle_force_tan_old_z[q] = fs_z;

        particle_torque_x[i] -= rparti * tx;
        particle_torque_y[i] -= rparti * ty;
        particle_torque_z[i] -= rparti * tz;

        particle_force_x[i] += fx;
        particle_force_y[i] += fy;
        particle_force_z[i] += fz;
        
        particle_number_of_time_steps[q] = number_of_time_steps;

        if (Gdirection==-1 && j < number_of_moving_wall_particles)
        {
            wall_force -= fy;
           // printf("fy %e, wall_force %e\n", fy,wall_force);
        }
        else if (Gdirection==1 && j < number_of_wall_particles && j >= number_of_moving_wall_particles)
        {
            wall_force -= fy;  //cheack this line later. Maybe mistake! No, looks right. Wall_force switchs signs
        }
        /*
        if (i==953 && j==572)
        {
            printf("fn %e i %i j %i\n", fn,i,j);
            printf("fx %e, particle_force_normal_old[q] %e\n",fx, particle_force_normal_old[q]);
            printf("\n");
        }
        */
        
	}
	


}
