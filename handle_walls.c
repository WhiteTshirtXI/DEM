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
handle_walls()
{
	int             i;
    //***********************************
		double          px, py, pz;
		for (i = 0; i < number_of_moving_wall_particles; i++) {
             //gdirection=-1 initially, force is negative. on the down direction, so have negative sign
            particle_angular_velocity_x[i] = particle_angular_velocity_y[i] =
            particle_angular_velocity_z[i]  = 0.0;
            particle_force_x[i] = particle_force_y[i] = particle_force_z[i] = 0.0;
            particle_torque_x[i] = particle_torque_y[i] = particle_torque_z[i] = 0.0;
            
            particle_velocity_z[i] = 0.0;
            if (Gdirection==-1)
            {
            wall_force = -wall_mass; 
			px = particle_x[i] += move_particle(wall_mass, particle_velocity_x[i], 0.0);
			if (px < 0.0)
				particle_x[i] = px += (double) Nx;
			else if (px >= Nx)
				particle_x[i]= px -= (double) Nx;

            py = particle_y[i] += move_particle(wall_mass, particle_velocity_y[i], wall_force);
            if (py < 0.0)
                particle_y[i]= py += (double) Ny;
            else if (py >= Ny)
                particle_y[i] = py -= (double) Ny;

			pz = particle_z[i] += move_particle(wall_mass, 0.0, 0.0);  //change position
			if (pz < 0.0)
				particle_z[i] = pz += (double) Nz;
			else if (pz >= Nz)
				particle_z[i] = pz -= (double) Nz;
			
			//particle[i].velocity_x = 0.0;
			particle_velocity_x[i] = wall_speed;
			//particle[i].velocity_y = 0.0;
			particle_velocity_y[i] += update_velocity(wall_mass, wall_force);
            }
    
            else if (Gdirection==1)
            {
                wall_force = wall_mass;
			px = particle_x[i] += move_particle(wall_mass, particle_velocity_x[i], 0.0);
			if (px < 0.0)
				particle_x[i] = px += (double) Nx;
			else if (px >= Nx)
				particle_x[i] = px -= (double) Nx;

            py = particle_y[i] += 0.0;
            if (py < 0.0)
                particle_y[i] = py += (double) Ny;
            else if (py >= Ny)
                particle_y[i] = py -= (double) Ny;
			
			pz = particle_z[i] +=0.0;
			//pz = particle[i].z += move_particle(wall_mass, 0.0, 0.0);
			if (pz < 0.0)
				particle_z[i] = pz += (double) Nz;
			else if (pz >= Nz)
				particle_z[i] = pz -= (double) Nz;

			particle_velocity_x[i] = 0.0;
			//particle[i].velocity_x = wall_speed;  //move towards the opposite direction
			particle_velocity_y[i]= 0.0;
            }
        }

		for (i = number_of_moving_wall_particles; i < number_of_wall_particles; i++) {
            particle_velocity_z[i] = 0.0;
            
            particle_angular_velocity_x[i] = particle_angular_velocity_y[i] =
            particle_angular_velocity_z[i]  = 0.0;
            particle_force_x[i] = particle_force_y[i] = particle_force_z[i] = 0.0;
            particle_torque_x[i] = particle_torque_y[i] = particle_torque_z[i] = 0.0;
            if (Gdirection==-1)
            {
                wall_force = -wall_mass;
                px = particle_x[i] += move_particle(wall_mass, particle_velocity_x[i], 0.0);
                if (px < 0.0)
                    particle_x[i] = px += (double) Nx;
                else if (px >= Nx)
                    particle_x[i] = px -= (double) Nx;

                py = particle_y[i] += 0.0;
                if (py < 0.0)
                    particle_y[i] = py += (double) Ny;
                else if (py >= Ny)
                    particle_y[i] = py -= (double) Ny;
                pz = particle_z[i] += 0.0;
                if (pz < 0.0)
                    particle_z[i] = pz += (double) Nz;
                else if (pz >= Nz)
                    particle_z[i] = pz -= (double) Nz;
                particle_velocity_x[i] = 0.0;
                //particle[i].velocity_x = -wall_speed/2.0;  //move towards the opposite direction
                particle_velocity_y[i] = 0.0;
        }
            else if (Gdirection==1)
            {
            wall_force = wall_mass; 
                px = particle_x[i]+= move_particle(wall_mass, particle_velocity_x[i], 0.0);
                if (px < 0.0)
                    particle_x[i] = px += (double) Nx;
                else if (px >= Nx)
                    particle_x[i] = px -= (double) Nx;
                
                py = particle_y[i] += move_particle(wall_mass, particle_velocity_y[i], -wall_force);
                if (py < 0.0)
                    particle_y[i] = py += (double) Ny;
                else if (py >= Ny)
                    particle_y[i] = py -= (double) Ny;
                
                pz = particle_z[i] += move_particle(wall_mass, 0.0, 0.0);
                if (pz < 0.0)
                    particle_z[i] = pz += (double) Nz;
                else if (pz >= Nz)
                    particle_z[i] = pz -= (double) Nz;
                
                //particle[i].velocity_x = 0.0;
                particle_velocity_x[i] = -wall_speed;
                //particle[i].velocity_y = 0.0;
                particle_velocity_y[i] += update_velocity(wall_mass, wall_force);

            }

		}
}

void 
handle_moving_walls()
{
	int             i;

	{
        //wall_force = -wall_mass;  //gdirection=-1 initially, force is negative. on the down direction, so have negative sign
        
        if (Gdirection==-1)
        {
            for (i = 0; i < number_of_moving_wall_particles; i++)
			{
                particle_velocity_x[i]= wall_speed;
                particle_velocity_y[i] += update_velocity(wall_mass, wall_force);
                particle_velocity_z[i] = 0.0;
                
                particle_angular_velocity_x[i] = particle_angular_velocity_y[i] =
                particle_angular_velocity_z[i]  = 0.0;
			}
			
		 for (i = number_of_moving_wall_particles; i < number_of_wall_particles; i++)
             {
                 particle_velocity_x[i] = 0.0;
                 particle_velocity_y[i] = 0.0;
                 particle_velocity_z[i] = 0.0;
                 
                 particle_angular_velocity_x[i] = particle_angular_velocity_y[i] =
                 particle_angular_velocity_z[i]  = 0.0;
             }
        }
        
        else if (Gdirection==1)
        {
            for (i = 0; i < number_of_moving_wall_particles; i++)
			{
                particle_velocity_x[i] = 0.0;
                particle_velocity_y[i] = 0.0;
                particle_velocity_z[i] = 0.0;
                
                particle_angular_velocity_x[i] = particle_angular_velocity_y[i] =
                particle_angular_velocity_z[i]  = 0.0;
			}
			
            for (i = number_of_moving_wall_particles; i < number_of_wall_particles; i++)
             {
                 particle_velocity_x[i] =-wall_speed;
                 particle_velocity_y[i] += update_velocity(wall_mass, wall_force);  // updated wall_force
                 particle_velocity_z[i] = 0.0;
                 
                 particle_angular_velocity_x[i] = particle_angular_velocity_y[i] =
                 particle_angular_velocity_z[i]  = 0.0;
            }
        }
	}
}
