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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>

#include "wall_position.c"

#define PARTICLE_RADIUS		(0.5*size_of_cells)
#define DISPERSION			0.025

int
main()
{

	int		i         , j, l, c, NX, NY, NZ , num = 0, num_small = 0, num_large = 0;
	double		xcounter, zcounter, r, X, Y, Z, x, y, z, dx, dy, dz, vx,
			vy           , vz;
	double		nx     , ny, nz, Rx, Ry, Rz, Dx, Dy, Dz, theta;
	double		avgvel_y = 0.0, avgvel_x = 0.0, avgvel_z = 0.0;
	double		XSTART, XEND, YSTART, YEND, ZSTART, ZEND;
	
	FILE           *fp;
	fp = fopen("position.in", "w");

	if (ZPERIODIC != 1 ) {
		ZSTART=(wall_z(number_of_particles+back, 0.0, 0.0, 0.0)+1.0);  //1.5+1=2.5
		ZEND=(wall_z(number_of_particles+front, 0.0, 0.0, 0.0)-1.0);  //Nz-1.5-1=3.0
	}
	else {
		ZSTART=0.5;
		ZEND=Nz-0.5;
	}
	
	if (YPERIODIC != 1 ) {  //not periodic
		YSTART=(wall_y(number_of_particles+bottom, 0.0, 0.0, 0.0)+1.0);  //1.5+1.0=2.5
		YEND=(wall_y(number_of_particles+top, 0.0, 0.0, 0.0)-1.0);      //Ny-1.5-1.0=7.5
	}
	else {
		YSTART=0.5;
		YEND=Nz-0.5;
	}
	
	if (XPERIODIC != 1 ) {
		XSTART=(wall_x(number_of_particles+left, 0.0, 0.0, 0.0)+1.0);
		XEND=(wall_x(number_of_particles+right, 0.0, 0.0, 0.0)-1.0);
	}
	else {
		XSTART=0.5;
		XEND=Nz-0.5;
	}
	
	printf("%e %e %e %e %e %e \n", XSTART, YSTART, ZSTART, XEND, YEND, ZEND);
	
	Dz = 2.0 * PARTICLE_RADIUS;
	nz = (ZEND - ZSTART) / Dz;
	NZ = (int)nz;
	Rz = 0.5 * (nz - (double)NZ ) / (double)NZ ;
	Dz += 2.0 * Rz;

	Dx = sin(60.0 * M_PI / 180.0) * 2.0 * PARTICLE_RADIUS;
	nx = (XEND - XSTART) / Dx;
	NX = (int)nx;
	Rx = 0.5 * (nx - (double)NX) / (double)NX;
	Dx += 2.0 * Rx;

	Dy = sin(60.0 * M_PI / 180.0) * 2.0 * PARTICLE_RADIUS;
	ny = (YEND - YSTART) / Dy;
	NY = (int)ny;
	Ry = 0.5 * (ny - (double)NY) / (double)NY;
	Dy += 2.0 * Ry;

	Y = YSTART + PARTICLE_RADIUS;
	xcounter = zcounter = 0.0;

	num = 0;
				
	for (l = 0; l < NY; l++) {
		xcounter += 1.0;

		if (fmod(xcounter, 2.0) == 0.0)
			X = XSTART + 0.5 * Dx;
		else
			X = XSTART;

		for (i = 0; i < NX; i++) {
			zcounter += 1.0;

			if (fmod(zcounter, 2.0) == 0.0)
				Z = ZSTART + 0.5 * Dz;
			else
				Z = ZSTART;

			for (j = 0; j < NZ; j++) {
				
				dx = (Rx) * (1.0 - 2.0*(double)arc4random_uniform(1000.0) / 1000.0);
				dy = (Ry) * (1.0 - 2.0*(double)arc4random_uniform(1000.0) / 1000.0);
				vy = 0.01 * 0.5 * (1.0 - 2.0*(double)arc4random_uniform(1000.0) / 1000.0);
				vx = 0.01 * 0.5 * (1.0 - 2.0*(double)arc4random_uniform(1000.0) / 1000.0);

				x = X + dx;
				y = Y + dy;

				vz = 0.0;
				z = Z;

				r = PARTICLE_RADIUS - DISPERSION + DISPERSION * (1.0 - 2.0*(double)arc4random_uniform(1000.0) / 1000.0);
				
				fprintf(fp, "%lf %lf %lf %lf %lf %lf 0.0 0.0 0.0 %lf 0.0 %d\n", x/size_of_cells, y/size_of_cells, z/size_of_cells, vx, vy, vz, r/size_of_cells, c);
				num++;
				
				Z += Dz;

			}

			X += Dx;

		}

		Y += Dy;

	}

	fclose(fp);
	
	printf("The number of total particles is %d\n", num);
	

}
