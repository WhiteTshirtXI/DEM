#include "pd.h"

void
contact_detect()
{
	int             i, j, k, l, m, t;
	int             ix, iy, iz, ipx, ipy, ipz;
	List            zlist[Nz], ylist[2][Ny],
	                xlist[2][3][Nx];

	t = 0;

	/* zero all headers and done flags for iz-list */

	for (i = 0; i < Nz; i++) {
		zlist[i].done = 0;
		zlist[i].head = -1;
	}

	/*
	 * zero all headers and done flags for iy lists and iz and iz-1
	 * columns
	 */

	for (i = 0; i < Ny; i++) {

		ylist[0][i].done = 0;
		ylist[0][i].head = -1;

		ylist[1][i].done = 0;
		ylist[1][i].head = -1;

	}

	/*
	 * zero headers and done flags for ix lists in the iy, iy+1, and iy-1
	 * columns
	 */
	for (i = 0; i < Nx; i++) {

		xlist[0][0][i].done = 0;
		xlist[0][0][i].head = -1;

		xlist[0][1][i].done = 0;
		xlist[0][1][i].head = -1;

		xlist[1][0][i].done = 0;
		xlist[1][0][i].head = -1;

		xlist[1][1][i].done = 0;
		xlist[1][1][i].head = -1;

		xlist[1][2][i].done = 0;
		xlist[1][2][i].head = -1;

	}
	/* sorts particles into iz cells before we begin */
	for (i = 0; i < number_of_particles; i++) {
		iz = (int) particle_z[i];
		particle_nextz[i] = zlist[iz].head;
		zlist[iz].head = i;
	}

	/* loop over all particles */
	for (i = number_of_particles - 1; i >= 0; i--) {

		if (i >= number_of_wall_particles || (particle_x[i] < max_x && particle_x[i] > min_x && particle_y[i] < max_y && particle_y[i] > min_y && particle_z[i] < max_z && particle_z[i] > min_z)) {

			/* check which iz-list particle i is in */
			iz = (int) (particle_z[i]);
			/*
			 * if iy lists have not yet been generated for this
			 * iz list, do it
			 */
			 
			if (zlist[iz].done != 1) {
				/* mark it as done, so we do not do it again */
				zlist[iz].done = 1;
				/*
				 * we need to start at the top (head) of
				 * iz-list (already generated)
				 */
				j = zlist[iz].head;
				
				while (j != -1) {
					/*
					 * loop over all particles in this iz
					 * list to place them in iy lists
					 */
					iy = (int) (particle_y[j]);
					particle_nexty[j] = ylist[0][iy].head;
					ylist[0][iy].head = j;
					j = particle_nextz[j];
				}
				
				/*
				 * if iz == 0 we must fix to loop through
				 * list iz-1
				 */
				 
				if (iz == 0)
					ipz = Nz - 1;
				else
					ipz = iz - 1;
					
				j = zlist[ipz].head;
				
				while (j != -1) {
					/*
					 * loop over all particles in iz-1
					 * list to place them in iy lists
					 */
					iy = (int) (particle_y[j]);
					particle_nexty[j] = ylist[1][iy].head;
					ylist[1][iy].head = j;
					j = particle_nextz[j];
				}
				
				/*
				 * now that iy lists are generated for iz and
				 * iz-1, start again to get the ix lists (for
				 * iz, iz-1, and all the related iys) --
				 * resetting iz now....
				 */
				 
				j = zlist[iz].head;
				
				while (j != -1) {
				
					/*
					 * check which iy-list particle j is
					 * in
					 */
					iy = (int) (particle_y[j]);
					
					/*
					 * if ix lists have not yet been
					 * generated for this iy list, do it
					 */
					if (ylist[0][iy].done != 1) {
					
						/*
						 * mark it as done, so we do
						 * not do it again
						 */
						 
						ylist[0][iy].done = 1;
						/*
						 * we need to start at the
						 * top (head) of iy-list
						 * (already generated)
						 */
						 
						k = ylist[0][iy].head;
						while (k != -1) {
							/*
							 * loop over all
							 * particles in this
							 * iy, iz list, place
							 * them in ix lists
							 */
							ix = (int) (particle_x[k]);
							particle_nextx[k] = xlist[0][0][ix].head;
							xlist[0][0][ix].head = k;
							k = particle_nexty[k];
						}
						
						k = ylist[1][iy].head;
						while (k != -1) {
							/*
							 * loop over
							 * particles in iy,
							 * (iz-1) list, place
							 * them in ix lists
							 */
							ix = (int) (particle_x[k]);
							particle_nextx[k] = xlist[1][0][ix].head;
							xlist[1][0][ix].head = k;
							k = particle_nexty[k];
						}
						
						if (iy == 0)
							ipy = Ny - 1;
						else
							ipy = iy - 1;
							
						k = ylist[0][ipy].head;
						
						while (k != -1) {
							/*
							 * loop over
							 * particles in iy-1,
							 * (iz) list, place
							 * them in ix lists
							 */
							 
							ix = (int) (particle_x[k]);
							particle_nextx[k] = xlist[0][1][ix].head;
							xlist[0][1][ix].head = k;
							k = particle_nexty[k];
						}
						
						k = ylist[1][ipy].head;
						while (k != -1) {
							/*
							 * loop over
							 * particles in iy-1,
							 * (iz-1) list, place
							 * them in ix lists
							 */
							ix = (int) (particle_x[k]);
							particle_nextx[k] = xlist[1][1][ix].head;
							xlist[1][1][ix].head = k;
							k = particle_nexty[k];
						}
						
						if (iy == (Ny - 1))
							ipy = 0;
						else
							ipy = iy + 1;
							
						k = ylist[1][ipy].head;
						while (k != -1) {
							/*
							 * loop over
							 * particles in iy+1,
							 * (iz-1) list, place
							 * them in ix lists
							 */
							ix = (int) (particle_x[k]);
							particle_nextx[k] = xlist[1][2][ix].head;
							xlist[1][2][ix].head = k;
							k = particle_nexty[k];
						}
						
						/*
						 * lists are now complete,
						 * let's check if anyone is
						 * touching!
						 */
						 
						/*
						 * starting with iz ([0])
						 * list and loop through iy
						 */
						k = ylist[0][iy].head;
						while (k != -1) {
							/*
							 * check which
							 * ix-list particle i
							 * is in
							 */
							ix = (int) (particle_x[k]);
							/*
							 * is this the first
							 * time checking
							 * collisions for
							 * this iz-iy-ix
							 */
							if (xlist[0][0][ix].done != 1) {
								/*
								 * mark it as
								 * done, so
								 * we do not
								 * do it
								 * again
								 */
								xlist[0][0][ix].done = 1;
								/*
								 * loop
								 * through
								 * all ixs in
								 * ix-iy-iz
								 * list
								 */
								l = xlist[0][0][ix].head;
								while (l != -1) {
									m = particle_nextx[l];
									/*
									 * che
									 * ck
									 * col
									 * lis
									 * ion
									 * s
									 * bet
									 * wee
									 * n
									 * ix-
									 * iy-
									 * iz
									 * and
									 *
									 * ix-
									 * iy-
									 * iz
									 */
									while (m != -1) {
										if (l < m && m >= number_of_wall_particles)
											calculate_contact_force(m, l);
										else if (l >= number_of_wall_particles)
											calculate_contact_force(l, m);
										m = particle_nextx[m];
									}
									m = xlist[0][1][ix].head;
									/*
									 * che
									 * ck
									 * col
									 * lis
									 * ion
									 * s
									 * bet
									 * wee
									 * n
									 * ix-
									 * iy-
									 * iz
									 * and
									 *
									 * ix-
									 * iy-
									 * 1-i
									 * z
									 */
									while (m != -1) {
										if (l < m && m >= number_of_wall_particles)
											calculate_contact_force(m, l);
										else if (l >= number_of_wall_particles)
											calculate_contact_force(l, m);
										m = particle_nextx[m];
									}
									m = xlist[1][0][ix].head;
									/*
									 * che
									 * ck
									 * col
									 * lis
									 * ion
									 * s
									 * bet
									 * wee
									 * n
									 * ix-
									 * iy-
									 * iz
									 * and
									 *
									 * ix-
									 * iy-
									 * iz-
									 * 1
									 */
									while (m != -1) {
										if (l < m && m >= number_of_wall_particles)
											calculate_contact_force(m, l);
										else if (l >= number_of_wall_particles)
											calculate_contact_force(l, m);
										m = particle_nextx[m];
									}
									m = xlist[1][1][ix].head;
									/*
									 * che
									 * ck
									 * col
									 * lis
									 * ion
									 * s
									 * bet
									 * wee
									 * n
									 * ix-
									 * iy-
									 * iz
									 * and
									 *
									 * ix-
									 * iy-
									 * 1-i
									 * z-1
									 *
									 */
									while (m != -1) {
										if (l < m && m >= number_of_wall_particles)
											calculate_contact_force(m, l);
										else if (l >= number_of_wall_particles)
											calculate_contact_force(l, m);
										m = particle_nextx[m];
									}
									m = xlist[1][2][ix].head;
									/*
									 * che
									 * ck
									 * col
									 * lis
									 * ion
									 * s
									 * bet
									 * wee
									 * n
									 * ix-
									 * iy-
									 * iz
									 * and
									 *
									 * ix-
									 * iy+
									 * 1-i
									 * z-1
									 *
									 */
									while (m != -1) {
										if (l < m && m >= number_of_wall_particles)
											calculate_contact_force(m, l);
										else if (l >= number_of_wall_particles)
											calculate_contact_force(l, m);
										m = particle_nextx[m];
									}
									
									if (ix == 0)
										ipx = Nx - 1;
									else
										ipx = ix - 1;
										
									m = xlist[0][0][ipx].head;
									/*
									 * che
									 * ck
									 * col
									 * lis
									 * ion
									 * s
									 * bet
									 * wee
									 * n
									 * ix-
									 * iy-
									 * iz
									 * and
									 *
									 * ix-
									 * 1-i
									 * y-i
									 * z
									 */
									while (m != -1) {
										if (l < m && m >= number_of_wall_particles)
											calculate_contact_force(m, l);
										else if (l >= number_of_wall_particles)
											calculate_contact_force(l, m);
										m = particle_nextx[m];
									}
									
									m = xlist[0][1][ipx].head;
									while (m != -1) {
										if (l < m && m >= number_of_wall_particles)
											calculate_contact_force(m, l);
										else if (l >= number_of_wall_particles)
											calculate_contact_force(l, m);
										m = particle_nextx[m];
									}
									
									m = xlist[1][0][ipx].head;
									while (m != -1) {
										if (l < m && m >= number_of_wall_particles)
											calculate_contact_force(m, l);
										else if (l >= number_of_wall_particles)
											calculate_contact_force(l, m);
										m = particle_nextx[m];
									}
									
									m = xlist[1][1][ipx].head;
									while (m != -1) {
										if (l < m && m >= number_of_wall_particles)
											calculate_contact_force(m, l);
										else if (l >= number_of_wall_particles)
											calculate_contact_force(l, m);
										m = particle_nextx[m];
									}
									
									m = xlist[1][2][ipx].head;
									while (m != -1) {
										if (l < m && m >= number_of_wall_particles)
											calculate_contact_force(m, l);
										else if (l >= number_of_wall_particles)
											calculate_contact_force(l, m);
										m = particle_nextx[m];
									}
									
									if (ix == (Nx - 1))
										ipx = 0;
									else
										ipx = ix + 1;
									m = xlist[0][1][ipx].head;
									while (m != -1) {
										if (l < m && m >= number_of_wall_particles)
											calculate_contact_force(m, l);
										else if (l >= number_of_wall_particles)
											calculate_contact_force(l, m);
										m = particle_nextx[m];
									}
									
									m = xlist[1][0][ipx].head;
									while (m != -1) {
										if (l < m && m >= number_of_wall_particles)
											calculate_contact_force(m, l);
										else if (l >= number_of_wall_particles)
											calculate_contact_force(l, m);
										m = particle_nextx[m];
									}
									
									m = xlist[1][1][ipx].head;
									while (m != -1) {
										if (l < m && m >= number_of_wall_particles)
											calculate_contact_force(m, l);
										else if (l >= number_of_wall_particles)
											calculate_contact_force(l, m);
										m = particle_nextx[m];
									}
									
									m = xlist[1][2][ipx].head;
									while (m != -1) {
										if (l < m && m >= number_of_wall_particles)
											calculate_contact_force(m, l);
										else if (l >= number_of_wall_particles)
											calculate_contact_force(l, m);
										m = particle_nextx[m];
									}
									
									l = particle_nextx[l];
								}
							}
							k = particle_nexty[k];
						}
						k = ylist[0][iy].head;
						/*
						 * loop over all ix-iy-iz
						 * lists to discard them
						 */
						while (k != -1) {
							ix = (int) (particle_x[k]);
							xlist[0][0][ix].head = -1;
							xlist[0][0][ix].done = 0;
							k = particle_nexty[k];
						}
						k = ylist[1][iy].head;
						/*
						 * loop over ix-iy-iz1 lists
						 * to discard them
						 */
						while (k != -1) {
							ix = (int) (particle_x[k]);
							xlist[1][0][ix].head = -1;
							xlist[1][0][ix].done = 0;
							k = particle_nexty[k];
						}
						if (iy == 0)
							ipy = Ny - 1;
						else
							ipy = iy - 1;
						k = ylist[0][ipy].head;
						/*
						 * loop over all ix-iy1-iz
						 * lists to discard them
						 */
						while (k != -1) {
							ix = (int) (particle_x[k]);
							xlist[0][1][ix].head = -1;
							xlist[0][1][ix].done = 0;
							k = particle_nexty[k];
						}
						
						k = ylist[1][ipy].head;
						/*
						 * loop over all ix-iy1-iz1
						 * lists to discard them
						 */
						while (k != -1) {
							ix = (int) (particle_x[k]);
							xlist[1][1][ix].head = -1;
							xlist[1][1][ix].done = 0;
							k = particle_nexty[k];
						}
						
						if (iy == (Ny - 1))
							ipy = 0;
						else
							ipy = iy + 1;
						k = ylist[1][ipy].head;
						/*
						 * loop over all ix-iy2-iz1
						 * lists to discard them
						 */
						while (k != -1) {
							ix = (int) (particle_x[k]);
							xlist[1][2][ix].head = -1;
							xlist[1][2][ix].done = 0;
							k = particle_nexty[k];
						}
					}
					j = particle_nextz[j];
				}
				j = zlist[iz].head;
				/*
				 * loop over iz-list  particles to discard
				 * iz-iy lists
				 */
				while (j != -1) {
					iy = (int) (particle_y[j]);
					ylist[0][iy].head = -1;
					ylist[0][iy].done = 0;
					j = particle_nextz[j];
				}
				if (iz == 0)
					ipz = Nz - 1;
				else
					ipz = iz - 1;
				j = zlist[ipz].head;
				while (j != -1) {
					iy = (int) (particle_y[j]);
					ylist[1][iy].head = -1;
					ylist[1][iy].done = 0;
					j = particle_nextz[j];
				}
			}
		}
	}
}
