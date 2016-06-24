
#include <stdexcept>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <cmath>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <math.h>

int number_of_particles=0000;
float Volumes=0;
double simulation_size_z = 8;
double simulation_size_x = 20;
double simulation_size_y = 40;

int particlenumber (int begin_file)
{
    int i_0, BIGNUMBER = 90000;
    float junk_0;
    float x=0,y=0,z=0,radius=0, v_x=0;
    int color=0;
    FILE *fp_0;
    
    char particlenumberfile[256];
    
    sprintf (particlenumberfile, "file%04d", begin_file);
    
    fp_0 = fopen (particlenumberfile, "r");
    
    for (i_0 = 0; i_0 < BIGNUMBER; i_0++) {
        if ((fscanf(fp_0, "%e %e %e %e %e %e %e %e %e %e %e %d\n", &x, &y, &z, &v_x, &junk_0, &junk_0, &junk_0, &junk_0, &junk_0, &radius, &junk_0, &color) == EOF))
        break;
    }
    fclose (fp_0);
    
    return (i_0);    
}

double ParticleVolumes (int begin_file)
{
    int i_0;
    float junk_0;
    float x=0,y=0,z=0,radius=0, v_x=0;
    int color=0;
    FILE *fp_0;
    int Number_Particle=0;
    double ParticleVolumes;
    
    char particlenumberfile[256];
    
    sprintf (particlenumberfile, "file%04d", begin_file);
    
    fp_0 = fopen (particlenumberfile, "r");
    
    for (i_0 = 0; i_0 < number_of_particles; i_0++) {
        fscanf(fp_0, "%e %e %e %e %e %e %e %e %e %e %e %d\n", &x, &y, &z, &v_x, &junk_0, &junk_0, &junk_0, &junk_0, &junk_0, &radius, &junk_0, &color);
        if (color != 0)
        {
            ParticleVolumes=ParticleVolumes+ 4.0/3.0*M_PI*pow(radius,3.0);
            Number_Particle++;
        }
            }
    fclose (fp_0);
    
    return (ParticleVolumes);
}

int main(int argc, char ** argv)
{
    int ft;
    char           infile[256];
    char           outfile[256];
    
    int begin_file=0000;
    int end_file=0000;
    double x=0,y=0,z=0,radius=0,junk_0=0,time=0.0, Y_position=0.0, reference_y_bottom=100.0, reference_y_top=0.0, wall_y, wall_up=0, wall_bottom=10, TotalVolume=0;
    int color=0;
    int blue=0;
    FILE *fp, *fpout;
    int i;
    double VolumeFraction;
    double Wall_distance;
    //int begin_file=0001;


    for (int i=1; i<argc;i++) {
        if (std::string(argv[i])=="-b") {
            i++;  //why needs i++ here?
            begin_file=atoi(argv[i]);
            std::cout<<"begin file number "<<begin_file<<std::endl;
        }
        
        if (std::string(argv[i])=="-e") {
            i++;  //why needs i++ here?
            end_file=atoi(argv[i]);
            std::cout<<"end file number "<<end_file<<std::endl;
        }
    }
    std::string usage= std::string("usage: -b to draw picture of a file number\n -e for the final file number");
    if (argc!=3) {
        std::cout << usage << std::endl;
    }
    //std::cout<<"flip_time_is ? "<<std::endl;
    //std::cin>>ft;

  //if using -b 1, begin_file will be file0001

    number_of_particles = particlenumber(begin_file);
    Volumes = ParticleVolumes(begin_file);
    
    //std::cout<<"number_of_particles "<<number_of_particles<<std::endl;
    			sprintf(outfile, "Wall_Distance");  //create output file
			    			if ((fpout = fopen(outfile, "w")) == NULL) {
							printf("Cannot Open output file\n");
							exit(1);
							}
    
    for (int g=begin_file; g<end_file; g++)
{
    wall_bottom=10;
    wall_up=0;
    
    std::cout<<g<<std::endl;
    sprintf (infile, "file%04d", g);
    if((fp=fopen(infile, "r"))==NULL)
    {
        std::cout<<"no file number "<<g<<std::endl;
        exit(1);
    }
    time=g/30.0;

    for (i=0;i<number_of_particles;i++)
    {
       fscanf(fp, "%le %le %le %le %le %le %le %le %le %le %le %d\n", &x, &y, &z, &junk_0, &junk_0, &junk_0, &junk_0, &junk_0, &junk_0, &radius, &junk_0, &color);

         {
             if (color ==0)
             
                {
                if (y<(simulation_size_y-10) && y>5)
                    wall_up=std::max(wall_up, y);  //number is kind of added in randomly, need to recheck if changed geometry
                else if (y<(simulation_size_y-20))
                    wall_bottom=std::min(wall_bottom,y);
                else if (y>(simulation_size_y-5))
                    wall_bottom=y;  //this is the periodic case, when bottom wall goes to the top
                }
        }
    }
    if (wall_bottom <(simulation_size_y-5.0))
    {
        //TotalVolume=(wall_up-wall_bottom)*simulation_size_z*simulation_size_x-(simulation_size_z*4.0+simulation_size_z*simulation_size_x)*(4.0/3.0*M_PI*(pow(0.5,3.0)));
        Wall_distance=wall_up-wall_bottom;
    }
    else if (wall_bottom > (simulation_size_y-5.0))
    {
        //TotalVolume=(wall_up-(wall_bottom-simulation_size_y))*simulation_size_z*simulation_size_x-(simulation_size_z*4.0+simulation_size_z*simulation_size_x)*(4.0/3.0*M_PI*(pow(0.5,3.0)));  //+1 for wall particle volumes *4 for fans
        Wall_distance=wall_up-(wall_bottom-simulation_size_y);
        std::cout<<"wall bottom - simulation size y "<<wall_bottom-simulation_size_y<<std::endl;
    }
    
    
    fprintf(fpout, "%e %e\n", time, Wall_distance);
    fclose(fp);
}
    fclose(fpout);
    std::cout<<" wall_bottom "<<wall_bottom<< " wall_top "<<wall_up<<" Wall distance  "<<Wall_distance <<std::endl;
    std::cout<< "Time & Wall_Distance "<<std::endl;
    
  return 0;
}
