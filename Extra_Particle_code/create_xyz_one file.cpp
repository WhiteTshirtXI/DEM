
#include <stdexcept>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <cmath>
#include <sstream>
#include <fstream>
#include <algorithm>


int particlenumber ()
{
    int i_0, BIGNUMBER = 90000;
    float junk_0;
    float x=0,y=0,z=0,radius=0, v_x=0;
    int color=0;
    FILE *fp_0;
    
    char particlenumberfile[256];
    
    sprintf (particlenumberfile, "mixing_compressed_half.in");//input file
    
    fp_0 = fopen (particlenumberfile, "r");
    
    for (i_0 = 0; i_0 < BIGNUMBER; i_0++) {
        if ((fscanf(fp_0, "%e %e %e %e %e %e %e %e %e %e %e %d\n", &x, &y, &z, &v_x, &junk_0, &junk_0, &junk_0, &junk_0, &junk_0, &radius, &junk_0, &color) == EOF))
        break;
    }
    fclose (fp_0);
    
    return (i_0);    
}

int main(int argc, char ** argv)
{
    int ft;
    char           infile[256];
    char           outfile[256];
    
    int begin_file=0000;
    int number_of_particles=0000;
    double x=0,y=0,z=0,radius=0,junk_0=0,time=0.0, Y_position=0.0, reference_y=100.0, wall_y;
    int color=0;
    int blue=0;
    FILE *fp, *fpout;
    int i;

    number_of_particles = particlenumber();

    			sprintf(outfile, "position.xyz");  //create output file
			    			if ((fpout = fopen(outfile, "w")) == NULL) {
							printf("Cannot Open a new File\n");
							exit(1);
							}

    fprintf(fpout, "%d\n", number_of_particles);
    fprintf(fpout, "using DMV to read position\n");
    sprintf (infile, "mixing_compressed_half.in");
    if((fp=fopen(infile, "r"))==NULL)
    {
        std::cout<<"no input file "<<std::endl;
        exit(1);
    }

    for (i=0;i<number_of_particles;i++)
    {
       fscanf(fp, "%le %le %le %le %le %le %le %le %le %le %le %d\n", &x, &y, &z, &junk_0, &junk_0, &junk_0, &junk_0, &junk_0, &junk_0, &radius, &junk_0, &color);

         {
             junk_0=0;
             if (color ==0)
             {
                fprintf(fpout, "H %le %le %le\n", x, y, z);
             }
             if (color ==3)
             {
                 fprintf(fpout, "O %le %le %le\n", x, y, z);
             }
             if (color ==7)
             {
                 fprintf(fpout, "C %le %le %le\n", x, y, z);
             }

        }
        
  }
    

    fclose(fp);
  return 0;
}
