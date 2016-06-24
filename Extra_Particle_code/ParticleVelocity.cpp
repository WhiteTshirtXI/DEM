
#include <stdexcept>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <cmath>
#include <sstream>
#include <fstream>


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

int main(int argc, char ** argv)
{

    char           infile[256];
    char           outfile[256];
    
    int begin_file=0000;
    int end_file=0000;
    int number_of_particles=0000;
    double x=0,y=0,z=0,radius=0,junk_0=0, v_x=0;
    int color=0;
    
    FILE *fp, *fpout;
    int i;
    //int begin_file=0001;


    for (int i=1; i<argc;i++) {
        if (std::string(argv[i])=="-f") {
            i++;  //why needs i++ here?
            begin_file=atoi(argv[i]);
            std::cout<<"file number "<<begin_file<<std::endl;
        }
    }
    std::string usage= std::string("usage: -f which file to use to record velocity\n");
    if (argc!=1) {
        std::cout << usage << std::endl;
    }

  //if using -b 1, begin_file will be file0001


    number_of_particles = particlenumber(begin_file);
    
    //std::cout<<"number_of_particles "<<number_of_particles<<std::endl;
    			sprintf(outfile, "file%04d x_Velocity", begin_file);
			    			if ((fpout = fopen(outfile, "w")) == NULL) {
							printf("Cannot Open File 290\n");
							exit(1);
							}
				
							
	

    int g=begin_file;
{
    std::cout<<g<<std::endl;
    sprintf (infile, "file%04d", g);
    if((fp=fopen(infile, "r"))==NULL)
    {
        std::cout<<"no file number "<<g<<std::endl;
        exit(1);
    }
    fprintf(fpout, "y v_x\n");
    for (i=0;i<number_of_particles;i++)
    {
       fscanf(fp, "%le %le %le %le %le %le %le %le %le %le %le %d\n", &x, &y, &z, &v_x, &junk_0, &junk_0, &junk_0, &junk_0, &junk_0, &radius, &junk_0, &color);

         {
				fprintf(fpout, "%e %e\n", y, v_x);

         }
        

  }

    fclose(fp);
    
}
    fclose(fpout);

    
  return 0;
}
