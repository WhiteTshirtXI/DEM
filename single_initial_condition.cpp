
#include <stdexcept>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <cmath>
#include <sstream>
#include <fstream>
#include <algorithm>

char           infile[256];
double Y_END=0.0, Y_Start=10.0;
int particlenumber ()
{
    int i_0, BIGNUMBER = 90000;
    float junk_0;
    float x=0,y=0,z=0,radius=0, v_x=0;
    int color=0;
    FILE *fp_0;
    
    char particlenumberfile[256];
    
    //sprintf (infile); 
    
    fp_0 = fopen (infile, "r");
    
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
    char           outfile[256];
    
    int begin_file=0000;
    int number_of_particles=0000;
    double x=0,y=0,z=0,radius=0,junk_0=0,time=0.0, Y_position=0.0, reference_y=100.0, wall_y;
    int color=0;
    int blue=0;
    FILE *fp, *fpout;
    int i;
    //int begin_file=0001;

    for (int i=1; i<argc;i++) {
        if (std::string(argv[i])=="-i") {
            i++;  //why needs i++ here?
            strcpy(infile, argv[i]); //cannot use atoi, which is used for convert string to int
            std::cout<<"begin file number "<<begin_file<<std::endl;
        }
        
    }
    std::string usage= std::string("usage: -i the file to be mixed");
    if (argc!=3) {
        std::cout << usage << std::endl;
    }

    number_of_particles = particlenumber();
    
    std::cout<<"number_of_particles "<<number_of_particles<<std::endl;
    			sprintf(outfile, "single_initial.in");                      //create output file
			    			if ((fpout = fopen(outfile, "w")) == NULL) {
							printf("Cannot Open a new File\n");
							exit(1);
							}

    //sprintf (infile, "thin_unmixing.in");  //initial file defiend by user
    if((fp=fopen(infile, "r"))==NULL)
    {
        std::cout<<"no input file "<<std::endl;
        exit(1);
    }

    for (i=0;i<number_of_particles;i++)
    {
        
        fscanf(fp, "%le %le %le %le %le %le %le %le %le %le %le %d\n", &x, &y, &z, &junk_0, &junk_0, &junk_0, &junk_0, &junk_0, &junk_0, &radius, &junk_0, &color);
        Y_END=std::max(y,Y_END);
        Y_Start=std::min(y,Y_Start);
        
    }
    fclose (fp); //need to close file since finished scan
    
    fp = fopen (infile, "r");
    std::cout<<"Y_start  Y_END  " <<Y_Start <<" "<< Y_END<<std::endl;
    std::cout<<"what the fuck going wrong"<<(Y_END+Y_Start)/2.0-0.5<<" "<< (Y_END+Y_Start)/2.0+0.5 <<std::endl;
    for (i=0;i<number_of_particles;i++)
    {
       fscanf(fp, "%le %le %le %le %le %le %le %le %le %le %le %d\n", &x, &y, &z, &junk_0, &junk_0, &junk_0, &junk_0, &junk_0, &junk_0, &radius, &junk_0, &color);

         {
             junk_0=0;
             if (color ==0)
             {
                    fprintf(fpout, "%le %le %le %le %le %le %le %le %le %le %le %d\n", x, y, z, junk_0, junk_0, junk_0, junk_0, junk_0, junk_0, radius, junk_0, color);
             }
				if (color==3 ||color==7 )
				{
                    if (y>=(Y_END+Y_Start)/2.0-0.5  && y<= (Y_END+Y_Start)/2.0+0.5  && x<1.0)   // 8 & 9 are calculated manually, trying to get a easy code. (Top-Bottom)/2
                    {
                        fprintf(fpout, "%le %le %le %le %le %le %le %le %le %le %le 3\n", x, y, z, junk_0, junk_0, junk_0, junk_0, junk_0, junk_0, radius, junk_0);
                        blue++;
                    }

                    else
                    {
                        fprintf(fpout, "%le %le %le %le %le %le %le %le %le %le %le 7\n", x, y, z, junk_0, junk_0, junk_0, junk_0, junk_0, junk_0, radius, junk_0);
                    }
			    }
        }
        
  }
    

    fclose(fp);
    std::cout<<"number of blue particles "<< blue<< std::endl;
  return 0;
}
