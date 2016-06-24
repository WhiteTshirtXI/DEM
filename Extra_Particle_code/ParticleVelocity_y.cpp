
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
    int ft;
    char           infile[256];
    char           outfile_1[256];
    char           outfile_2[256];
    char           outfile_3[256];
    char           outfile_4[256];
    char           outfile_5[256];
    
    int begin_file=0000;
    int end_file=0000;
    int number_of_particles=0000;
    double x=0,y=0,z=0,radius=0,junk_0=0, v_y;
    int color=0;
    int blue=0;
    FILE *fp, *fpout_1, *fpout_2, *fpout_3, *fpout_4, *fpout_5;
    int i;
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
    std::string usage= std::string("usage: -b to draw picture of a file number\n -e for the final file number\n printout radius y v_y");
    if (argc!=3) {
        std::cout << usage << std::endl;
    }
    std::cout<<"flip_time"<<std::endl;
    std::cin>>ft;
    
    //if using -b 1, begin_file will be file0001
    
    
    number_of_particles = particlenumber(begin_file);
    
    //std::cout<<"number_of_particles "<<number_of_particles<<std::endl;
    sprintf(outfile_1, "blue1_velocity_%d",ft);
    if ((fpout_1 = fopen(outfile_1, "w")) == NULL) {
        printf("Cannot Open File 290\n");
        exit(1);
    }
    
    sprintf(outfile_2, "blue2_velocity_%d",ft);
    if ((fpout_2 = fopen(outfile_2, "w")) == NULL) {
        printf("Cannot Open File 470\n");
        exit(1);
    }
    
    sprintf(outfile_3, "blue3_velocity_%d",ft);
    if ((fpout_3 = fopen(outfile_3, "w")) == NULL) {
        printf("Cannot Open File 650\n");
        exit(1);
    }
    
    sprintf(outfile_4, "blue4_velocity_%d",ft);
    if ((fpout_4 = fopen(outfile_4, "w")) == NULL) {
        printf("Cannot Open File 830\n");
        exit(1);
    }
    
    sprintf(outfile_5, "blue5_velocity_%d",ft);
    if ((fpout_5 = fopen(outfile_5, "w")) == NULL) {
        printf("Cannot Open File 1010\n");
        exit(1);
    }
				
    
    
    
    for (int g=begin_file; g<end_file; g++)
    {
        std::cout<<g<<std::endl;
        sprintf (infile, "file%04d", g);
        if((fp=fopen(infile, "r"))==NULL)
        {
            std::cout<<"no file number "<<g<<std::endl;
            exit(1);
        }
        
        for (i=0;i<number_of_particles;i++)
        {
            fscanf(fp, "%le %le %le %le %le %le %le %le %le %le %le %d\n", &x, &y, &z, &junk_0, &v_y, &junk_0, &junk_0, &junk_0, &junk_0, &radius, &junk_0, &color);
            
            {
                if (color==3 && blue==0)
                {
                    
                    fprintf(fpout_1, "%e %e %e\n", radius, y, v_y);
                    blue++;
                }
                
                else if (color==3 && blue==1)
                {
                    
                    fprintf(fpout_2, "%e %e %e\n", radius, y, v_y);
                    blue++;
                }
                
                else if (color==3 && blue==2)
                {
                    
                    fprintf(fpout_3, "%e %e %e\n", radius, y, v_y);
                    blue++;
                }
                
                else if (color==3 && blue==3)
                {
                    
                    fprintf(fpout_4, "%e %e %e\n", radius, y, v_y);
                    blue++;
                }
                
                else if (color==3 && blue==4)
                {
                    
                    fprintf(fpout_5, "%e %e %e\n", radius, y, v_y);
                    blue=0;
                }
            }
            
            
        }
        
        fclose(fp);
        
    }
    fclose(fpout_1);
    fclose(fpout_2);
    fclose(fpout_3);
    fclose(fpout_4);
    fclose(fpout_5);
    
    return 0;
}
