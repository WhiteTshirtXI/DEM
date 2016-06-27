
#include <stdexcept>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <cmath>
#include <sstream>
#include <fstream>

int i_0=0,k=0;
double vx=0.0,vy=0.0,vz=0.0;

int particlenumber (int begin_file)
{
    int BIGNUMBER = 90000;
    float junk_0;
    float x=0,y=0,z=0,radius=0, v_x=0;
    int color=0;
    FILE *fp_0;
    
    char particlenumberfile[256];
    
    sprintf (particlenumberfile, "file%04d", begin_file);
    
    if ((fp_0 = fopen (particlenumberfile, "r")) == NULL) {
        printf("Cannot read File!!!!\n");
        exit(1);
    }
    
    for (i_0 = 0; i_0 < BIGNUMBER; i_0++) {
        if ((fscanf(fp_0, "%e %e %e %e %e %e %e %e %e %e %e %d\n", &x, &y, &z, &v_x, &junk_0, &junk_0, &junk_0, &junk_0, &junk_0, &radius, &junk_0, &color) == EOF))
            break;
    }
    fclose (fp_0);
    
    return (i_0);
}

void particle_V_average (int begin_file)
{
    int BIGNUMBER = 90000, i;
    float junk_0;
    float x=0,y=0,z=0,radius=0, v_x=0, v_y=0, v_z=0;
    int color=0;
    FILE *fp_0;
    
    char particlenumberfile[256];
    
    sprintf (particlenumberfile, "file%04d", begin_file);
    
    if ((fp_0 = fopen (particlenumberfile, "r")) == NULL) {
        printf("Cannot read File\n");
        exit(1);
    }
    for (i = 0; i < i_0; i++) {//352 is the number of wall particles
        if ((fscanf(fp_0, "%e %e %e %e %e %e %e %e %e %e %e %d\n", &x, &y, &z, &v_x, &v_y, &v_z, &junk_0, &junk_0, &junk_0, &radius, &junk_0, &color) == EOF))
            break;
        else {
            if (color==7)
            {
            vx+=v_x;
            vy+=v_y;
            vz+=v_z;
            k++;
            }
        
    }
    }
    vx=vx/k;
    vy=vy/k;
    vz=vz/k;
    fclose (fp_0);
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
    double x=0,y=0,z=0,radius=0,junk_0=0, v_x=0,v_y=0,v_z=0;
    double vx_layer, vy_layer, vz_layer;// this only works for each slides particles number <400 need a notice sign
    int color=0;
    int blue=0;
    FILE *fp, *fpout_1;
    int i;
    //int begin_file=0001;
    
    
    for (int i=1; i<argc;i++) {  //file is the file number we are going to analyze.
        if (std::string(argv[i])=="-f") {
            i++;  //why needs i++ here?
            begin_file=atoi(argv[i]);
            std::cout<<"begin file number "<<begin_file<<"  output:Granular_Temperature "<<std::endl;
        }
        
    }
    std::string usage= std::string("usage: -f is the file number we are going to analyze.\n Get granular temperature on each slides on different y\n");
    if (argc!=3) {
        std::cout << usage << std::endl;
    }
    
    number_of_particles = particlenumber(begin_file);
    particle_V_average(begin_file);
    
    std::cout<<"Vx "<<vx<<"  Vy "<<vy<<"  Vz  "<<vz<<" k= "<<k<<std::endl;
    
    FILE *fp_0;
    char particlenumberfile[256];
    
    sprintf (particlenumberfile, "file%04d", begin_file);
    sprintf(outfile_1, "Granular_Temperature");
    if ((fpout_1 = fopen(outfile_1, "w")) == NULL) {
        printf("Cannot Open File\n");
        exit(1);
    }
    fprintf(fpout_1, "layer(y) d_vx d_vy d_vz T_g \n");
    
    for (int x_layer= 0; x_layer<15; x_layer++)
    {
    int i_layer=0; // number of particle in this layer
        if ((fp_0 = fopen (particlenumberfile, "r")) == NULL) {
            printf("Cannot read File\n");
            exit(1);
        }
    {
    for (int i = 0; i < i_0; i++) {
        if ((fscanf(fp_0, "%le %le %le %le %le %le %le %le %le %le %le %d\n", &x, &y, &z, &v_x, &v_y, &v_z, &junk_0, &junk_0, &junk_0, &radius, &junk_0, &color) == EOF))
            break;
        else {
            if (y>=x_layer && y<(x_layer+1) && color ==7)
            vx_layer+=v_x;
            vy_layer+=v_y;
            vz_layer+=v_z;
            i_layer++;
        }
    }
    fclose (fp_0);
    vx_layer=vx_layer/i_layer;
    vy_layer=vy_layer/i_layer;
    vz_layer=vz_layer/i_layer;
    fprintf(fpout_1, "%i %e %e %e %e\n", x_layer, (vx_layer-vx)*(vx_layer-vx), (vy_layer-vy)*(vy_layer-vy),(vz_layer-vz)*(vz_layer-vz), ((vx_layer-vx)*(vx_layer-vx)+(vy_layer-vy)*(vy_layer-vy)+ (vz_layer-vz)*(vz_layer-vz))/3.0); //velocity is averaged here
    }
    }
    
    return 0;
}
