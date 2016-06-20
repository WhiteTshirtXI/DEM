#include <stdio.h>
#include <stdlib.h>
#include <math.h>   
#include <unistd.h>   
#include <GL/glut.h>
#include <GL/glx.h>
#include <tiffio.h>     /* Sam Leffler's libtiff library. */

#ifndef M_PI
#define M_PI 3.14159265
#endif

int centered=1;
int thermal=0;
int circle=0;
int tagged=0;
int size=800;
int *c, n, l=0, picture_end=0;
int increment=1, offset=0;
double circle_radius=6.0, circle_x=27.0, circle_y=46.0;
double scaling, xscaling, yscaling, zscaling;
double max_d, max_d_test;
double max_x;
double max_y;
double max_z;

static int attributeList[] = { GLX_RGBA, GLX_ALPHA_SIZE, 4, GLX_DEPTH_SIZE, 4, None };
static Display* dpy;
static GC xgc;
static Pixmap pixmap;
static GLXPixmap glxPixmap;
static Window win;
XVisualInfo *vi;
XSetWindowAttributes swa;
GLXContext cx;
    
char command[50],input[10], *base;

GLdouble *x, *y, *z, *r, *temperature;

GLfloat shin = 50.0, ishin = 0.0;

GLfloat mat_amb_yellow[] = { 0.40, 0.35, 0.02, 1.0 };
GLfloat mat_dif_yellow[] = { 0.75, 0.71, 0.03, 1.0 };
GLfloat mat_spe_yellow[] = { 0.63, 0.66, 0.17, 1.0 };

GLfloat mat_amb_green[] = { 0.02, 0.12, 0.02, 1.0 };
GLfloat mat_dif_green[] = { 0.08, 0.61, 0.08, 1.0 };
GLfloat mat_spe_green[] = { 0.63, 0.73, 0.63, 1.0 };
   
GLfloat mat_amb_bluegreen[] = { 0.02, 0.08, 0.04, 1.0 };
GLfloat mat_dif_bluegreen[] = { 0.08, 0.41, 0.15, 1.0 };
GLfloat mat_spe_bluegreen[] = { 0.63, 0.73, 0.63, 1.0 };
   
GLfloat mat_amb_blue[] = { 0.02, 0.02, 0.62, 1.0 };
GLfloat mat_dif_blue[] = { 0.08, 0.08, 0.61, 1.0 };
GLfloat mat_spe_blue[] = { 0.63, 0.63, 0.73, 1.0 };
   
GLfloat mat_amb_darkblue[] = { 0.0125, 0.0125, 0.40, 1.0 };
GLfloat mat_dif_darkblue[] = { 0.05, 0.05, 0.40, 1.0 };
GLfloat mat_spe_darkblue[] = { 0.63, 0.63, 0.73, 1.0 };
   
GLfloat mat_amb_red[] = { 0.17, 0.01, 0.01, 1.0 };
GLfloat mat_dif_red[] = { 0.61, 0.04, 0.04, 1.0 };
GLfloat mat_spe_red[] = { 0.73, 0.63, 0.63, 1.0 };

GLfloat mat_amb_orange[] = { 0.25, 0.18, 0.01, 1.0 };
GLfloat mat_dif_orange[] = { 0.67, 0.50, 0.04, 1.0 };
GLfloat mat_spe_orange[] = { 0.68, 0.64, 0.30, 1.0 };  

GLfloat mat_amb_darkorange[] = { 0.20, 0.08, 0.01, 1.0 };
GLfloat mat_dif_darkorange[] = { 0.67, 0.25, 0.04, 1.0 };
GLfloat mat_spe_darkorange[] = { 0.68, 0.64, 0.30, 1.0 };

GLfloat mat_amb_grey[] = { 0.02, 0.02, 0.02, 1.0 };
GLfloat mat_dif_grey[] = { 0.08, 0.08, 0.08, 1.0 };
GLfloat mat_spe_grey[] = { 0.63, 0.63, 0.63, 1.0 };
   
GLfloat mat_amb_mgrey[] = { 0.07, 0.07, 0.07, 0.0 };
GLfloat mat_spe_mgrey[] = { 0.10, 0.10, 0.10, 0.0 };
GLfloat mat_dif_mgrey[] = { 0.68, 0.68, 0.68, 0.0 };
           
GLfloat mat_amb_igrey[] = { 0.91, 0.91, 0.91, 0.0 };
GLfloat mat_spe_igrey[] = { 0.93, 0.93, 0.93, 0.0 };
GLfloat mat_dif_igrey[] = { 0.93, 0.93, 0.93, 0.0 };
           
void usage()
{
  printf("\n");
  printf("usage:  movie_view -[nbeitsworXYCTh] \n\
\t movie viewer for PD data files.\n\
\n\
-n\t number of particles            \n\
-B\t file name base                \n\
-b\t beginning file number                \n\
-e\t ending file number                \n\
-i\t increment file number by...               \n\
-c\t do NOT center the image (i.e., if it is already centered)    \n\
-t (int)\t tag particles for mixing visualization using reference file number (int).               \n\
-s\t define size of picture to generate (square)               \n\
-T\t use temeratures to determine color               \n\
-C\t tag a blob of color for mixing               \n\
-X\t (double) x position of blob              \n\
-Y\t (double) y position of blob             \n\
-r\t (double) radius of blob              \n\
-w\t define scaling width for zooming image               \n\
-x\t define x scaling width for shiftinging image               \n\
-y\t define y scaling width for shifting image               \n\
-z\t define z scaling width for shifting image               \n\
\t \t this info is printed at the onset of image gernation and can \n\
\t \t be used to continue image capture \n\
-o\t offset the saved file number by...                \n\
-h\t show this message              \n\n");

}				
 
int
writetiff(char *filename, char *description,
  int x, int y, int width, int height, int compression)
{
  TIFF *file;
  GLubyte *image, *p;
  int i;

  file = TIFFOpen(filename, "w");
  if (file == NULL) {
    return 1;
  }
  image = (GLubyte *) malloc(width * height * sizeof(GLubyte) * 3);

  /* OpenGL's default 4 byte pack alignment would leave extra bytes at the
     end of each image row so that each full row contained a number of bytes
     divisible by 4.  Ie, an RGB row with 3 pixels and 8-bit componets would
     be laid out like "RGBRGBRGBxxx" where the last three "xxx" bytes exist
     just to pad the row out to 12 bytes (12 is divisible by 4). To make sure
     the rows are packed as tight as possible (no row padding), set the pack
     alignment to 1. */
  glPixelStorei(GL_PACK_ALIGNMENT, 1);

  glReadPixels(x, y, width, height, GL_RGB, GL_UNSIGNED_BYTE, image);
  TIFFSetField(file, TIFFTAG_IMAGEWIDTH, (uint32) width);
  TIFFSetField(file, TIFFTAG_IMAGELENGTH, (uint32) height);
  TIFFSetField(file, TIFFTAG_BITSPERSAMPLE, 8);
  TIFFSetField(file, TIFFTAG_COMPRESSION, compression);
  TIFFSetField(file, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
  TIFFSetField(file, TIFFTAG_SAMPLESPERPIXEL, 3);
  TIFFSetField(file, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
  TIFFSetField(file, TIFFTAG_ROWSPERSTRIP, 1);
  TIFFSetField(file, TIFFTAG_IMAGEDESCRIPTION, description);
  p = image;
  for (i = height - 1; i >= 0; i--) {
    if (TIFFWriteScanline(file, p, i, 0) < 0) {
      free(image);
      TIFFClose(file);
      return 1;
    }
    p += width * sizeof(GLubyte) * 3;
  }
  free(image);
  TIFFClose(file);
  return 0;
}
  
void init(void) 
{  
   GLfloat light_position[] = {0.5, 1.0, 1.0, 0.0};
   
   glMatrixMode(GL_PROJECTION);  
   
   glLightfv(GL_LIGHT0, GL_POSITION, light_position);
    
   glEnable(GL_LIGHT0);
   glEnable(GL_LIGHTING);
   glEnable(GL_DEPTH_TEST);
   glEnable(GL_ALPHA_TEST);

   gluPerspective( /* field of view in degree */ 30.0,
                   /* aspect ratio */ 1.0,
                   /* Z near */ 0.01, 
                   /* Z far */ 500.0);
		   
   gluLookAt(0.25, 0.25, 3.0,  /* eye is at (0,0,3) */
             0.0, 0.0, 0.0,  /* center is at (0,0,0) */
             0.0, 1.0, 0.0);  /* up is in positive Y direction */	     
   
   glClearColor (1.0, 1.0, 1.0, 0.0);   
}

void make_sphere(GLint gc, GLdouble gr, GLdouble gx, GLdouble gy, GLdouble gz)
{
   glPushMatrix();
   if (gc==0) {
     
	    glMaterialfv(GL_FRONT, GL_AMBIENT, mat_amb_grey);
	    glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_dif_grey);
	    glMaterialfv(GL_FRONT, GL_SPECULAR, mat_spe_grey);
	    glMaterialf(GL_FRONT, GL_SHININESS, shin);
       
     }
     
   else if (gc==2){
     
	    glMaterialfv(GL_FRONT, GL_AMBIENT, mat_amb_mgrey);
	    glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_dif_mgrey);
	    glMaterialfv(GL_FRONT, GL_SPECULAR, mat_spe_mgrey);
	    glMaterialf(GL_FRONT, GL_SHININESS, shin);
	
     } 
		  
   else if (gc==1){
     
	    glMaterialfv(GL_FRONT, GL_AMBIENT, mat_amb_igrey);
	    glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_dif_igrey);
	    glMaterialfv(GL_FRONT, GL_SPECULAR, mat_spe_igrey);
	    glMaterialf(GL_FRONT, GL_SHININESS, shin);
	
     } 
	
   else if (gc==6){
     
	    glMaterialfv(GL_FRONT, GL_AMBIENT, mat_amb_darkblue);
	    glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_dif_darkblue);
	    glMaterialfv(GL_FRONT, GL_SPECULAR, mat_spe_darkblue);
	    glMaterialf(GL_FRONT, GL_SHININESS, shin);
	
     } 
			  
   else if (gc==4){
     
	    glMaterialfv(GL_FRONT, GL_AMBIENT, mat_amb_blue);
	    glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_dif_blue);
	    glMaterialfv(GL_FRONT, GL_SPECULAR, mat_spe_blue);
	    glMaterialf(GL_FRONT, GL_SHININESS, shin);
	
     } 
			  
   else if (gc==5){
     
	    glMaterialfv(GL_FRONT, GL_AMBIENT, mat_amb_bluegreen);
	    glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_dif_bluegreen);
	    glMaterialfv(GL_FRONT, GL_SPECULAR, mat_spe_bluegreen);
	    glMaterialf(GL_FRONT, GL_SHININESS, shin);
	
     } 
     
   else if (gc==3){
     
	    glMaterialfv(GL_FRONT, GL_AMBIENT, mat_amb_green);
	    glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_dif_green);
	    glMaterialfv(GL_FRONT, GL_SPECULAR, mat_spe_green);
	    glMaterialf(GL_FRONT, GL_SHININESS, shin);
	
     } 
		  
   else if (gc==7){
     
	    glMaterialfv(GL_FRONT, GL_AMBIENT, mat_amb_yellow);
	    glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_dif_yellow);
	    glMaterialfv(GL_FRONT, GL_SPECULAR, mat_spe_yellow);
	    glMaterialf(GL_FRONT, GL_SHININESS, shin);
	
     } 
			  
     
   else if (gc==8){
     
	    glMaterialfv(GL_FRONT, GL_AMBIENT, mat_amb_orange);
	    glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_dif_orange);
	    glMaterialfv(GL_FRONT, GL_SPECULAR, mat_spe_orange);
	    glMaterialf(GL_FRONT, GL_SHININESS, shin);
	
     } 
		  
   else if (gc==9){
     
	    glMaterialfv(GL_FRONT, GL_AMBIENT, mat_amb_darkorange);
	    glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_dif_darkorange);
	    glMaterialfv(GL_FRONT, GL_SPECULAR, mat_spe_darkorange);
	    glMaterialf(GL_FRONT, GL_SHININESS, shin);
	
     } 
		  
   else if (gc==10){
     
	    glMaterialfv(GL_FRONT, GL_AMBIENT, mat_amb_red);
	    glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_dif_red);
	    glMaterialfv(GL_FRONT, GL_SPECULAR, mat_spe_red);
	    glMaterialf(GL_FRONT, GL_SHININESS, shin);
	    	
     } 
		  
   else {
        
	    glMaterialfv(GL_FRONT, GL_AMBIENT, mat_amb_igrey);
	    glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_dif_igrey);
	    glMaterialfv(GL_FRONT, GL_SPECULAR, mat_spe_igrey);
	    glMaterialf(GL_FRONT, GL_SHININESS, shin);
	    
     } 
    
   glTranslatef(gx, gy, gz);
   glutSolidSphere (gr, 30, 30);
   glPopMatrix();
}

int get_color(double temperature) 
{  
   if (3.0*temperature<=0.1) return(0);
   else if (3.0*temperature<=0.2) return(1);
   else if (3.0*temperature<=0.3) return(2);
   else if (3.0*temperature<=0.4) return(3);
   else if (3.0*temperature<=0.5) return(4);
   else if (3.0*temperature<=0.6) return(5);
   else if (3.0*temperature<=0.7) return(6);
   else if (3.0*temperature<=0.8) return(7);
   else if (3.0*temperature<=0.9) return(8);
   else  return(9);
   
}

void draw_image(void)
{
   int i, wx, wy, ww, wh, color, ijunk;
   double junk;
   char newfile[256];
	 
	 FILE *fp;
    
   sprintf (input, "%s%04d", base, l);
       
   if((fp=fopen(input, "r"))==NULL) {
     printf("Cannot Open File\n");
     exit(1);
   }
      
   for (i=0;i<n;i++) {
     if (!tagged && !circle) fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d", &x[i], &y[i], &z[i], &junk,
		     &junk, &junk, &junk, &junk, &junk, &r[i], &temperature[i], &c[i]);
     else fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d", &x[i], &y[i], &z[i], &junk,
		     &junk, &junk, &junk, &junk, &junk, &r[i], &temperature[i],
		     &ijunk);
   }
   
   fclose(fp);
   printf("getting ready to draw\n");
   glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
   
   for (i=0;i<n;i++) {
      int scolor;
      
      if(thermal) scolor=get_color(temperature[i]);
      else scolor=c[i];
      
      if (centered) make_sphere(scolor,r[i]/scaling, (x[i]-0.5*xscaling)/scaling,
      		(y[i]-0.5*yscaling)/scaling, (z[i]-0.5*zscaling)/scaling);
      else make_sphere(scolor,r[i]/scaling, x[i]/scaling, y[i]/scaling, z[i]/scaling);
   }
  
   printf("finished drawing\n");
   XCopyArea(dpy, pixmap, win, xgc, 0, 0, size, size, 0, 0);
   
   sprintf(newfile, "file%04d.tif",(l+offset));	
	 writetiff(newfile, "Stills for the movie", 0, 0, size, size, COMPRESSION_NONE);
	 
   sprintf(command, "convert file%04d.tif file%04d.jpg\n",l+offset,l+offset);
   system (command);
   sprintf(command, "rm file%04d.tif\n",l+offset);
   system (command);
   l+=increment;
 
}

int main(int argc, char** argv)
{  
   int reference=0;
   int i, option_letter=0, usage_flag=0;
   double junk;
   extern char *optarg;
   FILE *fp; 

   base="file";
    
/**   parse commmand line options   **/

   scaling=0.0;
   xscaling=0.0;
   yscaling=0.0;
   zscaling=0.0;
   
   while ((option_letter = getopt(argc, argv, "n:B:b:e:i:o:s:x:y:z:w:t:r:X:Y:TcC:h")) != -1 ){
  
    if (option_letter == 'n') {    
      n = atoi(optarg);
      usage_flag++;      
    }    
    else if (option_letter == 'B') base = optarg;
    else if (option_letter == 'b') l = atoi(optarg);
    else if (option_letter == 'e') picture_end = atoi(optarg);      
    else if (option_letter == 'i') increment = atoi(optarg);      
    else if (option_letter == 'o') offset = atoi(optarg);      
    else if (option_letter == 's') size = atoi(optarg);      
    else if (option_letter == 'w') scaling = atof(optarg);      
    else if (option_letter == 'x') xscaling = atof(optarg);      
    else if (option_letter == 'y') yscaling = atof(optarg);      
    else if (option_letter == 'z') zscaling = atof(optarg);      
    else if (option_letter == 'r') circle_radius = atof(optarg);      
    else if (option_letter == 'X') circle_x = atof(optarg);      
    else if (option_letter == 'Y') circle_y = atof(optarg);      
    else if (option_letter == 't') {
    	tagged = 1;
	reference = atoi(optarg); 
    }      
    else if (option_letter == 'C') {
    	circle = 1;  
	reference = atoi(optarg); 
    }          
    else if (option_letter == 'T') thermal = 1;      
    else if (option_letter == 'c') centered = 0;      
    else if (option_letter == 'h') {    
       usage();
       exit(1);      
     }  	    
   }   
   if (usage_flag < 1) {   
     usage();
     exit(1);      
   }

   c=(int *)malloc(n*sizeof(int));
   x=(double *)malloc(n*sizeof(double));
   y=(double *)malloc(n*sizeof(double));
   z=(double *)malloc(n*sizeof(double));
   r=(double *)malloc(n*sizeof(double));
   temperature=(double *)malloc(n*sizeof(double));
   
   if (!tagged && !circle) reference=l;
   sprintf (input, "%s%04d", base, reference);
       
   if((fp=fopen(input, "r"))==NULL) {
     printf("Cannot Open File %s\n", input);
     exit(1);
   }
      
   for (i=0;i<n;i++) {
     fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d", &x[i], &y[i], &z[i], &junk,
		     &junk, &junk, &junk,
		     &junk, &junk, &r[i], &temperature[i], &c[i]);
     max_d_test=sqrt(x[i]*x[i]+y[i]*y[i]+z[i]*z[i]);
     if(max_d_test>max_d) max_d=max_d_test;
     if(x[i]>max_x) max_x=x[i];
     if(y[i]>max_y) max_y=y[i];
     if(z[i]>max_z) max_z=z[i];
   }
   
   if (tagged)
   	for (i=0;i<n;i++) {
		if(!centered) {
			if(x[i]>0.0) c[i]=1;
			else c[i]=2;
		}
		else {
			if(x[i]>0.5*max_x) c[i]=1;
			else c[i]=2;
		}
	}
   	
   else if (circle)
   	for (i=0;i<n;i++) {
		double pos;
		
		if (c[i]!=0) {
			pos=sqrt((x[i]-circle_x)*(x[i]-circle_x)+(y[i]-circle_y)*(y[i]-circle_y));
			if(pos<circle_radius) c[i]=3;
			else c[i]=1;
		}
	}
   	
   fclose(fp);
   
   if (scaling==0.0) scaling=max_d;
   if (xscaling==0.0) xscaling=max_x;
   if (yscaling==0.0) yscaling=max_y;
   if (zscaling==0.0) zscaling=max_z;
   printf("scaling information\n use -w %lf -x %lf -y %lf -z %lf\n\tto continue these images\n",
    scaling, xscaling, yscaling, zscaling);
   
   /* get a connection */
    dpy = XOpenDisplay(0);
    if (!dpy) {
    	printf("Problem opening display!\n");
    	exit(1);
	}
	 printf("opened display\n");
   
    /* get an appropriate visual */
    vi = glXChooseVisual(dpy, DefaultScreen(dpy), attributeList);
    if (!vi) {
    	printf("Problem opening visuals!\n");
    	exit(1);
	}
 printf("opened visual\n");
   
    /* create a GLX context */
    cx = glXCreateContext(dpy, vi, 0, GL_FALSE);

    /* create a color map */
    swa.colormap = XCreateColormap(dpy, RootWindow(dpy, vi->screen),
                                   vi->visual, AllocNone);
printf("created colormap\n");
 
    /* create an imaginary window */
    swa.border_pixel = 0;
    swa.event_mask = ExposureMask | StructureNotifyMask | KeyPressMask;
    win = XCreateWindow(dpy, RootWindow(dpy, vi->screen), 0, 0, size, size,
   0, vi->depth, InputOutput, vi->visual,
   CWBorderPixel|CWColormap|CWEventMask, &swa);
    XStoreName(dpy, win, "pixmap rendering");
    
    printf("created window\n");
 
    XMapWindow(dpy, win);
    
    /* create the pixmap */
    xgc = XCreateGC(dpy, win, 0, NULL);
    pixmap = XCreatePixmap(dpy, win, size, size, vi->depth);
    glxPixmap = glXCreateGLXPixmap(dpy, vi, pixmap);
printf("created pixmap\n");
 
    /* connect the context to the pixmap, and set up transforms to match */
    glXMakeCurrent(dpy, glxPixmap, cx);
    glutInit(&argc, argv);
    printf("starting init\n");
    init();

    /* draw and save the series of images */
    while(l<=picture_end) draw_image();
   
    return 0;
}
