CC = cc
#CFLAGS = -I /usr/X11/include/ -I /usr/local/include

# if it complains about 64 bit architecture, use -m32
CFLAGS = -I /usr/X11/include/ -I /usr/local/include

#LIBS = -L/usr/X11R6/lib -L/usr/local/lib -lX11 -lglut -lGL -lGLU -ltiff -lm -Wl,-dylib_file,/System/Library/Frameworks/OpenGL.framework/Versions/A/Libraries/libGL.dylib:/System/Library/Frameworks/OpenGL.framework/Versions/A/Libraries/libGL.dylib

# if it complains about 64 bit architecture, use -m32
LIBS = -L/usr/X11R6/lib -L/usr/local/lib -lX11 -lglut -lGL -lGLU -ltiff -lm -Wl,-dylib_file,/System/Library/Frameworks/OpenGL.framework/Versions/A/Libraries/libGL.dylib:/System/Library/Frameworks/OpenGL.framework/Versions/A/Libraries/libGL.dylib


movie_view : main.o 
	$(CC) -o movie_view main.o $(LIBS)
	
main.o :  main.c
	$(CC) -c $(CFLAGS) main.c
	
clean:  
	rm -f movie_view *.o
