CC = gcc
LD = gcc
CFLAGS = -Wall -O3 -fopt-info-vec -march=native #-fopt-info-vec-missed #-march=native 
#LDFLAGS = -pg
#CFLAGS = -Wall -Ofast -march=haswell -ffast-math -fopt-info-vec #-fopt-info-vec-missed
INCLUDES=-I/opt/X11/include
LIBS = -L/opt/X11/lib -lX11 -lm
RM = /bin/rm -f
OBJS = graphics.o galsim.o
EXEC = galsim

$(EXEC): $(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) $(LIBS) -o $(EXEC)

graphics.o: graphics.c graphics.h
	$(CC) $(CFLAGS) $(INCLUDES) -c graphics.c

galsim.o: galsim.c graphics.h
	$(CC) $(CFLAGS) $(INCLUDES) -c galsim.c

clean:
	$(RM) $(EXEC) $(OBJS)
