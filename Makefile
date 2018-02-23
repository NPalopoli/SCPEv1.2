#Makefile to compile SCPE
#January 24 2002
# V 1.0


OBJS = scpe.o random.o diag.o
CC = gcc
LFLAGS = -lm 
#LFLAGS = -lm  -pg  uso con profile
#LFLAGS = -lm  -g  #debugging
CFLAGS = -O3 -fomit-frame-pointer

scpe.exe: $(OBJS)
	$(CC) $(CFLAGS) -o scpe.exe $(OBJS) $(LFLAGS)


clean: rm *.o 

