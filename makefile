 CC = g++
 CFLAGS  = -g -Wall
 main = src/Work-distribution_peak_finder
 bin = bin_Linux_x86_64/Work-distribution_peak_finder__Linux_x86_64
 util = src/functions

# Requires GSL to compile
 LIB = -lgsl -lgslcblas -lm



$(bin) :	Work-distribution_peak_finder.o functions.o
	$(CC) $(CFLAGS)	-o	$(bin)	Work-distribution_peak_finder.o	functions.o $(LIB)

Work-distribution_peak_finder.o :  $(main).cpp
	$(CC) $(CFLAGS) -c $(main).cpp

functions.o : $(util).cpp  $(util).h 
	$(CC) $(CFLAGS) -c $(util).cpp

clean: 
	$(RM) count *.o *~
