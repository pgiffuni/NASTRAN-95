# Declaration of variables
FC = gfortran
FFLAGS =  -w -fno-range-check -fno-automatic -fcray-pointer 
#-std=legacy 

# file names
EXEC=nasthelp.exe
SOURCES = $(wildcard *.f90)
OBJECTS = $(SOURCES:.f90=.o)

# Main target
$(EXEC):$(OBJECTS)
	$(FC) $(OBJECTS) $(FFLAGS) -g -o $(EXEC)
	
# To obtain the object files
%.o:%.f90
	$(FC) $(FFLAGS) -g -c $< -o $@

# To remove generated files

clean:
	rm -f $(EXEC) $(OBJECTS)
