# Declaration of variables
FC = gfortran
FFLAGS = -w -fno-automatic -fno-range-check -fcray-pointer
# -fdefault-integer-8 -fdefault-real-8  -std=legacy  

# file names
EXEC=nastran.exe
SOURCES = $(wildcard *.f90)
OBJECTS = $(SOURCES:.f90=.o)

# Main target
$(EXEC):$(OBJECTS)
	$(FC) $(OBJECTS) -g -o $(EXEC)
	
# To obtain the object files
%.o:%.f90
	$(FC) $(FFLAGS) -g -c $< -o $@

# To remove generated files

clean:
	rm -f $(EXEC) $(OBJECTS)


