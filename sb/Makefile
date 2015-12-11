# Declaration of variables
FC = gfortran
FFLAGS =  -w -fno-range-check -fno-automatic -fcray-pointer -std=legacy 

# file names
EXEC=nastran.exe
SOURCES = $(wildcard *.f)
OBJECTS = $(SOURCES:.f=.o)

# Main target
$(EXEC):$(OBJECTS)
	$(FC) $(OBJECTS) $(FFLAGS) -g -o $(EXEC)
	
# To obtain the object files
%.o:%.f
	$(FC) $(FFLAGS) -g -c $< -o $@

# To remove generated files

clean:
	rm -f $(EXEC) $(OBJECTS)
