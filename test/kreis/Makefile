HOMEDIR=../../

include $(HOMEDIR)Makefile.global

TOOLS_OBJ=$(HOMEDIR)$(TOOLS_DIR)/obj/
CODE_OBJ=$(HOMEDIR)obj/
INCLUDE += -I$(TOOLS_OBJ) -I$(CODE_OBJ)

all: compile run test

compile:
	@$(MAKE) -C $(HOMEDIR) --no-print-directory

gridgen: gridgen.F90 $(CODE_OBJ)/types.o $(CODE_OBJ)/file_io.o
	$(FC) $(FFLAGS) $^ -o $@ $(LIBS) $(INCLUDE)

run: gridgen
	@./$<
	@time $(HOMEDIR)/bin/adapt

post: run
	@$(HOMEDIR)/tools/bin/sol2para

clean:
	@rm -f *.dat *.vtk gridgen *.mod fort.*
test:
	@mv sol.dat sol_00000.dat
	for i in sol_*.dat; do ../../tools/bin/sol2para $$i $${i%dat}vtk; done
