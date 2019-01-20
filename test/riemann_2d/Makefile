HOMEDIR=../../

include $(HOMEDIR)Makefile.global

TOOLS_OBJ=$(HOMEDIR)$(TOOLS_DIR)/obj/
CODE_OBJ=$(HOMEDIR)obj/
INCLUDE += -I$(TOOLS_OBJ) -I$(CODE_OBJ)
SOL2CSV =$(HOMEDIR)tools/bin/sol2csv 1
SOL2PARA =$(HOMEDIR)tools/bin/sol2para

all: compile run
compile:
	@$(MAKE) -C $(HOMEDIR) --no-print-directory

gridgen: gridgen.F90 $(CODE_OBJ)/types.o $(CODE_OBJ)/file_io.o
	$(FC) $(FFLAGS) $^ -o $@ $(LIBS) $(INCLUDE)

run: gridgen
	@./$<
	@$(HOMEDIR)bin/adapt
	$(SOL2PARA)


clean:
	@rm -f *.dat *.vtk gridgen *.mod fort.* *.csv 

test:
	for i in sol_*.dat; do $(SOL2PARA) $$i $${i%dat}vtk; done
	for i in sol_*.dat; do $(SOL2CSV) $$i $${i%dat}csv; done
