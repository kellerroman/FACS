HOMEDIR=../../

include $(HOMEDIR)Makefile.global

TOOLS_OBJ=$(HOMEDIR)$(TOOLS_DIR)/obj/
CODE_OBJ=$(HOMEDIR)obj/
INCLUDE += -I$(TOOLS_OBJ) -I$(CODE_OBJ)

all: compile run post

compile:
	@$(MAKE) -C $(HOMEDIR) --no-print-directory

gridgen: gridgen.F90 $(CODE_OBJ)/types.o $(CODE_OBJ)/file_io.o
	$(FC) $(FFLAGS) $^ -o $@ $(LIBS) $(INCLUDE)

run: gridgen
	@./$<
	@cp sol.dat sol_00000.dat
	@$(HOMEDIR)/tools/bin/sol2csv
	@mv sol.csv start.csv
	@$(HOMEDIR)/bin/adapt
post: exact.csv
	@$(HOMEDIR)/tools/bin/sol2csv
	@gnuplot -p plot.gnuplot

clean:
	@rm -f *.dat *.vtk gridgen *.mod fort.* *.csv exact_riemann

exact.csv: exact_riemann
	@./$<

exact_riemann: exact_riemann.F90
	gfortran $< -o $@

test:
	for i in sol_*.dat; do ../../tools/bin/sol2para $$i $${i%dat}vtk; done
