HOMEDIR=../../../

include $(HOMEDIR)Makefile.global

TOOLS_OBJ=$(HOMEDIR)$(TOOLS_DIR)/obj/
CODE_OBJ=$(HOMEDIR)obj/
INCLUDE += -I$(TOOLS_OBJ) -I$(CODE_OBJ)
SOL2CSV =$(HOMEDIR)tools/bin/sol2csv 2
SOL2PARA =$(HOMEDIR)tools/bin/sol2para

para: compile run test
post: compile run plot

compile:
	@$(MAKE) -C $(HOMEDIR) --no-print-directory

gridgen: gridgen.F90 $(CODE_OBJ)/types.o $(CODE_OBJ)/file_io.o
	$(FC) $(FFLAGS) $^ -o $@ $(LIBS) $(INCLUDE)

run: gridgen
	@./$<
	@cp sol.dat sol_00000.dat
	@$(SOL2CSV)
	@mv sol.csv start.csv
	@$(HOMEDIR)bin/adapt
plot: exact.csv
	@$(SOL2CSV)
	@gnuplot -p plot.gnuplot

clean:
	@rm -f *.dat *.vtk gridgen *.mod fort.* *.csv 

exact.csv: ../exact_riemann
	@$<

../exact_riemann: ../exact_riemann.F90
	gfortran $< -o $@

test:
	for i in sol_*.dat; do $(SOL2PARA) $$i $${i%dat}vtk; done
	for i in sol_*.dat; do $(SOL2CSV) $$i $${i%dat}csv; done
