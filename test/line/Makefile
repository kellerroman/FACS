HOMEDIR=../../

include $(HOMEDIR)Makefile.global

TOOLS_OBJ=$(HOMEDIR)$(TOOLS_DIR)/obj/
CODE_OBJ=$(HOMEDIR)obj/
INCLUDE += -I$(TOOLS_OBJ) -I$(CODE_OBJ)

all: compile run

compile:
	@$(MAKE) -C $(HOMEDIR) --no-print-directory

gridgen: gridgen.F90 $(CODE_OBJ)/types.o $(CODE_OBJ)/file_io.o
	$(FC) $(FFLAGS) $^ -o $@ $(LIBS) $(INCLUDE)

run: gridgen
	@./$<
	@$(HOMEDIR)/tools/bin/sol2csv
	@mv sol.csv old.csv
	@$(HOMEDIR)/bin/adapt
	@$(HOMEDIR)/tools/bin/sol2csv
	@$(HOMEDIR)/tools/bin/sol2para
	@gnuplot -p plot.gnuplot

clean:
	rm -f *.dat *.vtk gridgen *.mod *.csv
