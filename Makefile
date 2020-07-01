MPIF90 ?= mpif90
SRCDIR = src
OBJDIR = obj
BINDIR = bin
FC = $(MPIF90)
FFLAGS = -O3 -fopenmp -fopt-info-all \
	-mcmodel=medium -fconvert=big-endian -g \
	-fbacktrace -fbounds-check -I$(OBJDIR) -J$(OBJDIR)

LEC: laser.o LEC.o
	@mkdir -p $(OBJDIR)
	@mkdir -p $(BINDIR)
	$(FC) $(FFLAGS) $(OBJDIR)/laser.o $(OBJDIR)/LEC.o -o $(BINDIR)/LEC

LEC.o: $(SRCDIR)/LEC.f
	$(FC) $(FFLAGS) -c -o $(OBJDIR)/LEC.o $(SRCDIR)/LEC.f 

laser.o: $(SRCDIR)/laser.f
	$(FC) $(FFLAGS) -c -o $(OBJDIR)/laser.o $(SRCDIR)/laser.f

clean:
	rm $(OBJDIR)/*.o
	rm -r $(BINDIR)
