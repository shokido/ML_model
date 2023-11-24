# Make file for compiling JCOPE-T code
FC = gfortran
FFLAGS = -O2
CFLAGS =

TARGET_NN=exec_1d_nn.out
all: $(TARGET_NN)
OBJ_CODES_NN= ml_param.o \
	ml_utils.o \
	nnf_sub.o \
	solve_diag.o \
	do_mlmodel_nn.o  

$(TARGET_NN):$(OBJ_CODES_NN) 
	$(FC) -o $(TARGET_NN) $(FFLAGS) $(OBJ_CODES_NN)

%.o: %.f90
	$(FC) -c $(FFLAGS) $<
clean:
	rm -f *.o *mod $(TARGET_NN) 
