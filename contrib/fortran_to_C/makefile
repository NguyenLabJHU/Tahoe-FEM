# $Id: makefile,v 1.1 2003-03-09 03:59:57 paklein Exp $
.SUFFIXES: .c .cpp .f .F

##################### Linux
#COMP_C  = cc
#COMP_CC = g++
#COMP_F  = g77
#LINK_C  = cc
#LINK_CC = g++
#LINK_F  = g77

##################### DARWIN
COMP_C   = cc
COMP_CC  = g++
COMP_F   = g77
C_FLAGS  = -D__GNU__
CC_FLAGS = -D__GNU__
F_FLAGS  = -D__GNU__

LINK_C     = $(COMP_C) 
LINK_CC    = $(COMP_CC)
LINK_F     = $(COMP_F)
L_FLAGS_C  = -Xlinker -Y -Xlinker 50 -lg2c
L_FLAGS_CC = -Xlinker -Y -Xlinker 50 -lg2c
L_FLAGS_F  = -Xlinker -Y -Xlinker 50

#################### SGI
#COMP_C  = cc
#COMP_CC = CC
#COMP_F  = f77 

# all but main's
F_TEST_OBJ  = main_f.o double_f.o sqrt_f.o sqrt_f_c.o sqrt_c.o
C_TEST_OBJ  = main_c.o double_f.o sqrt_f.o sqrt_f_c.o sqrt_c.o print_f.o
CC_TEST_OBJ = main_cpp.o double_f.o sqrt_f.o sqrt_f_c.o sqrt_c.o print_f.o

all: f_test C_test CC_test

f_test: $(F_TEST_OBJ)
	$(LINK_F) $(F_TEST_OBJ) -o f_test $(L_FLAGS_F)

C_test: $(C_TEST_OBJ)
	$(LINK_C) $(C_TEST_OBJ) -o C_test $(L_FLAGS_C)

CC_test: $(CC_TEST_OBJ)
	$(LINK_CC) $(CC_TEST_OBJ) -o CC_test $(L_FLAGS_CC)

clean:
	rm *.o C_test f_test CC_test

.f.o:
	$(COMP_F) $(F_FLAGS) -c $<

.c.o:
	$(COMP_C) $(C_FLAGS) -c $<

.cpp.o:
	$(COMP_CC) $(CC_FLAGS) -c $<
