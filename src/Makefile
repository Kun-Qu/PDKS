CPP  =mpiCC 
CC   =mpicc
RES  = 
OBJ  = grid_2d.o explicit.o bc.o lbm.o dks.o parallel.o $(RES)
LINKOBJ  = grid_2d.o explicit.o bc.o lbm.o dks.o parallel.o $(RES)
BIN  = Pdks
CXXINCS = -I/home/coolwater/lib/extctrl/include
CPUFLAGS =
CXXFLAGS = $(CXXINCS) $(CPUFLAGS) -O3  -DLAMINAR #-DINVISCID 
CFLAGS = $(INCS)   
CXXLIBS = -lextctrl -L/home/coolwater/lib/extctrl
LIBS = $(CXXLIBS)
RM = rm -f

.PHONY: all all-before all-after clean clean-custom

all: all-before $(BIN) all-after


clean: clean-custom
	${RM} $(OBJ) $(BIN)

$(BIN): $(OBJ)
	$(CPP) $(LINKOBJ) -o $(BIN) $(LIBS)

grid_2d.o: grid_2d.cpp
	$(CPP) -c grid_2d.cpp -o grid_2d.o $(CXXFLAGS)

explicit.o: explicit.cpp
	$(CPP) -c explicit.cpp -o explicit.o $(CXXFLAGS)

bc.o: bc.cpp
	$(CPP) -c bc.cpp -o bc.o $(CXXFLAGS)

lbm.o: lbm.cpp
	$(CPP) -c lbm.cpp -o lbm.o $(CXXFLAGS)

dks.o: dks.cpp
	$(CPP) -c dks.cpp -o dks.o $(CXXFLAGS)
parallel.o: parallel.cpp
	 $(CPP) -c parallel.cpp -o parallel.o $(CXXFLAGS)
