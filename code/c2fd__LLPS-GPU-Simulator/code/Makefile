### If these change, you need to run "make clean".
GRID = grid.cuda
STATE = state/phi2


################

NVCC =nvcc
GCC = g++
LD = ld
HDF5_PATH = ~/hdf5



NV_FLAGS = #-arch compute_20	# for m2070
#NV_FLAGS = -G -g

CFLAGS = -I . -I $(STATE)
CFLAGS += -I $(HDF5_PATH)/include

LDFLAGS += -lm -lcufft  #-limf
LDFLAGS += -L $(HDF5_PATH)/lib
LDFLAGS += -lhdf5 -lhdf5_hl -lsz -lz


################


APP = LLPS

OBJECTS = tools.o  misc.o  unary_function.o \
      $(GRID)/my_hdf5.o  $(GRID)/my_hdf5_viz.o parameters.o  \
      $(STATE)/nondim.o  \
      device.o  state/state.o  LLPS.o

GRID_OBJECTS = $(GRID)/grid_params.o  \
      $(GRID)/debug.o  \
      $(GRID)/math.o  \
      $(GRID)/io.o  \
      $(GRID)/bicgstab.o 


all: $(APP)

device.o: $(shell find . -name '*.inl') /dev/null

$(APP) : $(OBJECTS) $(GRID_OBJECTS) timestamp.cpp 
	$(NVCC) $(CFLAGS) $(LDFLAGS) $(NV_FLAGS) $^ -o $@ 
	$(COPY_BINARIES)


# Note: The second -e deletes a spurious line in .dep in nvcc version 3.2 beta
%.o : %.cpp
	$(NVCC) -c $< $(CFLAGS) $(NV_FLAGS) -o $@
	$(NVCC) -M $< $(CFLAGS) $(NV_FLAGS) > $@.dep
	sed -i -e "1 s:^:$(@D)/:" -e "s:.*/code/LLPS// .*:\\\\:" $@.dep

%.o : %.cu
	$(NVCC) -c $< $(CFLAGS) $(NV_FLAGS) -o $@
	$(NVCC) -M $< $(CFLAGS) $(NV_FLAGS) > $@.dep
	sed -i -e "1 s:^:$(@D)/:" -e "s:.*/code/LLPS// .*:\\\\:" $@.dep


#
# Pick up generated dependency files, and add /dev/null because gmake
# does not consider an empty list to be a list:
#
include $(shell find . -name '*.dep') /dev/null



################


clean:
	$(RM) $(APP) $(OBJECTS) $(GRID_OBJECTS) $(TEST_OBJECTS) *.dep */*.dep */*/*.dep *.gch */*.gch */*/*.gch test_thrust test_cusp test_grid a.out

clear:
	rm *~ */*~ */*/*~

etags:
	etags --langmap=c++:+.cu.inl -R --links=no --exclude=state/flow_test --exclude=grid.2d --exclude=state/phi3 --exclude='cusp*' --exclude='thrust*' --exclude=state/phi3 --exclude='fish*' .


.PHONY: all clean clear etags cp
