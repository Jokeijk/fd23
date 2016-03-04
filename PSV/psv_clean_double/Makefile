OBJ=box/box.o cpml/cpml.o fd/fd.o getpar/getpar.o field/field.o record/record.o record/isis.o material/material.o material/mktable.o model/model.o model/util.o param/param.o source/waveform.o

CUOBJ=box/cu_box.cu cpml/cu_cpml.cu field/cu_field.cu fd/cu_fd.cu material/cu_material.cu point_source/cu_point_source.cu record/cu_record.cu slide/cu_slide.cu thread/cu_thread.cu

CC=nvcc
PTHREAD=-lpthread
CFLAGS=-g -G 
CFLAGS=-O3 -DUSE_RESTRICT
%.o:%.cpp 
	${CC} ${CFLAGS} -c -o $@ $< 

all:psv2d_db psv2d_db_topabsorb
psv2d_db:${OBJ} ${CUOBJ} main.cpp
	${CC} -Dtopabsorb=false --ptxas-options=-v -m64 -arch sm_20 ${CFLAGS} ${PTHREAD} -o $@ $^ 

psv2d_db_topabsorb:${OBJ} ${CUOBJ} main.cpp
	${CC} -Dtopabsorb=true --ptxas-options=-v -m64 -arch sm_20 ${CFLAGS} ${PTHREAD} -o $@ $^ 


CUOBJ_Q=box/cu_box.o cpml/cu_cpml.o field/cu_field.o fd/cu_fd_Q.o material/cu_material.o point_source/cu_point_source.o record/cu_record.o slide/cu_slide.o thread/cu_thread.o

psv2d_q:${OBJ} ${CUOBJ_Q} main.o 
	${CC} --ptxas-options=-v -m64 -arch sm_20 ${CFLAGS} ${PTHREAD} -o $@ $^ 


point_theory:${OBJ} ${CUOBJ} theory/point_theory.o source/cu_source.cu
	${CC} --ptxas-options=-v -m64 -arch sm_20 ${CFLAGS} ${PTHREAD} -o $@ $^ 

ntdenew:
	cd source && make
test_material:
	cd test_material && make
clean:
	\rm -f *.o */*.o
