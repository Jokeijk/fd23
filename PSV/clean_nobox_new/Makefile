OBJ=cpml/cpml.o getpar/getpar.o field/field.o record/record.o record/isis.o material/material.o material/mktable.o model/model.o model/util.o param/param.o source/waveform.o

CUOBJ=cpml/cu_cpml.cu field/cu_field.cu fd/cu_fd.cu material/cu_material.cu point_source/cu_point_source.cu record/cu_record.cu slide/cu_slide.cu thread/cu_thread.cu

CC=/opt/cuda/bin/nvcc
PTHREAD=-lpthread
CFLAGS=-g -G
CFLAGS=-O3 -DUSE_RESTRICT
%.o:%.cpp 
	${CC} ${CFLAGS} -c -o $@ $< 
%.o:%.c 
	${CC} ${CFLAGS} -c -o $@ $< 

all:nbpsv2d nbpsv2d_topabsorb nbpsv2d_exp nbpsv2d_exp_topabsorb

nbpsv2d:${OBJ} ${CUOBJ} main.cpp
	${CC} -Dtopabsorb=false --ptxas-options=-v -m64 -arch sm_20 ${CFLAGS} ${PTHREAD} -o $@ $^ 

nbpsv2d_topabsorb:${OBJ} ${CUOBJ} main.cpp
	${CC} -Dtopabsorb=true --ptxas-options=-v -m64 -arch sm_20 ${CFLAGS} ${PTHREAD} -o $@ $^ 

nbpsv2d_exp:${OBJ} ${CUOBJ} main.cpp
	${CC} -Dtopabsorb=false -DEXPLOSION_SOURCE --ptxas-options=-v -m64 -arch sm_20 ${CFLAGS} ${PTHREAD} -o $@ $^ 

nbpsv2d_exp_topabsorb:${OBJ} ${CUOBJ} main.cpp
	${CC} -Dtopabsorb=true -DEXPLOSION_SOURCE --ptxas-options=-v -m64 -arch sm_20 ${CFLAGS} ${PTHREAD} -o $@ $^ 

zip:
	tar czf psv2d.tgz box fd getpar main.cpp material param cpml field model gpu.h Makefile point_source record slide thread source/waveform.* timer.h 
clean:
	\rm -f *.o */*.o
                                 
