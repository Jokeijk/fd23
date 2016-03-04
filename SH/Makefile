OBJ=cpml/cpml.o fd/fd.o getpar/getpar.o field/field.o record/record.o record/isis.o material/material.o material/mktable.o model/model.o model/util.o param/param.o point_source/waveform.o record/needsnap.o

CUFILE=cpml/cu_cpml.cu field/cu_field.cu fd/cu_fd.cu material/cu_material.cu point_source/cu_point_source.cu record/cu_record.cu slide/cu_slide.cu thread/cu_thread.cu


all:sh2d sh2d_topabsorb
CC=nvcc

PTHREAD=-lpthread

#CFLAGS=-g -G
CFLAGS=-O3 -DUSE_RESTRICT

%.o:%.cpp 
	${CC} ${CFLAGS} -c -o $@ $< 
sh2d:main.o ${OBJ} ${CUFILE}
	${CC} -Dtopabsorb=false --ptxas-options=-v -m64 -arch sm_20 ${CFLAGS} ${PTHREAD} -o $@ $^ 
sh2d_topabsorb:main.o ${OBJ} ${CUFILE}
	${CC} -Dtopabsorb=true --ptxas-options=-v -m64 -arch sm_20 ${CFLAGS} ${PTHREAD} -o $@ $^ 


clean:
	\rm -f *.o */*.o
