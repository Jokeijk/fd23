OBJ= ../param/param.o ../getpar/getpar.o ../material/material.o ../model/model.o  ../model/util.o test_material.o
CUOBJ=
CC=nvcc
PTHREAD=-lpthread
CFLAGS=-g -G
CFLAGS=-O3 -DUSE_RESTRICT
%.o:%.cpp 
	${CC} ${CFLAGS} -c -o $@ $< 
%.o:%.cu
	${CC} ${CFLAGS} -c -o $@ $< 
test_material:${OBJ} ${CUOBJ}
	${CC} --ptxas-options=-v -m64 -arch sm_20 ${CFLAGS} ${PTHREAD} -o $@ $^ 
clean:
	\rm -f ${OBJ} ${CUOBJ}
