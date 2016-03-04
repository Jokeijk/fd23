OBJ= ../param/param.o ../getpar/getpar.o waveform.o source.o
OBJ_line= ../param/param.o ../getpar/getpar.o waveform.o source_line.o
OBJ_line_complex= ../param/param.o ../getpar/getpar.o waveform.o source_line_complex.o
CUOBJ=cu_source.cu
CC=nvcc
PTHREAD=-lpthread
CFLAGS=-g -G
CFLAGS=-O3 -DUSE_RESTRICT
%.o:%.cpp 
	${CC} ${CFLAGS} -c -o $@ $< 
%.o:%.cu
	${CC} ${CFLAGS} -arch sm_20 -c -o $@ $< 

ntde_line_complex:${OBJ_line_complex} ${CUOBJ}
	${CC} --ptxas-options=-v -m64 -arch sm_20 ${CFLAGS} ${PTHREAD} -o $@ $^ 

ntdenew:${OBJ} ${CUOBJ}
	${CC} --ptxas-options=-v -m64 -arch sm_20 ${CFLAGS} ${PTHREAD} -o $@ $^ 

ntde_line:${OBJ_line} ${CUOBJ}
	${CC} --ptxas-options=-v -m64 -arch sm_20 ${CFLAGS} ${PTHREAD} -o $@ $^ 

clean:
	\rm -f *.o
