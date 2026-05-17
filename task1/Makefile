CC = mpicc
LIBS = -lm
SOURCES = main.c input.c output.c reflect_reverse_MPI.c
OBJS = main.o input.o output.o reflect_reverse_MPI.o
EXE = app

all : $(EXE)

%.o : %.c
	$(CC) -c $< $(LIBS)

$(EXE) : $(OBJS)
	$(CC) -o $@ $^ $(LIBS)

ifneq (clean, $(MAKECMDGOALS))
    -include deps.mk
endif

deps.mk: $(SOURCES)
	$(CC) -MM $^ > $@

clean:
	-rm $(OBJS)
	-rm $(EXE)
	-rm deps.mk


start:
	echo "mpirun -np p ./app n r s"
