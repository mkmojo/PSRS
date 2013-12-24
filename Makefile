CFLAGS = -Wall

run: CFLAGS += -O3 
run: psrs

debug: CFLAGS += -DDEBUG -g
debug: psrs

psrs: ss.o psrs.o main.o 
	mpicc ss.o psrs.o main.o -o psrsTest

main.o: main.c
	mpicc -c main.c $(CFLAGS)

psrs.o: psrs.c psrs.h
	mpicc -c psrs.c $(CFLAGS)

ss.o: sequential_sort.c sequential_sort.h
	mpicc -c sequential_sort.c -o ss.o $(CFLAGS)

clean:
	rm -f a.out psrsTest *.o

new:
	make clean 1> /dev/null && make run 1> /dev/null

newdebug:
	make clean 1> /dev/null && make debug 1> /dev/null
