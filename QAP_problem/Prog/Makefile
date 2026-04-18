
lsmcqap: qap.o timer.o stat.o
	gcc -O3 qap.o timer.o stat.o -o lsmcqap -lm -Wall -ansi

qap.o: qap.c
	gcc -O3 -c qap.c -o qap.o  -Wall -ansi

timer.o: timer.c timer.h
	gcc -O3 -c timer.c -o timer.o  -Wall -ansi

stat.o: stat.c stat.h
	gcc -O3 -c stat.c -o stat.o  -Wall -ansi
