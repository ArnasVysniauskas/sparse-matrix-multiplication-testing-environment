all:
	gcc -o app_1 -iquote . ./application/main.c -Wall -lgslcblas
	./app_1
