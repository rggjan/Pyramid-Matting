run: matting
	./matting test.ppm trimap.pnm results

gdb: matting_debug
	gdb ./matting_debug

valgrind: matting_debug
	valgrind ./matting_debug test.ppm trimap.pnm results

matting_debug: matting.cpp
	g++ -Wall -Wextra -O0 -g matting.cpp -o matting_debug

matting: matting.cpp
	g++ -Wall -Wextra -O2 matting.cpp -o matting -fopenmp
