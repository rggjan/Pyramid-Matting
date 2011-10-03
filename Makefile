run: matting
	./matting

gdb: matting_debug
	gdb ./matting_debug

valgrind: matting_debug
	valgrind ./matting_debug

matting_debug: matting.cpp
	g++ -Wall -Wextra -O0 -g matting.cpp -o matting_debug

matting: matting.cpp
	g++ -Wall -Wextra -O2 matting.cpp -o matting
