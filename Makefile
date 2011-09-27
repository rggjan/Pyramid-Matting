run: matting
	./matting

valgrind: matting_debug
	valgrind ./matting_debug

matting_debug: matting.cpp
	g++ -Wall -Wextra -g matting.cpp -o matting_debug

matting: matting.cpp
	g++ -Wall -Wextra -O3 matting.cpp -o matting
