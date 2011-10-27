run: matting
	./matting test.ppm trimap.pnm results --debug --gt GT02.pnm

batch: matting
	./batch.rb
	for input in test/input/*; do \
	  name=${input:%.ppm=%}; \
	  echo $$name; \
	  echo $$input; \
	done

gdb: matting_debug
	gdb ./matting_debug

valgrind: matting_debug
	valgrind ./matting_debug test.ppm trimap.pnm results

matting_debug: matting.cpp
	g++ -Wall -Wextra -O0 -g matting.cpp -o matting_debug

matting: matting.cpp
	g++ -Wall -Wextra -O2 matting.cpp -o matting -fopenmp
