all:
	gcc raytrace.c -o raytrace -lm

clean:
	rm -rf raytrace *~
