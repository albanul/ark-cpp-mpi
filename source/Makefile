all:
	mpicc -c *.cpp
	mpicc -O3 *.cpp -lstdc++ -o test
	rm *.o
	sbatch -n 32 -p gputest impi ./test
	#watch -n 1 squeue -u albanul92_251973