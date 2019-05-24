
app : main.o globals.o
	g++ ./main.o ./globals.o

main.o : ./main.cpp
	g++ -O3 -m64 -c -o "main.o" main.cpp

globals.o : ./globals.cpp
	g++ -O3 -m64 -c -o "globals.o" globals.cpp

clean :
	-rm ./a.out ./main.o ./globals.o
