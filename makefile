CC = g++
CFLAGS = -O3 --std=c++17 -stdlib=libc++ -I /opt/homebrew/Cellar/boost/1.82.0_1/include #-Wall

integrals.exe : random.o integrals.o                        
	$(CC) random.o integrals.o  -o integrals.exe
	
integrals.o : integrals.cpp 
	$(CC) -c integrals.cpp -o integrals.o $(CFLAGS)

random.o : random.cpp 
	$(CC) -c random.cpp -o random.o $(CFLAGS)

clean :
	rm *.o integrals.exe data.out

