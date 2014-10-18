#INCLUDES	= -I/Users/acadien/LAPACKPP/include
LIBS		=-lgfortran -llapack
particulate: particulate.o triangle.o
	g++ $(LIBS) -o particulate particulate.o triangle.o

particulate.o: particulate.cpp triangle.h constants.h
	g++ -c particulate.cpp

triangle.o: triangle.cpp triangle.h constants.h
	g++ -c triangle.cpp

clean:
	rm *.o particulate
