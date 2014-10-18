#INCLUDES	= -I/Users/acadien/LAPACKPP/include
LIBS		=-lgfortran -llapack
hyperbolic: hyperbolic.o triangle.o
	g++ $(LIBS) -o hyperbolic hyperbolic.o triangle.o

hyperbolic.o: hyperbolic.cpp triangle.h constants.h
	g++ -c hyperbolic.cpp

triangle.o: triangle.cpp triangle.h constants.h
	g++ -c triangle.cpp

clean:
	rm *.o hyperbolic
