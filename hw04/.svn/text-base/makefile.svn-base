#INCLUDES	= -I/Users/acadien/LAPACKPP/include
LIBS		=-lgfortran -llapack
interpolate: interpolate.o triangle.o
	g++ $(LIBS) -o interpolate interpolate.o triangle.o

interpolate.o: interpolate.cpp triangle.h constants.h
	g++ -c interpolate.cpp

triangle.o: triangle.cpp triangle.h constants.h
	g++ -c triangle.cpp

clean:
	rm *.o interpolate
