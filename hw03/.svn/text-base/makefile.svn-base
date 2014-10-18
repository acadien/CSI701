#INCLUDES	= -I/Users/acadien/LAPACKPP/include
LIBS		=-lgfortran -llapack
gradfield: gradfield.o triangle.o
	g++ $(LIBS) -o gradfield gradfield.o triangle.o

gradfield.o: gradfield.cpp triangle.h constants.h
	g++ -c gradfield.cpp

triangle.o: triangle.cpp triangle.h constants.h
	g++ -c triangle.cpp

clean:
	rm *.o gradfield
