trimesh: trimesh.o triangle.o
	g++ -o trimesh trimesh.o triangle.o

trimesh.o: trimesh.cpp constants.h triangle.h 
	g++ -c trimesh.cpp

triangle.o: triangle.cpp constants.h triangle.h
	g++ -c triangle.cpp

clean:
	rm *.o trimesh
