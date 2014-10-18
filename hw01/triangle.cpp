#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include "triangle.h"

using namespace std;

Triangle::Triangle(float* a, float* b, float* c){
  apnt=(float*)malloc(3*sizeof(float));
  bpnt=(float*)malloc(3*sizeof(float));
  cpnt=(float*)malloc(3*sizeof(float));
  apnt=a;
  bpnt=b;
  cpnt=c;
}

float Triangle::getx(int node){
  switch (node) {
  case A:
    return apnt[X];
  case B:
    return bpnt[X];
  case C:
    return cpnt[X];
  }
  cerr << "Invalid switch in Triangle::getx.\n";
  return -1;
}

float Triangle::gety(int node){
  switch (node) {
  case A:
    return apnt[Y];
  case B:
    return bpnt[Y];
  case C:
    return cpnt[Y];
  }
  cerr << "Invalid switch in Triangle::gety.\n";
  return -1;
}

float Triangle::getz(int node){
  switch (node) {
  case A:
    return apnt[Z];
  case B:
    return bpnt[Z];
  case C:
    return cpnt[Z];
  }
  cerr << "Invalid switch in Triangle::getz.\n";
  return -1;
}

float Triangle::calcArea(){
  //Assumes a 2D triangle. (Ignores the Z direction)
  float x_ba=bpnt[X]-apnt[X];
  float x_ca=cpnt[X]-apnt[X];
  float y_ba=bpnt[Y]-apnt[Y];
  float y_ca=cpnt[Y]-apnt[Y];
  
  return fabs(0.5*((x_ba*y_ca)-(y_ba*x_ca)));
}
