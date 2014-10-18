#include <iostream>
#include <stdio.h>

#include "constants.h"

class Triangle
{
  float* apnt;
  float* bpnt;
  float* cpnt;
public:
  Triangle(float* a,float* b,float* c);
  ~Triangle(){};
  float getx(int node);
  float gety(int node);
  float getz(int node);
  float calcArea();
};
