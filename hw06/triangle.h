#include <iostream>                                                             
#include <stdio.h>                                                              
                                                                                
#ifndef _CONST
#define _CONST
#include "constants.h"
#endif
                                                                                
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
