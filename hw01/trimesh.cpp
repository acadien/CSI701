#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <math.h>

#include "triangle.h"

using namespace std;

int main () {

  //======================
  //Read in the data into appropriate structures.
  //======================

  //Read in the pts file into the pts array
  FILE *pts_file;
  float** pts; //pts[pt#][dirn XYZ]
  int pts_len;
  if((pts_file = fopen ("pts","r"))==NULL){
    cerr << "Unable to open pts file." << endl;
  } else {
    fscanf(pts_file,"%i",&pts_len);
    pts=(float**)malloc(sizeof(float*)*pts_len);
    for(int i=0;i<pts_len;i++){
      pts[i]=(float*)malloc(sizeof(float)*3);
      fscanf(pts_file,"%f %f %f",&pts[i][X],&pts[i][Y],&pts[i][Z]);
    }
  }

  //Read in the ele file into the elems array
  FILE *ele_file;
  int** elems; //elems[ele#][node ABC]
  int elems_len;
  if((ele_file = fopen ("ele","r"))==NULL){
    cerr << "Unable to open ele file." << endl;
  } else {
    fscanf(ele_file,"%i",&elems_len);
    elems=(int**)malloc(sizeof(int*)*elems_len);
    for(int i=0;i<elems_len;i++){
      elems[i]=(int*)malloc(sizeof(int)*3);
      fscanf(ele_file,"%i %i %i",&elems[i][A],&elems[i][B],&elems[i][C]);
      elems[i][A]--; //Shift by 1 to properly index into pts array.
      elems[i][B]--;
      elems[i][C]--;
    }
  }

  //======================
  //Find the Bounding Box.
  //======================
  float xmin,xmax;
  float ymin,ymax;
  xmin=pts[0][X];
  xmax=pts[0][X];
  ymin=pts[0][Y];
  ymax=pts[0][Y];
  for(int i=1;i<pts_len;i++){
    if(pts[i][X]>xmax){
      xmax=pts[i][X];
    }
    if(pts[i][X]<xmin){
      xmin=pts[i][X];
    }
    if(pts[i][Y]>ymax){
      ymax=pts[i][Y];
    }
    if(pts[i][Y]<ymin){
      ymin=pts[i][Y];
    }
  }
  printf("\nBounding box is formed by:\nxmin = %f\nxmax = %f\nymin = %f\nymax = %f\n\n",xmin,xmax,ymin,ymax);

  //=====================
  //Create Triangle Objects and find the area.
  //=====================
  Triangle* tris=(Triangle*)malloc(sizeof(Triangle)*elems_len);
  float total_area=0;
  for(int i=0;i<elems_len;i++){
    tris[0]=Triangle(pts[elems[i][A]],pts[elems[i][B]],pts[elems[i][C]]);
    total_area+=tris[0].calcArea();
  }
  printf("Total Area of Triangles = %f\n\n",total_area);
}
