#ifndef _CONST
#define _CONST
#include "constants.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <math.h>
#include <vector>

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

  //Read in the fld file into the scl_fld array
  FILE *fld_file;
  float *scl_fld; //scl_fld[#fld]=scalar
  int ndim_fld;
  int n=0;
  if((fld_file = fopen ("fld","r"))==NULL){
    cerr << "Unable to open fld file." << endl;
  } else {
    fscanf(fld_file,"%i\n",&ndim_fld);
    if(ndim_fld != 1) {
      cerr << "Scalar field should have a ndim=1, actually ndim=" << ndim_fld << endl;
    }
    scl_fld=(float*)malloc(sizeof(float)*pts_len);
    //while(feof(fld_file)==0){
    for(int i=0;i<pts_len;i++){
      fscanf(fld_file,"%f",&scl_fld[i]);
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
      elems[i][A]=elems[i][A]-1; //Shift by 1 to properly index into pts array.
      elems[i][B]=elems[i][B]-1;
      elems[i][C]=elems[i][C]-1;
    }
  }

  //Create Triangles
  Triangle* tris=(Triangle*)malloc(sizeof(Triangle)*elems_len);
  for(int i=0;i<elems_len;i++){                                                 
    tris[i]=Triangle(pts[elems[i][A]],pts[elems[i][B]],pts[elems[i][C]]);  
  }

  //=====================
  //Calculate the Vector Field
  //=====================
  
  vector<float> Rx(pts_len);
  vector<float> Vx(pts_len);
  vector<float> Ry(pts_len);
  vector<float> Vy(pts_len);
  vector<float> M(pts_len);
  float yCA,yBA,xCA,xBA;
  float smx, smy;
  //The computation begins here.
  for(int i=0;i<elems_len;i++){
    yCA=tris[i].gety(C)-tris[i].gety(A); //differences used to calc grad(N)
    yBA=tris[i].gety(B)-tris[i].gety(A);
    xCA=tris[i].getx(C)-tris[i].getx(A);
    xBA=tris[i].getx(B)-tris[i].getx(A);
    smx=(yBA-yCA)*(scl_fld[elems[i][A]]/6.0)+(yCA)*(scl_fld[elems[i][B]]/6.0)+(-yBA)*(scl_fld[elems[i][C]]/6.0);
    smy=(xCA-xBA)*(scl_fld[elems[i][A]]/6.0)+(-xCA)*(scl_fld[elems[i][B]]/6.0)+(xBA)*(scl_fld[elems[i][C]]/6.0);
    Rx[elems[i][A]]+=smx;
    Rx[elems[i][B]]+=smx;
    Rx[elems[i][C]]+=smx;
    Ry[elems[i][A]]+=smy;
    Ry[elems[i][B]]+=smy;
    Ry[elems[i][C]]+=smy;
    M[elems[i][A]]+=tris[i].calcArea()/3.0;
    M[elems[i][B]]+=tris[i].calcArea()/3.0;
    M[elems[i][C]]+=tris[i].calcArea()/3.0;
  }
  for(int i=0;i<pts_len;i++){
    Vx[i]=Rx[i]/M[i];
    Vy[i]=Ry[i]/M[i];
  }
  printf("3\n");
  for(int i=0;i<pts_len;i++){
    printf("%f %f 0\n",Vx[i],Vy[i]);
  }

  //Create neighboring points and elements lists
  //Don't need this yet.
  /*
  vector<vector <int> > ptsnpts(pts_len,vector<int>()); // point -> neighboring points.
  vector<vector <int> > ptsnelems(pts_len,vector<int>()); //point[pt] -> neighboring elements.
  vector<int>::iterator iter;
  vector<int>::iterator iter2;
  
  for(int i=0;i<elems_len;i++){
    ptsnelems[elems[i][A]].push_back(i);
    ptsnelems[elems[i][B]].push_back(i);
    ptsnelems[elems[i][C]].push_back(i);
  }
  int exists=0;
  for(int i=0;i<pts_len;i++){
    //cout << i << ": ";
    for(iter=ptsnelems[i].begin();iter!=ptsnelems[i].end();iter++){
      for(int j=0;j<3;j++){
	if(elems[*iter][j]!=i) {
	  if(ptsnpts[i].size()==0){
	    ptsnpts[i].push_back(elems[*iter][j]);
	    //cout << elems[*iter][j] << " ";
	  } else {
	    exists=0;
	    for(iter2=ptsnpts[i].begin();iter2!=ptsnpts[i].end();iter2++){
	      if(elems[*iter][j]==*iter2) {
		exists++;
	      }
	    }
	    if(!exists){
	      ptsnpts[i].push_back(elems[*iter][j]);
	      //cout << elems[*iter][j] << " ";
	    }
	  } 
	}
      }
    }
    //cout << endl;
  }
  //ptsnpts[pt]=vector<int> of pts neighboring this point.
  */
  
}

