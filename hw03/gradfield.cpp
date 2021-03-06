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

extern "C" {
  void sgesv_(const int* nn1,const int* nrhs,float* kp,const int* nn2,int* piv, float* BV,const int* nn3, int* info);
}

#include "triangle.h"

using namespace std;

Triangle* Tris;
int** Elems; //elems[ele#][node ABC]
int Elems_len;
float** Pts;
int Pts_len;

void calc_laplace(float **, float *);
void calc_grad(float *,float*,float*);
int main () {


  //======================
  //Read in the data into appropriate structures.
  //======================

  //Read in the Pts file into the Pts array
  FILE *Pts_file;
  if((Pts_file = fopen ("pts","r"))==NULL){
    cerr << "Unable to open Pts file." << endl;
  } else {
    fscanf(Pts_file,"%i",&Pts_len);
    Pts=(float**)malloc(sizeof(float*)*Pts_len);
    for(int i=0;i<Pts_len;i++){
      Pts[i]=(float*)malloc(sizeof(float)*3);
      fscanf(Pts_file,"%f %f %f",&Pts[i][X],&Pts[i][Y],&Pts[i][Z]);
    }
  }

  //Read in the ele file into the Elems array
  FILE *ele_file;
  if((ele_file = fopen ("ele","r"))==NULL){
    cerr << "Unable to open ele file." << endl;
  } else {
    fscanf(ele_file,"%i",&Elems_len);
    Elems=(int**)malloc(sizeof(int*)*Elems_len);
    for(int i=0;i<Elems_len;i++){
      Elems[i]=(int*)malloc(sizeof(int)*3);
      fscanf(ele_file,"%i %i %i",&Elems[i][A],&Elems[i][B],&Elems[i][C]);
      Elems[i][A]=Elems[i][A]-1; //Shift by 1 to properly index into Pts array.
      Elems[i][B]=Elems[i][B]-1;
      Elems[i][C]=Elems[i][C]-1;
    }
  }

  //Create Triangles
  Tris=(Triangle*)malloc(sizeof(Triangle)*Elems_len);
  for(int i=0;i<Elems_len;i++){                                                 
    Tris[i]=Triangle(Pts[Elems[i][A]],Pts[Elems[i][B]],Pts[Elems[i][C]]);  
  }

  //=====================
  //Calculate the Global Matrix, K and Initialize the Bound Values, BV
  //=====================  
  float** K;
  K=(float**)malloc(sizeof(float*)*Pts_len);
  for(int i=0;i<Pts_len;i++){
    K[i]=(float*)malloc(sizeof(float)*Pts_len);
  }
  float* BV=new float[Pts_len];
  calc_laplace(K,BV);

  //====================
  //Solve the Laplacian by inverting the matrix and multiplying.
  //====================
  //Both RHS and LHS are created, solve for BV!
  float* phi=new float[Pts_len];
  int info;
  int* piv=new int[Pts_len];
  int nrhs=1;
  const int nn=Pts_len;
  float *kp=new float[Pts_len*Pts_len];
  for(int i=0;i<Pts_len;i++){
    for(int j=0;j<Pts_len;j++){
      kp[i*Pts_len+j]=K[j][i];
    }
  }
  sgesv_(&nn,&nrhs,kp,&nn,piv,BV,&nn,&info);
  
  for(int i=0;i<10;i++){
    cout <<BV[i]<<endl;
  }
    
  //===================
  //Calculate the Gradient Field
  //===================
  //Now that we have the scalar field, BV, solve for the vector field V.
  float* vx=new float[Pts_len];
  float* vy=new float[Pts_len];
  calc_grad(vx,vy,BV);
  

  //==================
  //Write Scalar field and Vector Field to files.
  //==================
  FILE *scl_file;
  FILE *vec_file;
  scl_file = fopen("scl_fld.dat","w");
  vec_file = fopen("vec_fld.dat","w");
  for(int i=0;i<Pts_len;i++){
    fprintf(scl_file,"%f\n",BV[i]);
    fprintf(vec_file,"%f %f\n",vx[i],vy[i]);
  }

}

void calc_grad(float* vx,float* vy, float* phi){
  //Calculate the gradient of the scalar field BV

  float* Rx=new float[Pts_len];  
  float* Ry=new float[Pts_len];  
  float* M=new float[Pts_len];  
  float smx, smy;
  float yCA,yBA,xCA,xBA;
  //The computation begins here.
  for(int i=0;i<Elems_len;i++){
    yCA=Tris[i].gety(C)-Tris[i].gety(A); //differences used to calc grad(N)
    yBA=Tris[i].gety(B)-Tris[i].gety(A);
    xCA=Tris[i].getx(C)-Tris[i].getx(A);
    xBA=Tris[i].getx(B)-Tris[i].getx(A);
    smx=(yBA-yCA)*(phi[Elems[i][A]]/6.0)+(yCA)*(phi[Elems[i][B]]/6.0)+(-yBA)*(phi[Elems[i][C]]/6.0);
    smy=(xCA-xBA)*(phi[Elems[i][A]]/6.0)+(-xCA)*(phi[Elems[i][B]]/6.0)+(xBA)*(phi[Elems[i][C]]/6.0);
    Rx[Elems[i][A]]+=smx;
    Rx[Elems[i][B]]+=smx;
    Rx[Elems[i][C]]+=smx;
    Ry[Elems[i][A]]+=smy;
    Ry[Elems[i][B]]+=smy;
    Ry[Elems[i][C]]+=smy;
    M[Elems[i][A]]+=Tris[i].calcArea()/3.0;
    M[Elems[i][B]]+=Tris[i].calcArea()/3.0;
    M[Elems[i][C]]+=Tris[i].calcArea()/3.0;
  }
  for(int i=0;i<Pts_len;i++){
    vx[i]=Rx[i]/M[i];
    vy[i]=Ry[i]/M[i];
  }
  
  delete [] M;
  delete [] Rx;
  delete [] Ry;
}

void calc_laplace(float **K,float *BV){

  float yCA,yBA,xCA,xBA;
  float NXA,NXB,NXC,NYA,NYB,NYC;
  float area;
  for(int i=0;i<Elems_len;i++){
    area=Tris[i].calcArea();
    yCA=Tris[i].gety(C)-Tris[i].gety(A); //differences used to calc grad(N)
    yBA=Tris[i].gety(B)-Tris[i].gety(A);
    xCA=Tris[i].getx(C)-Tris[i].getx(A);
    xBA=Tris[i].getx(B)-Tris[i].getx(A);
    NXA=(yBA-yCA)/2.0;
    NXB=(yCA)/2.0;
    NXC=(-yBA)/2.0;
    NYA=(xCA-xBA)/2.0;
    NYB=(-xCA)/2.0;
    NYC=(xBA)/2.0;
    K[Elems[i][A]][Elems[i][A]]+=(NXA*NXA+NYA*NYA)/area;
    K[Elems[i][B]][Elems[i][B]]+=(NXB*NXB+NYB*NYB)/area;
    K[Elems[i][C]][Elems[i][C]]+=(NXC*NXC+NYC*NYC)/area;
    K[Elems[i][A]][Elems[i][B]]+=(NXA*NXB+NYA*NYB)/area;
    K[Elems[i][A]][Elems[i][C]]+=(NXA*NXC+NYA*NYC)/area;
    K[Elems[i][B]][Elems[i][C]]+=(NXB*NXC+NYB*NYC)/area;
    K[Elems[i][B]][Elems[i][A]]+=(NXA*NXB+NYA*NYB)/area;
    K[Elems[i][C]][Elems[i][A]]+=(NXA*NXC+NYA*NYC)/area;
    K[Elems[i][C]][Elems[i][B]]+=(NXB*NXC+NYB*NYC)/area;
  }
  //Read in the bpin file and create the BV array
  FILE *bpin_file;
  int bp_in;
  int bpin_len;
  FILE *bpout_file;
  int bp_out;
  int bpout_len;
  bpout_file = fopen ("bp.out","r");
  if((bpin_file = fopen ("bp.in","r"))==NULL){
    cerr << "Unable to open bp.in file." << endl;
  } else {    
    fscanf(bpin_file,"%i",&bpin_len);    
    fscanf(bpout_file,"%i",&bpout_len);    
    for(int i=0;i<bpin_len;i++){
      fscanf(bpin_file,"%i",&bp_in);
      fscanf(bpout_file,"%i",&bp_out);
      bp_in--;
      bp_out--;
      for(int j=0;j<Pts_len;j++){
	K[bp_in][j]=0;  //might have this ordering wrong!
	K[bp_out][j]=0;  //might have this ordering wrong!
      }
      K[bp_in][bp_in]=1;
      K[bp_out][bp_out]=1;
      BV[bp_in]=1;
      BV[bp_out]=-1;
    }
  }
}
