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

//Global Data
Triangle* Tris;
int** Elems; //elems[ele#][node ABC]
float** Pts; //Unstructured Grid
float** Spts; //Structured Grid
int Elems_len;
int Pts_len;

//Functions
void initialize();
void calc_grad(float *,float*,float*); //Assignment 2
void calc_laplace(float**,float *); //Assignment 3
void interpolate(float*,float*,float*,float*); //Assignment 4
float fmax(float,float,float); //Assignment 4 helper
void hyperbolic(float**,float,int,float*,float*); //Assignment 5

int main () {
  
  //======================
  //Read computation in the data into appropriate structures.
  //======================
  initialize();
  
  //=====================
  //Calculate the Global Matrix, K and Initialize the Bound Values, BV
  //=====================  
  float** K;
  K=(float**)malloc(sizeof(float*)*Pts_len);
  for(int i=0;i<Pts_len;i++){
    K[i]=(float*)malloc(sizeof(float)*Pts_len);
  }
  float* sfld=new float[Pts_len];
  
  calc_laplace(K,sfld);

  //====================
  //Solve the Laplacian by inverting the matrix and multiplying.
  //====================
  //Both RHS and LHS are created, solve for sfld!
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
  sgesv_(&nn,&nrhs,kp,&nn,piv,sfld,&nn,&info);
      
  //===================
  //Calculate the Gradient Field
  //===================
  //Now that we have the scalar field, sfld, solve for the vector field V.
  float* vx=new float[Pts_len];
  float* vy=new float[Pts_len];
  calc_grad(vx,vy,sfld);
  
  //==================
  //Interpolate the vector field on to a structured grid.
  //==================
  float* ivx=new float[Spts_len2];
  float* ivy=new float[Spts_len2];
  interpolate(vx,vy,ivx,ivy);

  //==================
  //Solve the Hyperbolic function on the interpolated field
  //==================
  float** field=new float*[Spts_len1];
  float** nfield=new float*[Spts_len1];
  for(int i=0;i<Spts_len1;i++){field[i]=new float[Spts_len1];}
  //Initialize the field values: columns 10-20 should be 1
  for(int i=0;i<Spts_len1;i++){
    for(int j=0;j<Spts_len1;j++){
      if(j<20 && j>9)
	field[i][j]=1;
      else
	field[i][j]=0;
    }
  }
  //Initialize the delt and timesteps values.
  float delx=Spts[1][X]-Spts[0][X];
  float a=0; //Set a to the largest velocity value
  for(int i=0;i<Spts_len2;i++){
    if(fabs(ivx[i])>a)
      a=fabs(ivx[i]);
    if(fabs(ivy[i])>a)
      a=fabs(ivy[i]);
  }
  float delt=delx/a;
  int timesteps=8000;
  hyperbolic(field,delt,timesteps,ivx,ivy);
  
}

void initialize(){
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

  //Initialize Structured Points
  
  Spts=new float*[Spts_len2];
  for(int i=0;i<Spts_len2;i++){
    Spts[i]=new float[2];
    Spts[i][X]=-start+1.0/Spts_len1+(i%Spts_len1)*start*2/Spts_len1;
  }
  for(int i=0;i<Spts_len1;i++){
    for(int j=0;j<Spts_len1;j++){
      Spts[i*Spts_len1+j][Y]=-start+1.0/Spts_len1+(i%Spts_len1)*start*2/Spts_len1;
    }
  }
}

void calc_grad(float* vx,float* vy, float* sfld){
  //Calculate the gradient of the scalar field sfld

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
    smx=(yBA-yCA)*(sfld[Elems[i][A]]/6.0)+(yCA)*(sfld[Elems[i][B]]/6.0)+(-yBA)*(sfld[Elems[i][C]]/6.0);
    smy=(xCA-xBA)*(sfld[Elems[i][A]]/6.0)+(-xCA)*(sfld[Elems[i][B]]/6.0)+(xBA)*(sfld[Elems[i][C]]/6.0);
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
      BV[bp_in]=-1;
      BV[bp_out]=1;
    }
  }

  fclose(bpin_file);
  fclose(bpout_file);
}

void interpolate(float* vx,float* vy,float* ivx,float* ivy){
  vector<vector <int> > Ptsnelems(Pts_len,vector<int>()); // point -> neighboring points
  vector<int>::iterator iter1;
  vector<int>::iterator iter2;
  for(int i=0;i<Elems_len;i++){
    Ptsnelems[Elems[i][A]].push_back(i);
    Ptsnelems[Elems[i][B]].push_back(i);
    Ptsnelems[Elems[i][C]].push_back(i);
  }

  int e=0;
  float xpa,ypa,xpt,ypt;
  float g1x,g2x,g1y,g2y;
  float g11,g12,g22;
  float det;
  float eta,zeta,xi;
  int found,p1,p2;
  float dir;

  for(int i=0;i<Spts_len2;i++){
    xpt=Spts[i][X];
    ypt=Spts[i][Y];
    while(true){
      //Create the shape functions eta, zeta and xi
      g1x=Pts[Elems[e][B]][X]-Pts[Elems[e][A]][X];
      g2x=Pts[Elems[e][C]][X]-Pts[Elems[e][A]][X];
      xpa=xpt-Pts[Elems[e][A]][X];
      g1y=Pts[Elems[e][B]][Y]-Pts[Elems[e][A]][Y];
      g2y=Pts[Elems[e][C]][Y]-Pts[Elems[e][A]][Y];
      ypa=ypt-Pts[Elems[e][A]][Y];
      g11=g1x*g1x+g1y*g1y;
      g12=g1x*g2x+g1y*g2y;
      g22=g2x*g2x+g2y*g2y;
      det=g11*g22-g12*g12;
      eta=(g22*(xpa*g1x+ypa*g1y)-g12*(xpa*g2x+ypa*g2y))/det;
      zeta=(g11*(xpa*g2x+ypa*g2y)-g12*(xpa*g1x+ypa*g1y))/det;
      xi=1-eta-zeta;
      //Check if containing Element is found
      if(eta>0 && 1-eta>0 && zeta>0 && 1-zeta>0 && xi>0 && 1-xi>0){
	//cout << e << endl;
	ivx[i]=xi*vx[Elems[e][A]]+eta*vx[Elems[e][B]]+zeta*vx[Elems[e][C]];
	ivy[i]=xi*vy[Elems[e][A]]+eta*vy[Elems[e][B]]+zeta*vy[Elems[e][C]];
	break;
      }
      dir=fmax(eta,zeta,xi);
      if(dir==xi){
	p1=B;
	p2=C;
      }
      if(dir==eta){
	p1=A;
	p2=C;
      }
      if(dir==zeta){
	p1=A;
	p2=B;
      }
      if(dir==0){
	cout << "err";
	fflush(stdout);
      }
      found=0;
      for(iter1=Ptsnelems[Elems[e][p1]].begin();iter1!=Ptsnelems[Elems[e][p1]].end();iter1++){
	for(iter2=Ptsnelems[Elems[e][p2]].begin();iter2!=Ptsnelems[Elems[e][p2]].end();iter2++){
	  //cout << e << ":" << *iter1 << " " << *iter2 << endl;
	  if(*iter1==*iter2 && *iter1!=e){
	    found=1;
	    e=*iter1;
	    break;
	  }
	}
	if(found) break;
      }
      if(!found){
	//No element contains the point (attempted to step over boundary)
	//OR you are trying to cross the center circle, so try brute force before giving up.
	for(e=0;e<Elems_len;e++){
	  //Create the shape functions eta, zeta and xi
	  g1x=Pts[Elems[e][B]][X]-Pts[Elems[e][A]][X];
	  g2x=Pts[Elems[e][C]][X]-Pts[Elems[e][A]][X];
	  xpa=xpt-Pts[Elems[e][A]][X];
	  g1y=Pts[Elems[e][B]][Y]-Pts[Elems[e][A]][Y];
	  g2y=Pts[Elems[e][C]][Y]-Pts[Elems[e][A]][Y];
	  ypa=ypt-Pts[Elems[e][A]][Y];
	  g11=g1x*g1x+g1y*g1y;
	  g12=g1x*g2x+g1y*g2y;
	  g22=g2x*g2x+g2y*g2y;
	  det=g11*g22-g12*g12;
	  eta=(g22*(xpa*g1x+ypa*g1y)-g12*(xpa*g2x+ypa*g2y))/det;
	  zeta=(g11*(xpa*g2x+ypa*g2y)-g12*(xpa*g1x+ypa*g1y))/det;
	  xi=1-eta-zeta;
	  if(eta>0 && 1-eta>0 && zeta>0 && 1-zeta>0 && xi>0 && 1-xi>0){
	    //cout << e << endl;
	    ivx[i]=xi*vx[Elems[e][A]]+eta*vx[Elems[e][B]]+zeta*vx[Elems[e][C]];
	    ivy[i]=xi*vy[Elems[e][A]]+eta*vy[Elems[e][B]]+zeta*vy[Elems[e][C]];
	    found=1;
	    break;
	  }
	}
	if(!found){
	  e=0;
	  ivx[i]=0;
	  ivy[i]=0;
	}
    	break;
      }
    }
  }
}

float fmax(float a, float b, float c){
  if(a==fabs(a))
    a=0;
  else
    a=fabs(a);
  if(b==fabs(b))
    b=0;
  else
    b=fabs(b);
  if(c==fabs(c))
    c=0;
  else
    c=fabs(c);
  return (float)max(max(a,b),max(a,c))*-1.0;
}

void hyperbolic(float** fld, float delt, int steps,float* vx, float*vy){
  //Start by writing out the initial field to a file.
  FILE *sfld_file;
  sfld_file = fopen("start_fld.dat","w");
  for(int i=0;i<Spts_len1;i++){
    for(int j=0;j<Spts_len1;j++){
      fprintf(sfld_file,"%f ",fld[i][j]);
    }
    fprintf(sfld_file,"\n");
  }

  FILE *active_fld_file;
  active_fld_file = fopen("act_fld.dat","w");
  float** nfld=new float*[Spts_len1];
  for(int i=0;i<Spts_len1;i++){nfld[i]=new float[Spts_len1];}
  int count=0;
  int index,xi1,xi2,yj1,yj2;
  int vxmask,vymask;
  for(int t=0;t<steps;t++){
    for(int i=0;i<Spts_len1;i++){
      for(int j=0;j<Spts_len1;j++){
	index=i*Spts_len1+j;
	//Mask components based on bounds
	(i==0 || i==Spts_len1-1) ? vymask=0 : vymask=1;
	(j==0 || j==Spts_len1-1) ? vxmask=0 : vxmask=1;
	//Upwind differential, differentiate in direction of velocity
	//Multiply by mask to prevent indexing out of bounds
	nfld[i][j]=fld[i][j];
	if(vx[index]<0 && vxmask){
	  nfld[i][j]-=vx[index]*delt*(fld[i][j+1]-fld[i][j]);
	} 
	if(vx[index]>0 && vxmask){
	  nfld[i][j]-=vx[index]*delt*(fld[i][j]-fld[i][j-1]);
	}
	if(vy[index]<0 && vymask){
	  nfld[i][j]-=vy[index]*delt*(fld[i+1][j]-fld[i][j]);
	} 
	if(vy[index]>0 && vymask){
	  nfld[i][j]-=vy[index]*delt*(fld[i][j]-fld[i-1][j]);
	}

      }
    }
    for(int i=0;i<Spts_len1;i++){
      for(int j=0;j<Spts_len1;j++){
	//Propagate the field change.
	fld[i][j]=nfld[i][j];
      }
    }
    if(t%20==0){
      //Write out data. 
      for(int i=0;i<Spts_len1;i++){
	for(int j=0;j<Spts_len1;j++){
	  fprintf(active_fld_file,"%f ",nfld[i][j]);
	}
	fprintf(active_fld_file,"\n");
      }
    }
  }
  
  FILE *interp_xpts_file;
  FILE *interp_ypts_file;
  FILE *vec_file;
  interp_xpts_file = fopen("ixpts.dat","w");
  interp_ypts_file = fopen("iypts.dat","w");
  for(int i=0;i<Spts_len2;i++){
    fprintf(interp_xpts_file,"%f ",Spts[i][X]);
    fprintf(interp_ypts_file,"%f ",Spts[i][Y]);
    if(i!=0 && i%Spts_len1==99){
      fprintf(interp_xpts_file,"\n");
      fprintf(interp_ypts_file,"\n");
    }
  }
}

