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
#include <ctime>

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
float* vx; //interpolated velocity field
float* vy; //interpolated velocity field
int Elems_len;
int Pts_len;
vector< vector <int> > Ptsnelems; //Ptsnelems[Pt][Elems] the Elems that neighbor each Pt

//Functions
void initialize();
void calc_grad(float *,float*,float*); //Assignment 2
void calc_laplace(float**,float *); //Assignment 3
void interpolate(float,float,float*,float*,int*); //Adjusted for Assign6
float fmax(float,float,float); 
void particulate(float**); //Assignment 6

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
  vx=new float[Pts_len];
  vy=new float[Pts_len];
  calc_grad(vx,vy,sfld);
  
  //==================
  //Interpolate the vector field on to a structured grid.
  //==================
  float* ivx=new float[Spts_len2];
  float* ivy=new float[Spts_len2];
  float xpt, ypt;
  int elem=0;
  for(int i=0;i<Spts_len2;i++){
    xpt=Spts[i][X];
    ypt=Spts[i][Y];
    interpolate(xpt,ypt,&(ivx[i]),&(ivy[i]),&elem);
  }

  //==================
  //Initialize a set of particles and propegate em!
  //==================
  float** particles=new float*[Num_Ptcls];
  //SEED!!!!
  srand(time(0));
  for(int i=0;i<Num_Ptcls;i++){
    particles[i]=new float[3];
    particles[i][X]=(float)(drand48()*0.2-0.8); //Random number [-1,1)
    particles[i][Y]=(float)(drand48()*2-0.9999); //Random number [-1,1)
    particles[i][T]=0;//Timestamp for the particle
  }
  particulate(particles);
  
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

  // Initialize Pts -> Elem array, useful structure for interpolation

  Ptsnelems.resize(Pts_len,vector<int>());
  vector<int>::iterator iter1;
  vector<int>::iterator iter2;
  for(int i=0;i<Elems_len;i++){
    Ptsnelems[Elems[i][A]].push_back(i);
    Ptsnelems[Elems[i][B]].push_back(i);
    Ptsnelems[Elems[i][C]].push_back(i);
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

void interpolate(float xpt,float ypt,float* ivx,float* ivy,int* elm){
  //Given the point (xpt,ypt) interpolate the velocity field to get the
  //field values which are set to (ivx,ivy)
  //The int* elm is the seed element. The nearest neighbor search starts at 
  //element elm and once the parent element is found elm is set to that.

  int e;
  if(*elm<0){ //previously no parent element was found
    e=0;
  } else {
    e=*elm;
  }
  
  float xpa,ypa;
  float g1x,g2x,g1y,g2y;
  float g11,g12,g22;
  float det;
  float eta,zeta,xi;
  int found,p1,p2;
  float dir;
  vector<int>::iterator iter1;
  vector<int>::iterator iter2;
 
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
      *ivx=xi*vx[Elems[e][A]]+eta*vx[Elems[e][B]]+zeta*vx[Elems[e][C]];
      *ivy=xi*vy[Elems[e][A]]+eta*vy[Elems[e][B]]+zeta*vy[Elems[e][C]];
      *elm=e;
      return;
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
      cout << "Error interpolating point, max shape function=0, location:\n";
      cout << xpt << " " << ypt << endl;
      fflush(stdout);
    }
    found=0;
    for(iter1=Ptsnelems[Elems[e][p1]].begin();iter1!=Ptsnelems[Elems[e][p1]].end();iter1++){
      for(iter2=Ptsnelems[Elems[e][p2]].begin();iter2!=Ptsnelems[Elems[e][p2]].end();iter2++){
	if(*iter1==*iter2 && *iter1!=e){
	  found=1;
	  e=*iter1;
	  break;
	}
      }
      if(found) break; //BREAKS FOR LOOP, CONTINUES WHILE LOOP
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
	  *ivx=xi*vx[Elems[e][A]]+eta*vx[Elems[e][B]]+zeta*vx[Elems[e][C]];
	  *ivy=xi*vy[Elems[e][A]]+eta*vy[Elems[e][B]]+zeta*vy[Elems[e][C]];
	  *elm=e;
	  return;
	}
      }
      //If you get here, then no element contains the point.
      *elm=-1;
      *ivx=0;
      *ivy=0;
      return;
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

void particulate(float** particles){
  int Active_Ptcls=Num_Ptcls;
  float GT=0; //Global current time.
  float p_vx,p_vy,delt;//Individual particle's velocity and timestep
  float k1x,k2x,k3x,k4x;//RungeKutta Coefficients
  float k1y,k2y,k3y,k4y;//RungeKutta Coefficients
  int elem=0,ready=Num_Ptcls;//ready: number of particles ready to be sampled
  int* readies = new int[Num_Ptcls];
  float *ptcls_vx=new float[Num_Ptcls];
  float *ptcls_vy=new float[Num_Ptcls];

  FILE *fptcls_x, *fptcls_vx;
  FILE *fptcls_y, *fptcls_vy;
  fptcls_x=fopen("ptcls_x.dat","w");
  fptcls_y=fopen("ptcls_y.dat","w");
  fptcls_vx=fopen("ptcls_vx.dat","w");
  fptcls_vy=fopen("ptcls_vy.dat","w");

  int c=0;//loop counter
  while(GT<MAXT && Active_Ptcls>0){   
     
    if(ready==Active_Ptcls){ //Everyone's ready 
      //Take a timestep
      c++;
      ready=0;
      GT=c*SampleDelT;

      //Update the global timestep
      //Printout current particle locations
      for(int l=0;l<Active_Ptcls;l++){
	readies[l]=0;
	fprintf(fptcls_x,"%f ",particles[l][X]);
	fprintf(fptcls_y,"%f ",particles[l][Y]);
	fprintf(fptcls_vx,"%f ",ptcls_vx[l]);
	fprintf(fptcls_vy,"%f ",ptcls_vy[l]);
      }
      for(int l=Active_Ptcls;l<Num_Ptcls;l++){
	fprintf(fptcls_x,"-1 ");
	fprintf(fptcls_y,"-1 ");
	fprintf(fptcls_vx,"0 ");
	fprintf(fptcls_vy,"0 ");
      }
      fprintf(fptcls_x,"\n");
      fprintf(fptcls_y,"\n");
      fprintf(fptcls_vx,"\n");
      fprintf(fptcls_vy,"\n");
    }
    for(int i=0;i<Active_Ptcls;i++){
      if(readies[i]!=0){
	continue; //Waiting on other particles to catch up
      }
      interpolate(particles[i][X],particles[i][Y],&p_vx,&p_vy,&elem);
      //Calculate the timestep for this particle based on its speed and elem sz
      if(elem>=0){
	delt=sqrt(Tris[elem].calcArea())/sqrt(p_vx*p_vx+p_vy*p_vy);
	//cout <<" " << delt << " ";
	if(delt+particles[i][T]>SampleDelT){
	  delt=SampleDelT-particles[i][T];
	  particles[i][T]=0;
	} else {
	  particles[i][T]+=delt;
	}
	//Given the velocity and timestep propagate the particle
	//Uses a 4 stage runge-kutta integration method
	k1x=p_vx;
	k1y=p_vy;
        interpolate(particles[i][X]+k1x*delt/2,particles[i][Y]+k1y*delt/2,&p_vx,&p_vy,&elem);
	if(elem>=0){
	  k2x=p_vx;
	  k2y=p_vy;
	  interpolate(particles[i][X]+k2x*delt/2,particles[i][Y]+k2y*delt/2,&p_vx,&p_vy,&elem);
	  if(elem>=0){
	    k3x=p_vx;
	    k3y=p_vy;
	    interpolate(particles[i][X]+k3x*delt,particles[i][Y]+k3y*delt,&p_vx,&p_vy,&elem);
	    if(elem>=0){
	      k4x=p_vx;
	      k4y=p_vy;
	      ptcls_vx[i]=(k1x+2.0*k2x+2.0*k3x+k4x)/6.0;
	      particles[i][X]+=ptcls_vx[i]*delt;
	      ptcls_vy[i]=(k1y+2.0*k2y+2.0*k3y+k4y)/6.0;
	      particles[i][Y]+=ptcls_vy[i]*delt;
	      if(particles[i][T]==0){ 
//only increment readies only if all interpolations succeded!
		readies[i]++;
		ready++;
	      }
	    }
	  }
	}
      } 
      if(elem<0){ //Check if the last interpolation performed has failed
	Active_Ptcls--;
	particles[i][X]=particles[Active_Ptcls][X];
	particles[i][Y]=particles[Active_Ptcls][Y];
	particles[i][T]=particles[Active_Ptcls][T];
	readies[i]=readies[Active_Ptcls];
	ptcls_vx[i]=ptcls_vx[Active_Ptcls];
	ptcls_vy[i]=ptcls_vy[Active_Ptcls];
      }
    }
  }
}
