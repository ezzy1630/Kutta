#include "kutta/kernels/prelude.h"

KUTTA_KERNEL void fv_gradient(int i, KIn(int) offsets, KIn(int) adjacent,
    KIn(float) dx, KIn(float) dy, KIn(float) dz, KIn(float) inverse,
    KIn(float) cell, KIn(float) boundary,
    KOut(float) gx, KOut(float) gy, KOut(float) gz, int n) {
  float bx=0.0f,by=0.0f,bz=0.0f;
  for(int j=offsets[i];j<offsets[i+1];++j) {
    int other=adjacent[j];
    float d=(other<0 ? boundary[-other-1] : cell[other])-cell[i];
    bx+=dx[j]*d; by+=dy[j]*d; bz+=dz[j]*d;
  }
  gx[i]=inverse[i]*bx+inverse[i+n]*by+inverse[i+2*n]*bz;
  gy[i]=inverse[i+n]*bx+inverse[i+3*n]*by+inverse[i+4*n]*bz;
  gz[i]=inverse[i+2*n]*bx+inverse[i+4*n]*by+inverse[i+5*n]*bz;
}

KUTTA_KERNEL void fv_limiter(int i, KIn(int) offsets, KIn(int) adjacent,
    KIn(float) rx, KIn(float) ry, KIn(float) rz,
    KIn(float) cell, KIn(float) boundary,
    KIn(float) gx, KIn(float) gy, KIn(float) gz, KOut(float) limiter) {
  float lo=cell[i],hi=cell[i];
  for(int j=offsets[i];j<offsets[i+1];++j) {
    int other=adjacent[j];
    float v=other<0 ? boundary[-other-1] : cell[other];
    lo=kmin(lo,v); hi=kmax(hi,v);
  }
  float psi=1.0f;
  for(int j=offsets[i];j<offsets[i+1];++j) {
    float d=gx[i]*rx[j]+gy[i]*ry[j]+gz[i]*rz[j];
    if(d>0.0f) psi=kmin(psi,(hi-cell[i])/d);
    if(d<0.0f) psi=kmin(psi,(lo-cell[i])/d);
  }
  limiter[i]=kmax(0.0f,kmin(1.0f,psi));
}

KUTTA_KERNEL void fv_vector_gradient(int i,KIn(int) offsets,KIn(int) adjacent,
    KIn(float) dx,KIn(float) dy,KIn(float) dz,KIn(float) inverse,
    KIn(float) u,KIn(float) v,KIn(float) w,KIn(float) bu,KIn(float) bv,KIn(float) bw,
    KOut(float) ux,KOut(float) uy,KOut(float) uz,KOut(float) vx,KOut(float) vy,KOut(float) vz,
    KOut(float) wx,KOut(float) wy,KOut(float) wz,int n) {
  float ax=0.0f,ay=0.0f,az=0.0f,bx=0.0f,by=0.0f,bz=0.0f,cx=0.0f,cy=0.0f,cz=0.0f;
  for(int j=offsets[i];j<offsets[i+1];++j) {
    int other=adjacent[j];
    float a=(other<0?bu[-other-1]:u[other])-u[i];
    float b=(other<0?bv[-other-1]:v[other])-v[i];
    float c=(other<0?bw[-other-1]:w[other])-w[i];
    ax+=dx[j]*a;ay+=dy[j]*a;az+=dz[j]*a;
    bx+=dx[j]*b;by+=dy[j]*b;bz+=dz[j]*b;
    cx+=dx[j]*c;cy+=dy[j]*c;cz+=dz[j]*c;
  }
  float xx=inverse[i],xy=inverse[i+n],xz=inverse[i+2*n],yy=inverse[i+3*n],yz=inverse[i+4*n],zz=inverse[i+5*n];
  ux[i]=xx*ax+xy*ay+xz*az;uy[i]=xy*ax+yy*ay+yz*az;uz[i]=xz*ax+yz*ay+zz*az;
  vx[i]=xx*bx+xy*by+xz*bz;vy[i]=xy*bx+yy*by+yz*bz;vz[i]=xz*bx+yz*by+zz*bz;
  wx[i]=xx*cx+xy*cy+xz*cz;wy[i]=xy*cx+yy*cy+yz*cz;wz[i]=xz*cx+yz*cy+zz*cz;
}
