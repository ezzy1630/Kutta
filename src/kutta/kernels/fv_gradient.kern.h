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
