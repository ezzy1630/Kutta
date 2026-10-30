#include "kutta/kernels/prelude.h"

// Boundary kind: zero=Dirichlet, one=prescribed normal derivative. The latter
// supplies its flux through the gathered cell source, not through diffusion.
KUTTA_KERNEL void fv_hybrid_boundary(int i,KIn(int) owner,KIn(int) kind,
    KIn(float) cell,KIn(float) prescribed,KIn(float) normal_delta,
    KOut(float) boundary,int internal) {
  boundary[i]=kind[i]==0?prescribed[i]:cell[owner[i+internal]]+normal_delta[i];
}

KUTTA_KERNEL void fv_hybrid_velocity(int i,KIn(float) gx,KIn(float) gy,KIn(float) gz,
    KOut(float) u,KOut(float) v,KOut(float) w,float ux,float uy,float uz,float speed,int constant_speed) {
  float x=ux+gx[i],y=uy+gy[i],z=uz+gz[i],m=ksqrt(x*x+y*y+z*z);
  float scale=constant_speed!=0 && m>0.0f?speed/m:1.0f;
  u[i]=scale*x;v[i]=scale*y;w[i]=scale*z;
}

KUTTA_KERNEL void fv_hybrid_flux(int i,KIn(int) owner,KIn(int) neighbour,KIn(int) kind,
    KIn(float) potential,KIn(float) diffusion,KIn(float) nonorthogonal,KIn(float) prescribed_flux,
    KIn(float) boundary,KIn(float) sx,KIn(float) sy,KIn(float) sz,
    KOut(float) mass,int internal,float rho,float ux,float uy,float uz) {
  int p=owner[i];
  float base=ux*sx[i]+uy*sy[i]+uz*sz[i];
  if(i<internal) mass[i]=rho*(base+diffusion[i]*(potential[neighbour[i]]-potential[p])+nonorthogonal[i]);
  else mass[i]=rho*(kind[i-internal]==0?base+diffusion[i]*(boundary[i-internal]-potential[p])+nonorthogonal[i]:prescribed_flux[i-internal]);
}

// Area-averaged inlet recipes are represented as nonnegative polynomials of
// local speed. Averaging their volume-weighted outputs is a separate reduction.
KUTTA_KERNEL void fv_hybrid_turbulence(int i,KIn(float) u,KIn(float) v,KIn(float) w,
    KIn(float) volume,KOut(float) k,KOut(float) omega,KOut(float) vk,KOut(float) vo,
    float k0,float k2,float w0,float w1,float w2,float kfloor,float wfloor) {
  float speed2=u[i]*u[i]+v[i]*v[i]+w[i]*w[i];
  k[i]=kmax(kfloor,k0+k2*speed2);
  omega[i]=kmax(wfloor,w0+w1*ksqrt(speed2)+w2*speed2);
  vk[i]=volume[i]*k[i];vo[i]=volume[i]*omega[i];
}
