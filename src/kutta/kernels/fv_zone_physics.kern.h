#include "kutta/kernels/prelude.h"

KUTTA_KERNEL void fv_frame_flux(int i,KIn(float) frame,KIn(int) kind,KOut(float) mass,int internal) {
  if(i>=internal && (kind[i-internal]==0 || kind[i-internal]==3)) return;
  mass[i]-=frame[i];
}

// TG 2-10: div(rho*u_relative*u_absolute) + rho*Omega x (u - frame_translation).
// Porous: -(mu*D + rho*|u|*C/2)*u. A scalar spectral bound is implicit;
// its matching explicit correction retains every off-axis tensor term.
KUTTA_KERNEL void fv_zone_momentum(int i,KIn(int) zone,KIn(float) a,KIn(float) reference,KIn(float) u,KIn(float) v,KIn(float) w,
    KIn(float) volume,KOut(float) diagonal,KOut(float) rhs,int cells,int n,int component,int moving_resistance,float rho,float mu) {
  int c=zone[i];if(c<0) return;
  float x=u[i],y=v[i],z=w[i];
  if(moving_resistance!=0) {x-=reference[i];y-=reference[cells+i];z-=reference[2*cells+i];}
  float inertia=0.5f*rho*ksqrt(x*x+y*y+z*z);
  float xx=mu*a[c]+inertia*a[6*n+c],xy=mu*a[n+c]+inertia*a[7*n+c];
  float xz=mu*a[2*n+c]+inertia*a[8*n+c],yy=mu*a[3*n+c]+inertia*a[9*n+c];
  float yz=mu*a[4*n+c]+inertia*a[10*n+c],zz=mu*a[5*n+c]+inertia*a[11*n+c];
  float sink=mu*a[18*n+c]+inertia*a[19*n+c];
  float ux=u[i]-a[15*n+c],uy=v[i]-a[16*n+c],uz=w[i]-a[17*n+c];
  float force=0.0f,q=0.0f;
  if(component==0) {force=-(xx*x+xy*y+xz*z)-rho*(a[13*n+c]*uz-a[14*n+c]*uy);q=u[i];}
  if(component==1) {force=-(xy*x+yy*y+yz*z)-rho*(a[14*n+c]*ux-a[12*n+c]*uz);q=v[i];}
  if(component==2) {force=-(xz*x+yz*y+zz*z)-rho*(a[12*n+c]*uy-a[13*n+c]*ux);q=w[i];}
  diagonal[i]+=volume[i]*sink;rhs[i]+=volume[i]*(force+sink*q);
}
