#include "kutta/kernels/prelude.h"

KUTTA_KERNEL void fv_pressure_inlet_velocity(int i,KIn(int) kind,
    KIn(float) sx,KIn(float) sy,KIn(float) sz,KIn(float) dx,KIn(float) dy,KIn(float) dz,
    KIn(float) mass,KOut(float) u,KOut(float) v,KOut(float) w,int internal,float rho) {
  if(kind[i]!=1 || mass[i+internal]>=0.0f) return;
  int f=i+internal;
  float projection=dx[i]*sx[f]+dy[i]*sy[f]+dz[i]*sz[f];
  float speed=mass[f]/(rho*projection);
  u[i]=speed*dx[i];v[i]=speed*dy[i];w[i]=speed*dz[i];
}

KUTTA_KERNEL void fv_pressure_inlet_pressure(int i,KIn(int) kind,KIn(float) mass,
    KIn(float) total,KIn(float) u,KIn(float) v,KIn(float) w,KOut(float) pressure,int internal,float rho) {
  if(kind[i]==1 && mass[i+internal]<0.0f)
    pressure[i]=total[i]-0.5f*rho*(u[i]*u[i]+v[i]*v[i]+w[i]*w[i]);
}

// Existing wall projection is in the stationary frame. Restore the frame's
// normal velocity so only motion relative to the adjacent zone is projected.
KUTTA_KERNEL void fv_frame_wall_velocity(int i,KIn(int) kind,KIn(float) normal_flux,
    KIn(float) sx,KIn(float) sy,KIn(float) sz,KOut(float) u,KOut(float) v,KOut(float) w,
    int internal,float rho) {
  if(kind[i]!=0) return;
  int f=i+internal;float a2=sx[f]*sx[f]+sy[f]*sy[f]+sz[f]*sz[f];
  float normal=normal_flux[f]/(rho*a2);
  u[i]+=normal*sx[f];v[i]+=normal*sy[f];w[i]+=normal*sz[f];
}
