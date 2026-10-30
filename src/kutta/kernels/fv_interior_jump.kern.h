#include "kutta/kernels/prelude.h"

KUTTA_KERNEL void fv_jump_update(int i,KIn(int) faces,KIn(int) kind,KIn(int) count,KIn(float) p,KIn(float) area,
    KIn(float) mass,KOut(float) jump,KOut(float) slope,KOut(float) previous,int n,float rho,float mu) {
  float velocity=mass[faces[i]]/(rho*area[i]);previous[i]=mass[faces[i]];jump[i]=0.0f;slope[i]=0.0f;
  if(kind[i]==1) {
    float raw=kabs(velocity),x=kmax(p[3*n+i],kmin(raw,p[4*n+i]));
    float value=0.0f,derivative=0.0f;
    for(int j=count[i]-1;j>=0;--j) {derivative=derivative*x+value;value=value*x+p[(5+j)*n+i];}
    jump[i]=p[i]*value;
    if(raw>=p[3*n+i] && raw<=p[4*n+i]) slope[i]=p[i]*derivative*(velocity<0.0f?-1.0f:1.0f)/(rho*area[i]);
  }
  if(kind[i]==2) {
    jump[i]=-(mu*p[n+i]+0.5f*rho*p[2*n+i]*kabs(velocity))*velocity;
    slope[i]=-(mu*p[n+i]+rho*p[2*n+i]*kabs(velocity))/(rho*area[i]);
  }
}

KUTTA_KERNEL void fv_jump_pressure_gradient(int i,KIn(int) cells,KIn(int) offsets,KIn(int) incidence,KIn(int) entries,
    KIn(float) dx,KIn(float) dy,KIn(float) dz,KIn(float) inverse,KIn(float) jump,
    KOut(float) gx,KOut(float) gy,KOut(float) gz,int n) {
  int c=cells[i];float x=0.0f,y=0.0f,z=0.0f;
  for(int j=offsets[i];j<offsets[i+1];++j) {
    int entry=entries[j],a=entry>0?entry-1:-entry-1,k=incidence[j];
    float d=(entry>0?-1.0f:1.0f)*jump[a];x+=dx[k]*d;y+=dy[k]*d;z+=dz[k]*d;
  }
  gx[c]+=inverse[c]*x+inverse[n+c]*y+inverse[2*n+c]*z;
  gy[c]+=inverse[n+c]*x+inverse[3*n+c]*y+inverse[4*n+c]*z;
  gz[c]+=inverse[2*n+c]*x+inverse[4*n+c]*y+inverse[5*n+c]*z;
}

KUTTA_KERNEL void fv_jump_momentum(int i,KIn(int) cells,KIn(int) offsets,KIn(int) entries,KIn(int) faces,
    KIn(float) sx,KIn(float) sy,KIn(float) sz,KIn(float) volume,KIn(float) jump,
    KOut(float) u,KOut(float) v,KOut(float) w) {
  int c=cells[i];float x=0.0f,y=0.0f,z=0.0f;
  for(int j=offsets[i];j<offsets[i+1];++j) {
    int entry=entries[j],a=entry>0?entry-1:-entry-1,f=faces[a];float pressure=0.5f*jump[a];
    x+=pressure*sx[f];y+=pressure*sy[f];z+=pressure*sz[f];
  }
  u[c]+=x/volume[c];v[c]+=y/volume[c];w[c]+=z/volume[c];
}

KUTTA_KERNEL void fv_jump_flux(int i,KIn(int) faces,KIn(int) owner,KIn(int) neighbour,KIn(float) volume,
    KIn(float) diagonal,KIn(float) diffusion,KIn(float) jump,KIn(float) slope,KIn(float) previous,
    KOut(float) mass,KOut(float) coefficient,float rho) {
  int f=faces[i],p=owner[f],n=neighbour[f];
  float rc=0.5f*(volume[p]/kmax(diagonal[p],1.0e-20f)+volume[n]/kmax(diagonal[n],1.0e-20f));
  float a=rho*rc*diffusion[f],derivative=kmin(slope[i],0.0f),denominator=1.0f-a*derivative;
  mass[f]=(mass[f]+a*(jump[i]-derivative*previous[i]))/denominator;
  coefficient[f]/=denominator;
}
