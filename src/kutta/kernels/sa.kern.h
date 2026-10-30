#pragma once
#include "kutta/kernels/prelude.h"

KUTTA_FN float sa_fv1(float chi,float cv1) {
  // Ratio form avoids chi^3 overflow. No eddy-viscosity cap changes damping.
  float t=kmin(chi/cv1,1e12f),t3=t*t*t;
  return t3/(1.0f+t3);
}
KUTTA_KERNEL void sa_coefficients(int i,KIn(float) rho,KIn(float) mu,KIn(float) q,
    KIn(float) distance,KIn(float) strain,KIn(float) rotation,
    KOut(float) mut,KOut(float) gamma,KOut(float) fv1,KOut(float) fv2,
    KOut(float) rate,KOut(float) fw,int alternate,int tmr,
    float sigma,float cv1,float cw2,float cw3,float kappa,float max_ratio) {
  float v=kmax(q[i],0.0f),chi=v/(mu[i]/rho[i]);
  float t=kmin(chi/cv1,1e12f),one_minus_f1=1.0f/(1.0f+t*t*t);
  float f1=sa_fv1(chi,cv1),f2=(1.0f-chi*one_minus_f1)/(1.0f+chi*f1);
  float base=rotation[i]+(alternate?2.0f*kmin(0.0f,strain[i]-rotation[i]):0.0f);
  float invd=1.0f/distance[i],wall_rate=(v*invd)*invd/(kappa*kappa);
  float s=base+wall_rate*f2;
  // TG leaves the sign of modified S intact. A negative production or
  // destruction is split by sign in sa_sources, never made into a negative sink.
  // NASA's standard implementation bounds S below by .3*Omega and r above by 10.
  if(tmr) s=kmax(s,0.3f*rotation[i]);
  float r=wall_rate==0.0f?0.0f:(s==0.0f?10.0f:wall_rate/s);
  r=kmax(-1e4f,kmin(r,tmr?10.0f:1e4f));
  float r2=r*r,r6=r2*r2*r2,g=r+cw2*(r6-r);
  float scale=kmax(kabs(g),cw3),a=g/scale,b=cw3/scale;
  float a2=a*a,b2=b*b;
  float c6=kpow(cw3,6.0f);
  float wall=a*kpow((1.0f+c6)/(a2*a2*a2+b2*b2*b2),1.0f/6.0f);
  mut[i]=kmin(rho[i]*v*f1,max_ratio*mu[i]);
  gamma[i]=(mu[i]+rho[i]*v)/sigma;
  fv1[i]=f1;fv2[i]=f2;rate[i]=s;fw[i]=wall;
}
KUTTA_KERNEL void sa_sources(int i,KIn(float) rho,KIn(float) mu,KIn(float) q,
    KIn(float) distance,KIn(float) gx,KIn(float) gy,KIn(float) gz,
    KIn(float) rate,KIn(float) fw,KOut(float) source,KOut(float) sink,
    int tmr,float cb1,float cb2,float sigma,float kappa) {
  float v=kmax(q[i],0.0f),chi=v/(mu[i]/rho[i]);
  float ft2=tmr?1.2f*kexp(-0.5f*kmin(chi,100.0f)*kmin(chi,100.0f)):0.0f;
  float cw1=cb1/(kappa*kappa)+(1.0f+cb2)/sigma;
  float p=cb1*(1.0f-ft2)*rate[i];
  float invd=1.0f/distance[i];
  float destruction=(cw1*fw[i]-cb1*ft2/(kappa*kappa))*(v*invd)*invd;
  float growth=p-destruction;
  source[i]=rho[i]*(kmax(growth,0.0f)*v+cb2/sigma*(gx[i]*gx[i]+gy[i]*gy[i]+gz[i]*gz[i]));
  sink[i]=rho[i]*kmax(-growth,0.0f);
}
KUTTA_KERNEL void sa_face_diffusion(int i,KIn(int) owner,KIn(int) neighbour,
    KIn(float) weight,KIn(int) wall,KIn(float) mu,KIn(float) gamma,
    KOut(float) face,int internal,float sigma) {
  int c=owner[i];
  face[i]=i<internal?weight[i]*gamma[c]+(1.0f-weight[i])*gamma[neighbour[i]]:
    (wall[i-internal]?mu[c]/sigma:gamma[c]);
}
// Monotone inversion of chi*fv1(chi)=mu_t/mu. Direct nu_tilde is a different BC.
KUTTA_FN float sa_ratio_to_chi(float ratio,float cv1) {
  if(ratio==0.0f) return 0.0f;
  float lo=0.0f,hi=ratio+cv1;
  for(int j=0;j<40;++j) {float m=0.5f*(lo+hi);if(m*sa_fv1(m,cv1)<ratio) lo=m;else hi=m;}
  return 0.5f*(lo+hi);
}
KUTTA_KERNEL void sa_boundary(int i,KIn(int) owner,KIn(int) kind,KIn(int) spec,
    KIn(float) value,KIn(float) length,KIn(float) u,KIn(float) v,KIn(float) w,
    KIn(float) mass,KIn(float) q,KOut(float) boundary,KOut(int) transport,
    int internal,float nu,float cv1) {
  int f=i+internal,c=owner[f],type=kind[i];
  boundary[i]=q[c];transport[i]=1;
  if(type==0) {boundary[i]=0.0f;transport[i]=0;}
  else if(type==1 || ((type==2 || type==4) && mass[f]<0.0f)) {
    boundary[i]=spec[i]==0?value[i]:nu*sa_ratio_to_chi(value[i],cv1);
    // UG 7-74 assumes fv1=1 for intensity/length, unlike ratio inversion.
    if(spec[i]==2) boundary[i]=0.6708203932f*value[i]*length[i]*ksqrt(u[i]*u[i]+v[i]*v[i]+w[i]*w[i]);
    transport[i]=0;
  }
}
KUTTA_KERNEL void sa_momentum_viscosity(int i,KIn(int) owner,KIn(int) neighbour,KIn(int) kind,
    KIn(float) weight,KIn(float) mut,KOut(float) mu,int internal,float molecular) {
  float turbulent=mut[owner[i]];
  if(i<internal) turbulent=weight[i]*turbulent+(1.0f-weight[i])*mut[neighbour[i]];
  mu[i]=molecular+(i>=internal && kind[i-internal]==0?0.0f:turbulent);
}
// Input: force on body. Output: stress on fluid and friction-velocity y+.
KUTTA_KERNEL void sa_wall_diagnostics(int i,KIn(float) ax,KIn(float) ay,KIn(float) az,
    KIn(float) distance,KOut(float) x,KOut(float) y,KOut(float) z,KOut(float) yp,
    int internal,float rho,float mu) {
  int f=i+internal;float area=ksqrt(ax[f]*ax[f]+ay[f]*ay[f]+az[f]*az[f]);
  x[i]=-x[i]/area;y[i]=-y[i]/area;z[i]=-z[i]/area;
  float tau=ksqrt(x[i]*x[i]+y[i]*y[i]+z[i]*z[i]);
  yp[i]=rho/mu*distance[i]*ksqrt(tau/rho);
}
