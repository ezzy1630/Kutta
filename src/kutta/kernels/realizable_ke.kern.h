#pragma once
#include "kutta/kernels/prelude.h"

// cos(acos(x)/3) is the root of 4z^3-3z=x on [.5,1]. Bounded bisection
// needs no backend-specific transcendental wrappers and is float-accurate.
KUTTA_FN float rke_cos_third(float x) {
  float lo=0.5f,hi=1.0f;
  for(int j=0;j<24;++j) {float z=0.5f*(lo+hi);if(z*(4.0f*z*z-3.0f)<x) lo=z;else hi=z;}
  return 0.5f*(lo+hi);
}
KUTTA_KERNEL void rke_viscosity(int i,KIn(float) rho,KIn(float) mu,KIn(float) k,KIn(float) epsilon,
    KIn(float) ux,KIn(float) uy,KIn(float) uz,KIn(float) vx,KIn(float) vy,KIn(float) vz,
    KIn(float) wx,KIn(float) wy,KIn(float) wz,KOut(float) mut,KOut(float) cmu,KOut(float) strain,
    float a0,float k_floor,float e_floor,float max_ratio) {
  float a=ux[i],d=vy[i],f=wz[i],b=0.5f*(uy[i]+vx[i]),c=0.5f*(uz[i]+wx[i]),e=0.5f*(vz[i]+wy[i]);
  float ss=a*a+d*d+f*f+2.0f*(b*b+c*c+e*e),s=ksqrt(ss);
  // W = tr(S^3)/|S|^3. Normalize before cubing to retain tiny gradients.
  float invariant=0.0f;
  if(s>0.0f) {
    a/=s;b/=s;c/=s;d/=s;e/=s;f/=s;
    invariant=a*a*a+d*d*d+f*f*f+3.0f*((a+d)*b*b+(a+f)*c*c+(d+f)*e*e)+6.0f*b*c*e;
  }
  float as=2.4494897428f*rke_cos_third(kmax(-1.0f,kmin(1.0f,2.4494897428f*invariant)));
  float xy=uy[i]-vx[i],xz=uz[i]-wx[i],yz=vz[i]-wy[i];
  float ustar=ksqrt(ss+0.5f*(xy*xy+xz*xz+yz*yz));
  float kinetic=kmax(k[i],k_floor),eps=kmax(epsilon[i],e_floor);
  float coefficient=1.0f/(a0+as*(kinetic/eps)*ustar);
  cmu[i]=coefficient;mut[i]=kmin(rho[i]*coefficient*kinetic*(kinetic/eps),max_ratio*mu[i]);strain[i]=1.4142135624f*s;
}
KUTTA_KERNEL void rke_sources(int i,KIn(float) rho,KIn(float) mu,KIn(float) k,KIn(float) epsilon,
    KIn(float) mut,KIn(float) strain,KIn(int) wall_mask,KIn(float) wall_prod,
    KIn(float) wall_epsilon,KOut(float) sk,KOut(float) dk,KOut(float) se,KOut(float) de,
    KOut(float) gamma_k,KOut(float) gamma_e,float c2,float sigk,float sige,
    float molecular_prod,float k_floor,float e_floor) {
  float kinetic=kmax(k[i],k_floor),eps=kmax(epsilon[i],e_floor),s=strain[i];
  float eta=s*(kinetic/eps),c1=kmax(0.43f,eta/(eta+5.0f));
  float production=(mut[i]+molecular_prod*mu[i])*s*s;
  if(wall_mask[i]) {production=wall_prod[i];eps=wall_epsilon[i];}
  sk[i]=production;dk[i]=rho[i]*eps/kinetic;
  se[i]=rho[i]*c1*s*eps;
  de[i]=rho[i]*c2*eps/(kinetic+ksqrt(mu[i]/rho[i]*eps));
  gamma_k[i]=mu[i]+mut[i]/sigk;gamma_e[i]=mu[i]+mut[i]/sige;
}
KUTTA_KERNEL void rke_scalable_wall(int i,KIn(int) owner,KIn(int) wall,KIn(float) dy,
    KIn(float) ax,KIn(float) ay,KIn(float) az,KIn(float) u,KIn(float) v,KIn(float) w,
    KIn(float) bu,KIn(float) bv,KIn(float) bw,KIn(float) rho,KIn(float) mu,KIn(float) k,
    KOut(float) coefficient,KOut(float) ystar,KOut(float) tau,KOut(float) epsilon,
    KOut(float) production,KOut(float) tx,KOut(float) ty,KOut(float) tz,int internal,float k_floor,float e_floor) {
  coefficient[i]=0.0f;ystar[i]=0.0f;tau[i]=0.0f;epsilon[i]=0.0f;production[i]=0.0f;tx[i]=0.0f;ty[i]=0.0f;tz[i]=0.0f;
  if(!wall[i]) return;
  int c=owner[i+internal],f=i+internal;
  kvec3 normal=kv3(ax[f],ay[f],az[f]);normal=(1.0f/knorm(normal))*normal;
  kvec3 relative=kv3(u[c]-bu[i],v[c]-bv[i],w[c]-bw[i]);
  kvec3 tangent=relative-kdot(relative,normal)*normal;
  float speed=knorm(tangent),us=0.5477225575f*ksqrt(kmax(k[c],k_floor));
  float y=rho[c]/mu[c]*us*dy[i],effective=kmax(y,11.225f);
  float up=klog(9.793f*effective)/0.4187f;
  float coeff=rho[c]*us/up,stress=coeff*speed;
  // Replace y* in EVERY wall-function formula, TG 4-357, not only U+.
  float effective_distance=effective*mu[c]/(rho[c]*us);
  float dissipation=us*us*us/(0.4187f*effective_distance);
  coefficient[i]=coeff;ystar[i]=y;tau[i]=stress;epsilon[i]=kmax(dissipation,e_floor);
  production[i]=stress*stress/(0.4187f*rho[c]*us*effective_distance);
  tx[i]=-coeff*tangent.x;ty[i]=-coeff*tangent.y;tz[i]=-coeff*tangent.z;
}
KUTTA_KERNEL void rke_wall_cells(int i,KIn(int) offsets,KIn(int) faces,KIn(int) wall,
    KIn(float) distance,KIn(float) epsilon,KIn(float) production,KOut(int) mask,
    KOut(int) selected,KOut(float) fixed,KOut(float) prod,int internal) {
  int best=-1;float nearest=1e38f;
  for(int j=offsets[i];j<offsets[i+1];++j) {
    int b=faces[j]-internal;
    if(b>=0 && wall[b] && (distance[b]<nearest || (distance[b]==nearest && (best<0 || b<best)))) {best=b;nearest=distance[b];}
  }
  selected[i]=best;mask[i]=best>=0?1:0;fixed[i]=best>=0?epsilon[best]:0.0f;prod[i]=best>=0?production[best]:0.0f;
}
KUTTA_KERNEL void rke_face_diffusion(int i,KIn(int) owner,KIn(int) neighbour,KIn(float) weight,
    KIn(int) wall,KIn(float) mu,KIn(float) gk,KIn(float) ge,KOut(float) fk,KOut(float) fe,int internal) {
  int c=owner[i];float k=gk[c],e=ge[c];
  if(i<internal) {float a=weight[i];k=a*k+(1.0f-a)*gk[neighbour[i]];e=a*e+(1.0f-a)*ge[neighbour[i]];}
  else if(wall[i-internal]) {k=0.0f;e=mu[c];}
  fk[i]=k;fe[i]=e;
}
// TG 4-364..373. These coefficients alone do not implement EWT transport.
KUTTA_FN float rke_one_minus_exp(float x) {
  return x<0.05f?x*(1.0f-x*(0.5f-x*(1.0f/6.0f-x*(1.0f/24.0f-x/120.0f)))):1.0f-kexp(-x);
}
KUTTA_KERNEL void rke_ewt_algebra(int i,KIn(float) rho,KIn(float) mu,KIn(float) k,
    KIn(float) distance,KIn(float) outer,KOut(float) viscosity,KOut(float) epsilon,
    KOut(float) weight,float width,float cmu) {
  float y=distance[i],root=ksqrt(kmax(k[i],0.0f)),re=rho[i]*y*root/mu[i];
  if(1.0f/y==0.0f) {weight[i]=1.0f;viscosity[i]=outer[i];epsilon[i]=0.0f;return;}
  float cl=0.4187f*kpow(cmu,-0.75f),ae=2.0f*cl;
  float lm=y*cl*rke_one_minus_exp(re/70.0f);
  float le=y*cl*rke_one_minus_exp(re/ae);
  float x=(re-200.0f)/(width/2.297559925f);
  float tail=kexp(-2.0f*kabs(x));
  float blend=x>=0.0f?1.0f/(1.0f+tail):tail/(1.0f+tail);
  weight[i]=blend;viscosity[i]=blend*outer[i]+(1.0f-blend)*rho[i]*cmu*lm*root;
  // Small-Re analytic limit avoids 0/0 when k vanishes.
  float ratio=re/ae,correction=ratio>0.0f?ratio/rke_one_minus_exp(ratio):1.0f;
  epsilon[i]=ratio<0.05f?2.0f*(mu[i]/rho[i])*kmax(k[i],0.0f)/(y*y)*correction:root*root*root/le;
}
// method: 0=k/epsilon; 1=intensity/length; 2=intensity/viscosity ratio;
// 3=intensity/hydraulic diameter. Boundary Cmu=.09 is empirical (UG 7-76/77).
KUTTA_KERNEL void rke_boundary(int i,KIn(int) owner,KIn(int) kind,KIn(int) method,
    KIn(float) a,KIn(float) b,KIn(float) mass,KIn(float) u,KIn(float) v,KIn(float) w,
    KIn(float) k,KIn(float) epsilon,KOut(float) bk,KOut(float) be,KOut(int) transport,
    int internal,float rho,float mu) {
  int f=i+internal,c=owner[f];float kinetic=k[c],eps=epsilon[c];transport[i]=1;
  if(kind[i]==1 || ((kind[i]==2 || kind[i]==4) && mass[f]<0.0f)) {
    kinetic=a[i];eps=b[i];transport[i]=0;
    if(method[i]!=0) {
      kinetic=1.5f*a[i]*a[i]*(u[i]*u[i]+v[i]*v[i]+w[i]*w[i]);
      float length=method[i]==3?0.07f*b[i]:b[i];
      eps=method[i]==2?rho*0.09f*kinetic*(kinetic/(mu*b[i])):0.1643167673f*kinetic*ksqrt(kinetic)/length;
    }
  }
  bk[i]=kinetic;be[i]=eps;
}
// Fully developed channel: streamwise velocity is tangential to both walls.
// Assemble the same implicit wall traction as fv_boundary_momentum, scalar form.
KUTTA_KERNEL void rke_channel_wall(int i,KIn(int) offsets,KIn(int) faces,KIn(int) wall,
    KIn(float) ax,KIn(float) ay,KIn(float) az,KIn(float) coefficient,KOut(float) diagonal,int internal) {
  float value=0.0f;
  for(int j=offsets[i];j<offsets[i+1];++j) {
    int f=faces[j],b=f-internal;
    if(b>=0 && wall[b]) value+=coefficient[b]*ksqrt(ax[f]*ax[f]+ay[f]*ay[f]+az[f]*az[f]);
  }
  diagonal[i]=value;
}
// TG 4-374..376, zero pressure-gradient / thermal effects only.
KUTTA_KERNEL void rke_ewt_law(int i,KIn(float) yp,KOut(float) up,KOut(float) derivative) {
  float y=yp[i];
  if(y<=0.0f) {up[i]=0.0f;derivative[i]=1.0f;return;}
  float gamma=-0.01f*y*y*y*y/(1.0f+5.0f*y);
  float lam=kexp(gamma),turb=kexp(1.0f/gamma);
  up[i]=lam*y+turb*klog(9.793f*y)/0.4187f;
  derivative[i]=lam+turb/(0.4187f*y);
}
