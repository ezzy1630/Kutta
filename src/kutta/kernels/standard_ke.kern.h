#include "kutta/kernels/prelude.h"

// Standard k-epsilon (Fluent TG 4.3.1, Launder-Spalding) with standard wall
// functions (TG "Wall Treatment for epsilon-based Models"). The realizable
// model's transport, boundary and wall-cell kernels are shared; only the
// eddy viscosity, the epsilon source and the wall law differ.

// mu_t = rho Cmu k^2/epsilon (TG 4-50), constant Cmu. strain is |S|=sqrt(2SijSij).
KUTTA_KERNEL void ske_viscosity(int i,KIn(float) rho,KIn(float) mu,KIn(float) k,KIn(float) epsilon,
    KIn(float) ux,KIn(float) uy,KIn(float) uz,KIn(float) vx,KIn(float) vy,KIn(float) vz,
    KIn(float) wx,KIn(float) wy,KIn(float) wz,KOut(float) mut,KOut(float) cmu,KOut(float) strain,
    float c_mu,float k_floor,float e_floor,float max_ratio) {
  float a=ux[i],d=vy[i],f=wz[i],b=0.5f*(uy[i]+vx[i]),c=0.5f*(uz[i]+wx[i]),e=0.5f*(vz[i]+wy[i]);
  float ss=a*a+d*d+f*f+2.0f*(b*b+c*c+e*e);
  float kinetic=kmax(k[i],k_floor),eps=kmax(epsilon[i],e_floor);
  cmu[i]=c_mu;mut[i]=kmin(rho[i]*c_mu*kinetic*(kinetic/eps),max_ratio*mu[i]);strain[i]=1.4142135624f*ksqrt(ss);
}

// TG 4-48/4-49 sources in source-sink*phi form: k gains G_k and loses rho*eps;
// epsilon gains C1 (eps/k) G_k and loses C2 rho eps^2/k. Wall-adjacent cells
// take the wall-function production and their epsilon row is fixed.
KUTTA_KERNEL void ske_sources(int i,KIn(float) rho,KIn(float) mu,KIn(float) k,KIn(float) epsilon,
    KIn(float) mut,KIn(float) strain,KIn(int) wall_mask,KIn(float) wall_prod,
    KIn(float) wall_epsilon,KOut(float) sk,KOut(float) dk,KOut(float) se,KOut(float) de,
    KOut(float) gamma_k,KOut(float) gamma_e,float c1,float c2,float sigk,float sige,
    float molecular_prod,float k_floor,float e_floor) {
  float kinetic=kmax(k[i],k_floor),eps=kmax(epsilon[i],e_floor),s=strain[i];
  float production=(mut[i]+molecular_prod*mu[i])*s*s;
  if(wall_mask[i]) {production=wall_prod[i];eps=wall_epsilon[i];}
  sk[i]=production;dk[i]=rho[i]*eps/kinetic;
  se[i]=c1*(eps/kinetic)*production;
  de[i]=rho[i]*c2*eps/kinetic;
  gamma_k[i]=mu[i]+mut[i]/sigk;gamma_e[i]=mu[i]+mut[i]/sige;
}

// Standard wall functions. Velocity: log law U* = ln(E y*)/kappa for
// y* > 11.225, else the laminar U* = y*. G_k and epsilon at the wall cell use
// the equilibrium formulas with the actual y_P for every y* (unlike scalable
// wall functions, which replace y* by max(y*,11.225) everywhere).
KUTTA_KERNEL void ske_standard_wall(int i,KIn(int) owner,KIn(int) wall,KIn(float) dy,
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
  float y=rho[c]/mu[c]*us*dy[i];
  // Laminar branch: rho u*/y* = mu/y_P, which stays defined as k -> 0.
  float coeff=y>11.225f?rho[c]*us*0.4187f/klog(9.793f*y):mu[c]/dy[i],stress=coeff*speed;
  float dissipation=us*us*us/(0.4187f*dy[i]);
  coefficient[i]=coeff;ystar[i]=y;tau[i]=stress;epsilon[i]=kmax(dissipation,e_floor);
  production[i]=stress*stress/(0.4187f*rho[c]*us*dy[i]);
  tx[i]=-coeff*tangent.x;ty[i]=-coeff*tangent.y;tz[i]=-coeff*tangent.z;
}
