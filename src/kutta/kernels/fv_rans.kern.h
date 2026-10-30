#include "kutta/kernels/prelude.h"

// UG 8.4.9.1.3 and Eq. 8-88: p_s=p_0-rho*|v_cell|^2/2 on backflow.
KUTTA_KERNEL void fv_backflow_pressure(int i,KIn(int) owner,KIn(int) total,
    KIn(float) mass,KIn(float) u,KIn(float) v,KIn(float) w,
    KOut(float) boundary,int internal,float rho) {
  int f=i+internal,c=owner[f];
  if(total[i]!=0 && mass[f]<0.0f)
    boundary[i]-=0.5f*rho*(u[c]*u[c]+v[c]*v[c]+w[c]*w[c]);
}

KUTTA_KERNEL void fv_effective_viscosity(int i,KIn(int) owner,KIn(int) neighbour,
    KIn(int) kind,KIn(float) weight,KIn(float) mut,KOut(float) mu,int internal,float molecular) {
  int c=owner[i];float turbulent=mut[c];
  if(i<internal) turbulent=weight[i]*turbulent+(1.0f-weight[i])*mut[neighbour[i]];
  mu[i]=i>=internal && kind[i-internal]==0?0.0f:molecular+turbulent;
}

// Transpose of the velocity gradient, contracted with the area vector.
// Wall-function traction replaces BOTH parts of the viscous stress on walls.
KUTTA_KERNEL void fv_transpose_stress(int i,KIn(int) owner,KIn(int) neighbour,
    KIn(float) weight,KIn(float) mu,KIn(float) sx,KIn(float) sy,KIn(float) sz,
    KIn(float) ux,KIn(float) uy,KIn(float) uz,KIn(float) vx,KIn(float) vy,KIn(float) vz,
    KIn(float) wx,KIn(float) wy,KIn(float) wz,KOut(float) x,KOut(float) y,KOut(float) z,int internal) {
  int c=owner[i];float a=ux[c],b=uy[c],d=uz[c],e=vx[c],g=vy[c],h=vz[c],j=wx[c],k=wy[c],l=wz[c];
  if(i<internal) {
    int n=neighbour[i];float t=weight[i],s=1.0f-t;
    a=t*a+s*ux[n];b=t*b+s*uy[n];d=t*d+s*uz[n];
    e=t*e+s*vx[n];g=t*g+s*vy[n];h=t*h+s*vz[n];
    j=t*j+s*wx[n];k=t*k+s*wy[n];l=t*l+s*wz[n];
  }
  x[i]=mu[i]*(a*sx[i]+e*sy[i]+j*sz[i]);
  y[i]=mu[i]*(b*sx[i]+g*sy[i]+k*sz[i]);
  z[i]=mu[i]*(d*sx[i]+h*sy[i]+l*sz[i]);
}
KUTTA_KERNEL void fv_stress_source(int i,KIn(int) offsets,KIn(int) faces,KIn(float) sign,
    KIn(float) volume,KIn(float) x,KIn(float) y,KIn(float) z,
    KOut(float) su,KOut(float) sv,KOut(float) sw) {
  float a=0.0f,b=0.0f,c=0.0f;
  for(int j=offsets[i];j<offsets[i+1];++j) {
    int f=faces[j];a+=sign[j]*x[f];b+=sign[j]*y[f];c+=sign[j]*z[f];
  }
  su[i]+=a/volume[i];sv[i]+=b/volume[i];sw[i]+=c/volume[i];
}
KUTTA_KERNEL void fv_wall_momentum(int i,KIn(int) offsets,KIn(int) faces,KIn(int) kind,
    KIn(float) sx,KIn(float) sy,KIn(float) sz,KIn(float) coefficient,
    KIn(float) u,KIn(float) v,KIn(float) w,KIn(float) bu,KIn(float) bv,KIn(float) bw,
    KOut(float) diagonal,KOut(float) rhs,int internal,int component) {
  float ap=0.0f,b=0.0f;
  for(int j=offsets[i];j<offsets[i+1];++j) {
    int f=faces[j];if(f<internal) continue;
    int bf=f-internal;if(kind[bf]!=0) continue;
    float area=ksqrt(sx[f]*sx[f]+sy[f]*sy[f]+sz[f]*sz[f]);
    float nx=sx[f]/area,ny=sy[f]/area,nz=sz[f]/area;
    float nc=component==0?nx:(component==1?ny:nz);
    float uc=component==0?u[i]:(component==1?v[i]:w[i]);
    float wc=component==0?bu[bf]:(component==1?bv[bf]:bw[bf]);
    float c=area*coefficient[bf];
    ap+=c*(1.0f-nc*nc);
    b+=c*(wc-nc*(nx*bu[bf]+ny*bv[bf]+nz*bw[bf])
      +nc*(nx*u[i]+ny*v[i]+nz*w[i]-nc*uc));
  }
  diagonal[i]=ap;rhs[i]=b;
}
// Trace/3 is rotationally invariant. Both diagonals include wall drag; neither
// physical diagonal includes URF/pseudo time. SIMPLEC averages component inverses.
KUTTA_KERNEL void fv_momentum_diagonal(int i,KIn(float) ax,KIn(float) ay,KIn(float) az,
    KIn(float) ix,KIn(float) iy,KIn(float) iz,KOut(float) physical,KOut(float) inverse) {
  physical[i]=(ax[i]+ay[i]+az[i])/3.0f;inverse[i]=(ix[i]+iy[i]+iz[i])/3.0f;
}
KUTTA_KERNEL void fv_sst_force(int i,KIn(float) selected,KIn(float) sx,KIn(float) sy,KIn(float) sz,
    KIn(float) tx,KIn(float) ty,KIn(float) tz,KOut(float) x,KOut(float) y,KOut(float) z,int internal) {
  int f=i+internal;float a=-selected[i]*ksqrt(sx[f]*sx[f]+sy[f]*sy[f]+sz[f]*sz[f]);
  x[i]=a*tx[i];y[i]=a*ty[i];z[i]=a*tz[i];
}
// Specs: 0=direct k/omega, 1=intensity/length, 2=intensity/viscosity ratio.
// UG 8-75, 8-79 (C_mu^(1/4), C_mu=0.09), 8-80. Intensity is a fraction.
KUTTA_KERNEL void fv_turbulence_boundary(int i,KIn(int) owner,KIn(int) kind,KIn(int) spec,
    KIn(float) a,KIn(float) b,KIn(float) mass,KIn(float) bu,KIn(float) bv,KIn(float) bw,
    KIn(float) k,KIn(float) omega,KIn(float) wall_omega,
    KOut(float) bk,KOut(float) bo,KOut(int) transport_kind,int internal,float rho,float mu) {
  int f=i+internal,c=owner[f],type=kind[i];
  float kinetic=k[c],dissipation=omega[c];transport_kind[i]=1;
  if(type==1 || ((type==2 || type==4) && mass[f]<0.0f)) {
    kinetic=a[i];dissipation=b[i];
    if(spec[i]!=0) {
      kinetic=1.5f*a[i]*a[i]*(bu[i]*bu[i]+bv[i]*bv[i]+bw[i]*bw[i]);
      dissipation=spec[i]==1?ksqrt(kinetic)/(0.5477225575f*b[i]):rho*kinetic/(mu*b[i]);
    }
    transport_kind[i]=0;
  }
  bk[i]=kinetic;bo[i]=type==0?wall_omega[i]:dissipation;
}
