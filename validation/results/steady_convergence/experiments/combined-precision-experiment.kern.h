#include "kutta/kernels/prelude.h"

struct kpressure_pair {float hi,lo;};
KUTTA_FN kpressure_pair kpressure_add(kpressure_pair a,kpressure_pair b) {
  float s=a.hi+b.hi,v=s-a.hi;
  float e=((a.hi-(s-v))+(b.hi-v))+a.lo+b.lo;
  float h=s+e;kpressure_pair r={h,e-(h-s)};return r;
}
KUTTA_FN kpressure_pair kpressure_product(float a,float b) {
  float h=a*b;kpressure_pair r={h,kfma(a,b,-h)};return r;
}
KUTTA_FN kpressure_pair kpressure_scale(kpressure_pair a,float b) {
  auto p=kpressure_product(a.hi,b);kpressure_pair low={a.lo*b,0};return kpressure_add(p,low);
}
KUTTA_FN kpressure_pair kpressure_reconstruct(float p,float low,float gx,float gy,float gz,
    float x,float y,float z) {
  kpressure_pair a={p,low};a=kpressure_add(a,kpressure_product(gx,x));
  a=kpressure_add(a,kpressure_product(gy,y));return kpressure_add(a,kpressure_product(gz,z));
}

// Same LSQ observations and matrix, retaining sub-ulp pressure increments.
KUTTA_KERNEL void fv_pressure_gradient_precise(int i,KIn(int) offsets,KIn(int) adjacent,
    KIn(float) dx,KIn(float) dy,KIn(float) dz,KIn(float) inverse,
    KIn(float) cell,KIn(float) low,KIn(float) boundary,
    KOut(float) gx,KOut(float) gy,KOut(float) gz,int n) {
  kpressure_pair bx={0,0},by={0,0},bz={0,0};
  for(int j=offsets[i];j<offsets[i+1];++j) {
    int other=adjacent[j];kpressure_pair a={other<0?boundary[-other-1]:cell[other],other<0?0.0f:low[other]},b={-cell[i],-low[i]};
    auto d=kpressure_add(a,b);float delta=d.hi+d.lo;
    bx=kpressure_add(bx,kpressure_product(dx[j],delta));by=kpressure_add(by,kpressure_product(dy[j],delta));bz=kpressure_add(bz,kpressure_product(dz[j],delta));
  }
  gx[i]=inverse[i]*(bx.hi+bx.lo)+inverse[i+n]*(by.hi+by.lo)+inverse[i+2*n]*(bz.hi+bz.lo);
  gy[i]=inverse[i+n]*(bx.hi+bx.lo)+inverse[i+3*n]*(by.hi+by.lo)+inverse[i+4*n]*(bz.hi+bz.lo);
  gz[i]=inverse[i+2*n]*(bx.hi+bx.lo)+inverse[i+4*n]*(by.hi+by.lo)+inverse[i+5*n]*(bz.hi+bz.lo);
}

// Algebraically identical to fv_pressure_faces; only cancellation precision
// changes. The pressure face interpolation itself remains the original kernel.
KUTTA_KERNEL void fv_pressure_jump_precise(int i,KIn(int) owner,KIn(int) neighbour,
    KIn(float) px,KIn(float) py,KIn(float) pz,KIn(float) nx,KIn(float) ny,KIn(float) nz,
    KIn(float) p,KIn(float) low,KIn(float) gx,KIn(float) gy,KIn(float) gz,KIn(float) boundary,
    KOut(float) jump,int internal) {
  int c=owner[i];auto pp=kpressure_reconstruct(p[c],low[c],gx[c],gy[c],gz[c],px[i],py[i],pz[i]);
  if(i<internal) {
    int n=neighbour[i];auto pn=kpressure_reconstruct(-p[n],-low[n],-gx[n],-gy[n],-gz[n],nx[i],ny[i],nz[i]);
    auto difference=kpressure_add(pp,pn);float side=difference.hi+difference.lo;
    float dx=px[i]-nx[i],dy=py[i]-ny[i],dz=pz[i]-nz[i];
    float sx=0.5f*(px[i]+nx[i]),sy=0.5f*(py[i]+ny[i]),sz=0.5f*(pz[i]+nz[i]);
    float d2=dx*dx+dy*dy+dz*dz,s2=sx*sx+sy*sy+sz*sz;
    float trust=kmax(0.0f,kmin(1.0f,(4.0f*d2-s2)/(3.0f*d2)));
    kpressure_pair pc={p[c],low[c]},negative={-p[n],-low[n]};auto shared=kpressure_add(pc,negative);
    shared=kpressure_add(shared,kpressure_product(0.5f*(gx[c]+gx[n]),dx));
    shared=kpressure_add(shared,kpressure_product(0.5f*(gy[c]+gy[n]),dy));
    shared=kpressure_add(shared,kpressure_product(0.5f*(gz[c]+gz[n]),dz));
    float link=shared.hi+shared.lo;jump[i]=trust==1.0f?side:link+trust*(side-link);
  } else {
    kpressure_pair b={-boundary[i-internal],0};auto d=kpressure_add(pp,b);jump[i]=d.hi+d.lo;
  }
}

KUTTA_KERNEL void fv_correct_velocity_precise_pressure(int i,KIn(float) volume,KIn(float) inverse,KIn(float) diagonal,
    KIn(float) gx,KIn(float) gy,KIn(float) gz,KIn(float) correction,
    KOut(float) u,KOut(float) v,KOut(float) w,KOut(float) p,KOut(float) low,float alpha) {
  float d=volume[i]*inverse[i];u[i]-=d*gx[i];v[i]-=d*gy[i];w[i]-=d*gz[i];
  float a=kmin(alpha,0.7f*diagonal[i]*inverse[i]);kpressure_pair old={p[i],low[i]};
  auto q=kpressure_add(old,kpressure_product(a,correction[i]));p[i]=q.hi;low[i]=q.lo;
}

KUTTA_KERNEL void fv_rhie_chow_precise(int i,KIn(int) owner,KIn(int) neighbour,KIn(int) kind,
    KIn(float) sx,KIn(float) sy,KIn(float) sz,KIn(float) diffusion,KIn(float) volume,
    KIn(float) diagonal,KIn(float) u,KIn(float) v,KIn(float) w,KIn(float) jump,
    KOut(float) flux,KOut(float) low,int internal,float rho) {
  if(i>=internal && (kind[i-internal]==0 || kind[i-internal]==3)) {flux[i]=low[i]=0;return;}
  auto normal=kpressure_add(kpressure_product(u[i],sx[i]),kpressure_product(v[i],sy[i]));
  normal=kpressure_scale(kpressure_add(normal,kpressure_product(w[i],sz[i])),rho);
  if(i<internal || kind[i-internal]!=1) {
    int p=owner[i];float rc=volume[p]/kmax(diagonal[p],1e-20f);
    if(i<internal) {int n=neighbour[i];rc=0.5f*(rc+volume[n]/kmax(diagonal[n],1e-20f));}
    normal=kpressure_add(normal,kpressure_product(rho*rc*diffusion[i],jump[i]));
  }
  flux[i]=normal.hi;low[i]=normal.lo;
}

KUTTA_KERNEL void fv_correct_flux_precise(int i,KIn(int) owner,KIn(int) neighbour,
    KIn(float) correction,KIn(float) coefficient,KOut(float) flux,KOut(float) low,int internal) {
  int p=owner[i];kpressure_pair a={correction[p],0},b={i<internal?-correction[neighbour[i]]:0.0f,0};
  auto delta=kpressure_scale(kpressure_add(a,b),coefficient[i]);kpressure_pair old={flux[i],low[i]};
  auto value=kpressure_add(old,delta);flux[i]=value.hi;low[i]=value.lo;
}

KUTTA_KERNEL void fv_continuity_precise(int i,KIn(int) offsets,KIn(int) faces,KIn(float) sign,
    KIn(float) flux,KIn(float) low,KOut(float) result,KOut(float) rhs,int reference) {
  kpressure_pair sum={0,0};
  for(int j=offsets[i];j<offsets[i+1];++j) {
    int f=faces[j];kpressure_pair value={sign[j]*flux[f],sign[j]*low[f]};sum=kpressure_add(sum,value);
  }
  result[i]=kabs(sum.hi+sum.lo);rhs[i]=i==reference?0.0f:-(sum.hi+sum.lo);
}
