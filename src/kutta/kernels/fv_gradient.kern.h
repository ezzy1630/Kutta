#include "kutta/kernels/prelude.h"

// Row weights: dx/dy/dz per cell-face entry for internal faces, and
// fx/fy/fz per boundary face, so fits that differ only in boundary rows share
// the internal ones.
KUTTA_KERNEL void fv_gradient(int i, KIn(int) offsets, KIn(int) adjacent,
    KIn(float) dx, KIn(float) dy, KIn(float) dz,
    KIn(float) fx, KIn(float) fy, KIn(float) fz, KIn(float) inverse,
    KIn(float) cell, KIn(float) boundary,
    KOut(float) gx, KOut(float) gy, KOut(float) gz, int n) {
  float bx=0.0f,by=0.0f,bz=0.0f;
  for(int j=offsets[i];j<offsets[i+1];++j) {
    int other=adjacent[j];
    float d=(other<0 ? boundary[-other-1] : cell[other])-cell[i];
    float wx=other<0 ? fx[-other-1] : dx[j];
    float wy=other<0 ? fy[-other-1] : dy[j];
    float wz=other<0 ? fz[-other-1] : dz[j];
    bx+=wx*d; by+=wy*d; bz+=wz*d;
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

KUTTA_KERNEL void fv_vector_gradient(int i,KIn(int) offsets,KIn(int) adjacent,
    KIn(float) dx,KIn(float) dy,KIn(float) dz,KIn(float) fx,KIn(float) fy,KIn(float) fz,KIn(float) inverse,
    KIn(float) u,KIn(float) v,KIn(float) w,KIn(float) bu,KIn(float) bv,KIn(float) bw,
    KOut(float) ux,KOut(float) uy,KOut(float) uz,KOut(float) vx,KOut(float) vy,KOut(float) vz,
    KOut(float) wx,KOut(float) wy,KOut(float) wz,int n) {
  float ax=0.0f,ay=0.0f,az=0.0f,bx=0.0f,by=0.0f,bz=0.0f,cx=0.0f,cy=0.0f,cz=0.0f;
  for(int j=offsets[i];j<offsets[i+1];++j) {
    int other=adjacent[j];
    float a=(other<0?bu[-other-1]:u[other])-u[i];
    float b=(other<0?bv[-other-1]:v[other])-v[i];
    float c=(other<0?bw[-other-1]:w[other])-w[i];
    float rx=other<0?fx[-other-1]:dx[j],ry=other<0?fy[-other-1]:dy[j],rz=other<0?fz[-other-1]:dz[j];
    ax+=rx*a;ay+=ry*a;az+=rz*a;
    bx+=rx*b;by+=ry*b;bz+=rz*b;
    cx+=rx*c;cy+=ry*c;cz+=rz*c;
  }
  float xx=inverse[i],xy=inverse[i+n],xz=inverse[i+2*n],yy=inverse[i+3*n],yz=inverse[i+4*n],zz=inverse[i+5*n];
  ux[i]=xx*ax+xy*ay+xz*az;uy[i]=xy*ax+yy*ay+yz*az;uz[i]=xz*ax+yz*ay+zz*az;
  vx[i]=xx*bx+xy*by+xz*bz;vy[i]=xy*bx+yy*by+yz*bz;vz[i]=xz*bx+yz*by+zz*bz;
  wx[i]=xx*cx+xy*cy+xz*cz;wy[i]=xy*cx+yy*cy+yz*cz;wz[i]=xz*cx+yz*cy+zz*cz;
}

KUTTA_KERNEL void fv_scale_gradient(int i,KIn(float) gx,KIn(float) gy,KIn(float) gz,
    KIn(float) limiter,KOut(float) x,KOut(float) y,KOut(float) z) {
  x[i]=limiter[i]*gx[i];y[i]=limiter[i]*gy[i];z[i]=limiter[i]*gz[i];
}

KUTTA_KERNEL void fv_limiter_scale(int i,KIn(int) offsets,KIn(int) adjacent,
    KIn(float) rx,KIn(float) ry,KIn(float) rz,
    KIn(float) cell,KIn(float) boundary,
    KIn(float) gx,KIn(float) gy,KIn(float) gz,KOut(float) limiter,
    KOut(float) x,KOut(float) y,KOut(float) z) {
  fv_limiter(i,offsets,adjacent,rx,ry,rz,cell,boundary,gx,gy,gz,limiter);
  fv_scale_gradient(i,gx,gy,gz,limiter,x,y,z);
}

// Single pass only: iterating this correction feeds jumps back into themselves.
KUTTA_KERNEL void fv_jump_gradient(int i,KIn(int) offsets,KIn(int) faces,KIn(int) adjacent,KIn(float) sign,
    KIn(float) dx,KIn(float) dy,KIn(float) dz,KIn(float) inverse,
    KIn(float) rx,KIn(float) ry,KIn(float) rz,
    KIn(float) px,KIn(float) py,KIn(float) pz,KIn(float) nx,KIn(float) ny,KIn(float) nz,
    KIn(float) cell,KIn(float) boundary,KIn(float) gx,KIn(float) gy,KIn(float) gz,
    KOut(float) x,KOut(float) y,KOut(float) z,int n) {
  float bx=0.0f,by=0.0f,bz=0.0f;
  for(int j=offsets[i];j<offsets[i+1];++j) {
    int f=faces[j],other=adjacent[j];
    float own=cell[i]+gx[i]*rx[j]+gy[i]*ry[j]+gz[i]*rz[j];
    float next=0.0f;
    if(other<0) next=boundary[-other-1];
    else {
      float a=sign[j]>0.0f?nx[f]:px[f],b=sign[j]>0.0f?ny[f]:py[f],c=sign[j]>0.0f?nz[f]:pz[f];
      next=cell[other]+gx[other]*a+gy[other]*b+gz[other]*c;
    }
    float jump=next-own;bx+=dx[j]*jump;by+=dy[j]*jump;bz+=dz[j]*jump;
  }
  x[i]=gx[i]+0.5f*(inverse[i]*bx+inverse[i+n]*by+inverse[i+2*n]*bz);
  y[i]=gy[i]+0.5f*(inverse[i+n]*bx+inverse[i+3*n]*by+inverse[i+4*n]*bz);
  z[i]=gz[i]+0.5f*(inverse[i+2*n]*bx+inverse[i+4*n]*by+inverse[i+5*n]*bz);
}
