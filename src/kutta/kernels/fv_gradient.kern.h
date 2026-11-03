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

// Apply the same least-squares rows to the physical pair in one reduction.
// Independent float gradients of hi and lo lose the high reduction error,
// making pressure forces depend on arbitrary rebasing of the pressure pair.
KUTTA_KERNEL void fv_paired_gradient(int i,KIn(int) offsets,KIn(int) adjacent,
    KIn(float) dx,KIn(float) dy,KIn(float) dz,
    KIn(float) fx,KIn(float) fy,KIn(float) fz,KIn(float) inverse,
    KIn(float) cell,KIn(float) low,KIn(float) boundary,KIn(float) boundary_low,
    KOut(float) gx,KOut(float) gy,KOut(float) gz,KOut(float) lx,KOut(float) ly,KOut(float) lz,int n) {
  kpair bx={0.0f,0.0f},by={0.0f,0.0f},bz={0.0f,0.0f};
  for(int j=offsets[i];j<offsets[i+1];++j) {
    int other=adjacent[j];
    float hi=other<0?boundary[-other-1]:cell[other];
    float lo=other<0?boundary_low[-other-1]:low[other];
    kpair d=ktwo_sum(hi,-cell[i]);
    d=kpair_product_add(d,lo,1.0f);d=kpair_product_add(d,-low[i],1.0f);
    float wx=other<0?fx[-other-1]:dx[j],wy=other<0?fy[-other-1]:dy[j],wz=other<0?fz[-other-1]:dz[j];
    bx=kpair_product_add(bx,wx,d.hi);bx=kpair_product_add(bx,wx,d.lo);
    by=kpair_product_add(by,wy,d.hi);by=kpair_product_add(by,wy,d.lo);
    bz=kpair_product_add(bz,wz,d.hi);bz=kpair_product_add(bz,wz,d.lo);
  }
  kpair x={0.0f,0.0f},y={0.0f,0.0f},z={0.0f,0.0f};
  x=kpair_product_add(x,inverse[i],bx.hi);x=kpair_product_add(x,inverse[i],bx.lo);
  x=kpair_product_add(x,inverse[i+n],by.hi);x=kpair_product_add(x,inverse[i+n],by.lo);
  x=kpair_product_add(x,inverse[i+2*n],bz.hi);x=kpair_product_add(x,inverse[i+2*n],bz.lo);
  y=kpair_product_add(y,inverse[i+n],bx.hi);y=kpair_product_add(y,inverse[i+n],bx.lo);
  y=kpair_product_add(y,inverse[i+3*n],by.hi);y=kpair_product_add(y,inverse[i+3*n],by.lo);
  y=kpair_product_add(y,inverse[i+4*n],bz.hi);y=kpair_product_add(y,inverse[i+4*n],bz.lo);
  z=kpair_product_add(z,inverse[i+2*n],bx.hi);z=kpair_product_add(z,inverse[i+2*n],bx.lo);
  z=kpair_product_add(z,inverse[i+4*n],by.hi);z=kpair_product_add(z,inverse[i+4*n],by.lo);
  z=kpair_product_add(z,inverse[i+5*n],bz.hi);z=kpair_product_add(z,inverse[i+5*n],bz.lo);
  gx[i]=x.hi;gy[i]=y.hi;gz[i]=z.hi;lx[i]=x.lo;ly[i]=y.lo;lz[i]=z.lo;
}

// Limiter factor that face increment d (= grad . r, nonzero) allows, given
// the room delta (same sign as d) to the neighbour extremum. eps2 < 0 is
// Barth-Jespersen, delta/d. Otherwise Venkatakrishnan's smooth function of
// delta/d (Fluent's differentiable limiter is a modified Venkatakrishnan
// limiter, Theory Guide 23.3.4.3): increments small against sqrt(eps2) pass
// unlimited, so cells near a flat extremum stop switching between 0 and 1.
KUTTA_FN float klimit_face(float d,float delta,float eps2) {
  if(eps2<0.0f) return delta/d;
  float a=delta*delta+eps2,b=d*delta;
  return (a+2.0f*b)/(a+2.0f*d*d+b);
}

// Shifted copies whose max_abs give a field's extremes (operators.cpp
// field_ranges), for any bound >= max|field|.
KUTTA_KERNEL void fv_range_shift(int i, KIn(float) field, KOut(float) above, KOut(float) below, float bound) {
  above[i]=field[i]+bound; below[i]=bound-field[i];
}

KUTTA_KERNEL void fv_limiter(int i, KIn(int) offsets, KIn(int) adjacent,
    KIn(float) rx, KIn(float) ry, KIn(float) rz,
    KIn(float) cell, KIn(float) boundary,
    KIn(float) gx, KIn(float) gy, KIn(float) gz, KOut(float) limiter, float eps2) {
  float lo=cell[i],hi=cell[i];
  for(int j=offsets[i];j<offsets[i+1];++j) {
    int other=adjacent[j];
    float v=other<0 ? boundary[-other-1] : cell[other];
    lo=kmin(lo,v); hi=kmax(hi,v);
  }
  float psi=1.0f;
  for(int j=offsets[i];j<offsets[i+1];++j) {
    float d=gx[i]*rx[j]+gy[i]*ry[j]+gz[i]*rz[j];
    if(d>0.0f) psi=kmin(psi,klimit_face(d,hi-cell[i],eps2));
    if(d<0.0f) psi=kmin(psi,klimit_face(d,lo-cell[i],eps2));
  }
  limiter[i]=kmax(0.0f,kmin(1.0f,psi));
}

// fv_limiter of the float pair cell + low (SolverSettings::paired_state).
// Near a local extremum the limited face value is pinned to the neighbour
// extremum, so a limiter of hi alone jumps whenever hi is rebased, and the
// iteration can settle on a two-cycle across that kink. The neighbour
// differences of the pair are accurate in float; with zero low words this
// is fv_limiter bit for bit. A boundary value of kind 2 is the cell value
// and carries its low word; other boundary values are prescribed.
KUTTA_KERNEL void fv_paired_limiter(int i, KIn(int) offsets, KIn(int) adjacent,
    KIn(float) rx, KIn(float) ry, KIn(float) rz,
    KIn(float) cell, KIn(float) low, KIn(float) boundary, KIn(int) kind,
    KIn(float) gx, KIn(float) gy, KIn(float) gz, KOut(float) limiter, float eps2) {
  float lo=0.0f,hi=0.0f;
  for(int j=offsets[i];j<offsets[i+1];++j) {
    int other=adjacent[j];
    float v=other<0 ? boundary[-other-1] : cell[other];
    float l=other<0 ? (kind[-other-1]==2 ? low[i] : 0.0f) : low[other];
    float d=(v-cell[i])+(l-low[i]);
    lo=kmin(lo,d); hi=kmax(hi,d);
  }
  float psi=1.0f;
  for(int j=offsets[i];j<offsets[i+1];++j) {
    float d=gx[i]*rx[j]+gy[i]*ry[j]+gz[i]*rz[j];
    if(d>0.0f) psi=kmin(psi,klimit_face(d,hi,eps2));
    if(d<0.0f) psi=kmin(psi,klimit_face(d,lo,eps2));
  }
  limiter[i]=kmax(0.0f,kmin(1.0f,psi));
}

// Limit the complete paired gradient against the complete paired extrema.
// Using only grad(hi) misses an extremum when variation is carried by lo.
KUTTA_KERNEL void fv_paired_gradient_limiter(int i,KIn(int) offsets,KIn(int) adjacent,
    KIn(float) rx,KIn(float) ry,KIn(float) rz,
    KIn(float) cell,KIn(float) low,KIn(float) boundary,KIn(float) boundary_low,
    KIn(float) gx,KIn(float) gy,KIn(float) gz,KIn(float) lx,KIn(float) ly,KIn(float) lz,
    KOut(float) limiter,float eps2) {
  float lo=0.0f,hi=0.0f;
  for(int j=offsets[i];j<offsets[i+1];++j) {
    int other=adjacent[j];
    float v=other<0?boundary[-other-1]:cell[other];
    float l=other<0?boundary_low[-other-1]:low[other];
    kpair a={v,l},b={cell[i],low[i]};float d=kpair_delta(a,b);
    lo=kmin(lo,d);hi=kmax(hi,d);
  }
  float psi=1.0f;
  for(int j=offsets[i];j<offsets[i+1];++j) {
    float d=(gx[i]+lx[i])*rx[j]+(gy[i]+ly[i])*ry[j]+(gz[i]+lz[i])*rz[j];
    if(d>0.0f) psi=kmin(psi,klimit_face(d,hi,eps2));
    if(d<0.0f) psi=kmin(psi,klimit_face(d,lo,eps2));
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
    KOut(float) x,KOut(float) y,KOut(float) z,float eps2) {
  fv_limiter(i,offsets,adjacent,rx,ry,rz,cell,boundary,gx,gy,gz,limiter,eps2);
  fv_scale_gradient(i,gx,gy,gz,limiter,x,y,z);
}

KUTTA_KERNEL void fv_paired_limiter_scale(int i,KIn(int) offsets,KIn(int) adjacent,
    KIn(float) rx,KIn(float) ry,KIn(float) rz,
    KIn(float) cell,KIn(float) low,KIn(float) boundary,KIn(int) kind,
    KIn(float) gx,KIn(float) gy,KIn(float) gz,KOut(float) limiter,
    KOut(float) x,KOut(float) y,KOut(float) z,float eps2) {
  fv_paired_limiter(i,offsets,adjacent,rx,ry,rz,cell,low,boundary,kind,gx,gy,gz,limiter,eps2);
  fv_scale_gradient(i,gx,gy,gz,limiter,x,y,z);
}

KUTTA_KERNEL void fv_paired_gradient_limiter_scale(int i,KIn(int) offsets,KIn(int) adjacent,
    KIn(float) rx,KIn(float) ry,KIn(float) rz,
    KIn(float) cell,KIn(float) low,KIn(float) boundary,KIn(float) boundary_low,
    KIn(float) gx,KIn(float) gy,KIn(float) gz,KIn(float) lx,KIn(float) ly,KIn(float) lz,
    KOut(float) limiter,KOut(float) x,KOut(float) y,KOut(float) z,float eps2) {
  fv_paired_gradient_limiter(i,offsets,adjacent,rx,ry,rz,cell,low,boundary,boundary_low,gx,gy,gz,lx,ly,lz,limiter,eps2);
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

// Both directions are evaluated, rather than negating one stored difference:
// subtraction of equal traces must retain the original signed zero on each
// side. The mesh face offsets are the same rounded centroid differences as
// the cell-face offsets consumed by fv_jump_gradient.
KUTTA_KERNEL void fv_face_jump(int i,KIn(int) owner,KIn(int) neighbour,
    KIn(float) px,KIn(float) py,KIn(float) pz,KIn(float) nx,KIn(float) ny,KIn(float) nz,
    KIn(float) cell,KIn(float) boundary,KIn(float) gx,KIn(float) gy,KIn(float) gz,
    KOut(float) forward,KOut(float) reverse,int internal) {
  int p=owner[i];
  float own=cell[p]+gx[p]*px[i]+gy[p]*py[i]+gz[p]*pz[i];
  float next=0.0f;
  if(i<internal) {
    int n=neighbour[i];
    next=cell[n]+gx[n]*nx[i]+gy[n]*ny[i]+gz[n]*nz[i];
    reverse[i]=own-next;
  } else next=boundary[i-internal];
  forward[i]=next-own;
}

KUTTA_KERNEL void fv_gather_jump_gradient(int i,KIn(int) offsets,KIn(int) faces,KIn(float) sign,
    KIn(float) dx,KIn(float) dy,KIn(float) dz,KIn(float) inverse,
    KIn(float) forward,KIn(float) reverse,KIn(float) gx,KIn(float) gy,KIn(float) gz,
    KOut(float) x,KOut(float) y,KOut(float) z,int n) {
  float bx=0.0f,by=0.0f,bz=0.0f;
  for(int j=offsets[i];j<offsets[i+1];++j) {
    int f=faces[j];float jump=sign[j]>0.0f?forward[f]:reverse[f];
    bx+=dx[j]*jump;by+=dy[j]*jump;bz+=dz[j]*jump;
  }
  x[i]=gx[i]+0.5f*(inverse[i]*bx+inverse[i+n]*by+inverse[i+2*n]*bz);
  y[i]=gy[i]+0.5f*(inverse[i+n]*bx+inverse[i+3*n]*by+inverse[i+4*n]*bz);
  z[i]=gz[i]+0.5f*(inverse[i+2*n]*bx+inverse[i+4*n]*by+inverse[i+5*n]*bz);
}

// The face-jump derivative correction of the physical pair. Round only its
// final derivative, rather than rounding the nearly constant cell field first.
KUTTA_KERNEL void fv_paired_jump_gradient(int i,KIn(int) offsets,KIn(int) faces,KIn(int) adjacent,KIn(float) sign,
    KIn(float) dx,KIn(float) dy,KIn(float) dz,KIn(float) inverse,
    KIn(float) px,KIn(float) py,KIn(float) pz,KIn(float) nx,KIn(float) ny,KIn(float) nz,
    KIn(float) cell,KIn(float) low,KIn(float) boundary,KIn(float) boundary_low,
    KIn(float) gx,KIn(float) gy,KIn(float) gz,KIn(float) lx,KIn(float) ly,KIn(float) lz,
    KOut(float) x,KOut(float) y,KOut(float) z,int n) {
  kpair bx={0.0f,0.0f},by={0.0f,0.0f},bz={0.0f,0.0f};
  for(int j=offsets[i];j<offsets[i+1];++j) {
    int f=faces[j],other=adjacent[j];
    float a=sign[j]>0.0f?px[f]:nx[f],b=sign[j]>0.0f?py[f]:ny[f],c=sign[j]>0.0f?pz[f]:nz[f];
    kpair own=kpair_trace(cell[i],low[i],gx[i],gy[i],gz[i],lx[i],ly[i],lz[i],a,b,c);
    kpair next={0.0f,0.0f};
    if(other<0) {next.hi=boundary[-other-1];next.lo=boundary_low[-other-1];}
    else {
      float a=sign[j]>0.0f?nx[f]:px[f],b=sign[j]>0.0f?ny[f]:py[f],c=sign[j]>0.0f?nz[f]:pz[f];
      next=kpair_trace(cell[other],low[other],gx[other],gy[other],gz[other],lx[other],ly[other],lz[other],a,b,c);
    }
    kpair jump=ktwo_sum(next.hi,-own.hi);jump=kpair_product_add(jump,next.lo,1.0f);jump=kpair_product_add(jump,-own.lo,1.0f);
    bx=kpair_product_add(bx,dx[j],jump.hi);bx=kpair_product_add(bx,dx[j],jump.lo);
    by=kpair_product_add(by,dy[j],jump.hi);by=kpair_product_add(by,dy[j],jump.lo);
    bz=kpair_product_add(bz,dz[j],jump.hi);bz=kpair_product_add(bz,dz[j],jump.lo);
  }
  kpair a=ktwo_sum(gx[i],lx[i]),b=ktwo_sum(gy[i],ly[i]),c=ktwo_sum(gz[i],lz[i]);
  a=kpair_product_add(a,.5f*inverse[i],bx.hi);a=kpair_product_add(a,.5f*inverse[i],bx.lo);
  a=kpair_product_add(a,.5f*inverse[i+n],by.hi);a=kpair_product_add(a,.5f*inverse[i+n],by.lo);
  a=kpair_product_add(a,.5f*inverse[i+2*n],bz.hi);a=kpair_product_add(a,.5f*inverse[i+2*n],bz.lo);
  b=kpair_product_add(b,.5f*inverse[i+n],bx.hi);b=kpair_product_add(b,.5f*inverse[i+n],bx.lo);
  b=kpair_product_add(b,.5f*inverse[i+3*n],by.hi);b=kpair_product_add(b,.5f*inverse[i+3*n],by.lo);
  b=kpair_product_add(b,.5f*inverse[i+4*n],bz.hi);b=kpair_product_add(b,.5f*inverse[i+4*n],bz.lo);
  c=kpair_product_add(c,.5f*inverse[i+2*n],bx.hi);c=kpair_product_add(c,.5f*inverse[i+2*n],bx.lo);
  c=kpair_product_add(c,.5f*inverse[i+4*n],by.hi);c=kpair_product_add(c,.5f*inverse[i+4*n],by.lo);
  c=kpair_product_add(c,.5f*inverse[i+5*n],bz.hi);c=kpair_product_add(c,.5f*inverse[i+5*n],bz.lo);
  x[i]=a.hi;y[i]=b.hi;z[i]=c.hi;
}

// Copy stored inverse words; the arithmetic in the gradient kernels is unchanged.
KUTTA_KERNEL void fv_lsq_restore_inverse(int i,KIn(int) cells,KIn(float) values,
    KOut(float) inverse,int count,int n) {
  for(int k=0;k<6;++k) inverse[cells[i]+k*n]=values[i+k*count];
}
