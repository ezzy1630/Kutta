#include "kutta/kernels/prelude.h"

// Same bounded side reconstruction as fv_pressure_faces, without first
// rounding two O(p) traces and then subtracting away their small difference.
KUTTA_KERNEL void fv_compensated_pressure_jump(int i,KIn(int) owner,KIn(int) neighbour,
    KIn(float) px,KIn(float) py,KIn(float) pz,KIn(float) nx,KIn(float) ny,KIn(float) nz,
    KIn(float) p,KIn(float) gx,KIn(float) gy,KIn(float) gz,KIn(float) boundary,
    KOut(float) jump,int internal) {
  int c=owner[i];
  kpair full=ktwo_sum(p[c],i<internal?-p[neighbour[i]]:-boundary[i-internal]);
  full=kpair_product_add(full,gx[c],px[i]);
  full=kpair_product_add(full,gy[c],py[i]);
  full=kpair_product_add(full,gz[c],pz[i]);
  if(i>=internal) {jump[i]=full.hi+full.lo;return;}
  int n=neighbour[i];
  full=kpair_product_add(full,-gx[n],nx[i]);
  full=kpair_product_add(full,-gy[n],ny[i]);
  full=kpair_product_add(full,-gz[n],nz[i]);
  float side=full.hi+full.lo;
  float dx=px[i]-nx[i],dy=py[i]-ny[i],dz=pz[i]-nz[i];
  float sx=0.5f*(px[i]+nx[i]),sy=0.5f*(py[i]+ny[i]),sz=0.5f*(pz[i]+nz[i]);
  float d2=dx*dx+dy*dy+dz*dz,s2=sx*sx+sy*sy+sz*sz;
  float trust=kmax(0.0f,kmin(1.0f,(4.0f*d2-s2)/(3.0f*d2)));
  if(trust==1.0f) {jump[i]=side;return;}
  kpair shared=ktwo_sum(p[c],-p[n]);
  shared=kpair_product_add(shared,gx[c],0.5f*px[i]);shared=kpair_product_add(shared,gx[n],0.5f*px[i]);
  shared=kpair_product_add(shared,gx[c],-0.5f*nx[i]);shared=kpair_product_add(shared,gx[n],-0.5f*nx[i]);
  shared=kpair_product_add(shared,gy[c],0.5f*py[i]);shared=kpair_product_add(shared,gy[n],0.5f*py[i]);
  shared=kpair_product_add(shared,gy[c],-0.5f*ny[i]);shared=kpair_product_add(shared,gy[n],-0.5f*ny[i]);
  shared=kpair_product_add(shared,gz[c],0.5f*pz[i]);shared=kpair_product_add(shared,gz[n],0.5f*pz[i]);
  shared=kpair_product_add(shared,gz[c],-0.5f*nz[i]);shared=kpair_product_add(shared,gz[n],-0.5f*nz[i]);
  float link=shared.hi+shared.lo;
  jump[i]=link+trust*(side-link);
}

// Keep the face pressure canonical even when the cell pair uses hysteresis.
// A rounded O(p) trace can lose curvature that the compensated Rhie-Chow
// jump retains; carrying only the low cell trace would leave the momentum
// force and mass-flux pressure operators inconsistent on thin wall cells.
KUTTA_FN kpair fv_pressure_trace(float p,float low,float gx,float gy,float gz,
    float lx,float ly,float lz,float x,float y,float z) {
  kpair value=ktwo_sum(p,low);
  value=kpair_product_add(value,gx,x);value=kpair_product_add(value,lx,x);
  value=kpair_product_add(value,gy,y);value=kpair_product_add(value,ly,y);
  value=kpair_product_add(value,gz,z);value=kpair_product_add(value,lz,z);
  return value;
}
KUTTA_FN float fv_pressure_trace_difference(kpair a,kpair b) {
  kpair d=ktwo_sum(a.hi,-b.hi);
  d=kpair_product_add(d,a.lo,1.0f);d=kpair_product_add(d,-b.lo,1.0f);
  return d.hi+d.lo;
}
// Cancel full pressure pairs before rounding a small Rhie-Chow curvature.
KUTTA_KERNEL void fv_paired_pressure_jump(int i,KIn(int) owner,KIn(int) neighbour,
    KIn(float) px,KIn(float) py,KIn(float) pz,KIn(float) nx,KIn(float) ny,KIn(float) nz,
    KIn(float) p,KIn(float) low,KIn(float) gx,KIn(float) gy,KIn(float) gz,
    KIn(float) lx,KIn(float) ly,KIn(float) lz,KIn(float) boundary,KIn(float) boundary_low,
    KOut(float) jump,int internal) {
  int c=owner[i];
  kpair own=fv_pressure_trace(p[c],low[c],gx[c],gy[c],gz[c],lx[c],ly[c],lz[c],px[i],py[i],pz[i]);
  if(i>=internal) {kpair b={boundary[i-internal],boundary_low[i-internal]};jump[i]=fv_pressure_trace_difference(own,b);return;}
  int n=neighbour[i];
  kpair next=fv_pressure_trace(p[n],low[n],gx[n],gy[n],gz[n],lx[n],ly[n],lz[n],nx[i],ny[i],nz[i]);
  float side=fv_pressure_trace_difference(own,next);
  float dx=px[i]-nx[i],dy=py[i]-ny[i],dz=pz[i]-nz[i];
  float sx=0.5f*(px[i]+nx[i]),sy=0.5f*(py[i]+ny[i]),sz=0.5f*(pz[i]+nz[i]);
  float d2=dx*dx+dy*dy+dz*dz,s2=sx*sx+sy*sy+sz*sz;
  float trust=kmax(0.0f,kmin(1.0f,(4.0f*d2-s2)/(3.0f*d2)));
  if(trust==1.0f) {jump[i]=side;return;}
  own=fv_pressure_trace(p[c],low[c],gx[c],gy[c],gz[c],lx[c],ly[c],lz[c],.5f*dx,.5f*dy,.5f*dz);
  next=fv_pressure_trace(p[n],low[n],gx[n],gy[n],gz[n],lx[n],ly[n],lz[n],-.5f*dx,-.5f*dy,-.5f*dz);
  float shared=fv_pressure_trace_difference(own,next);jump[i]=shared+trust*(side-shared);
}

// Retain high-word dot/product errors as well as the low face increment.
// Adding F(hi) and F'(lo) after each has rounded loses a moving mass defect.
KUTTA_KERNEL void fv_compensated_rhie_chow(int i,KIn(int) owner,KIn(int) neighbour,KIn(int) kind,
    KIn(float) sx,KIn(float) sy,KIn(float) sz,KIn(float) diffusion,KIn(float) volume,KIn(float) diagonal,
    KIn(float) u,KIn(float) v,KIn(float) w,KIn(float) lu,KIn(float) lv,KIn(float) lw,
    KIn(float) jump,KOut(float) flux,KOut(float) low,int internal,float rho) {
  if(i>=internal && (kind[i-internal]==0 || kind[i-internal]==3)) {flux[i]=0.0f;low[i]=0.0f;return;}
  int c=owner[i];float rc=volume[c]/kmax(diagonal[c],1.0e-20f);
  if(i<internal) {int n=neighbour[i];rc=.5f*(rc+volume[n]/kmax(diagonal[n],1.0e-20f));}
  kpair value={0.0f,0.0f};
  value=kpair_product_add(value,u[i],sx[i]);value=kpair_product_add(value,lu[i],sx[i]);
  value=kpair_product_add(value,v[i],sy[i]);value=kpair_product_add(value,lv[i],sy[i]);
  value=kpair_product_add(value,w[i],sz[i]);value=kpair_product_add(value,lw[i],sz[i]);
  if(i<internal || kind[i-internal]!=1) value=kpair_product_add(value,rc*diffusion[i],jump[i]);
  kpair scaled={0.0f,0.0f};scaled=kpair_product_add(scaled,rho,value.hi);scaled=kpair_product_add(scaled,rho,value.lo);
  flux[i]=scaled.hi;low[i]=scaled.lo;
}

// Complete the already rounded high boundary trace. Neumann Rhie-Chow
// jumps must cancel the same full owner trace, including its rounding error.
KUTTA_KERNEL void fv_paired_pressure_boundary_low(int i,KIn(int) owner,KIn(int) kind,
    KIn(float) px,KIn(float) py,KIn(float) pz,KIn(float) p,KIn(float) low,
    KIn(float) gx,KIn(float) gy,KIn(float) gz,KIn(float) lx,KIn(float) ly,KIn(float) lz,
    KIn(float) boundary,KOut(float) boundary_low,int internal) {
  if(kind[i]==2 || kind[i]==4) {boundary_low[i]=0.0f;return;}
  int f=i+internal,c=owner[f];
  kpair value=fv_pressure_trace(p[c],low[c],gx[c],gy[c],gz[c],lx[c],ly[c],lz[c],px[f],py[f],pz[f]);
  value=kpair_product_add(value,-boundary[i],1.0f);boundary_low[i]=value.hi+value.lo;
}
KUTTA_KERNEL void fv_paired_pressure_faces(int i,KIn(int) owner,KIn(int) neighbour,KIn(int) kind,KIn(int) poor,
    KIn(float) px,KIn(float) py,KIn(float) pz,KIn(float) nx,KIn(float) ny,KIn(float) nz,
    KIn(float) p,KIn(float) low,KIn(float) gx,KIn(float) gy,KIn(float) gz,KIn(float) lx,KIn(float) ly,KIn(float) lz,
    KIn(float) boundary,KIn(float) boundary_low,KIn(float) diagonal,KOut(float) face,KOut(float) face_low,int internal,int standard) {
  int c=owner[i];kpair value=ktwo_sum(p[c],low[c]);
  if(i<internal) {
    int n=neighbour[i];
    if(standard) {
      float wc=1.0f/kmax(diagonal[c],1.0e-20f),wn=1.0f/kmax(diagonal[n],1.0e-20f);
      kpair numerator={0.0f,0.0f};
      numerator=kpair_product_add(numerator,wc,p[c]);numerator=kpair_product_add(numerator,wc,low[c]);
      numerator=kpair_product_add(numerator,wn,p[n]);numerator=kpair_product_add(numerator,wn,low[n]);
      kpair denominator=ktwo_sum(wc,wn);float q=numerator.hi/denominator.hi;
      numerator=kpair_product_add(numerator,-q,denominator.hi);
      numerator=kpair_product_add(numerator,-q,denominator.lo);
      value=ktwo_sum(q,(numerator.hi+numerator.lo)/denominator.hi);
    } else {
      kpair own=poor[c]?ktwo_sum(p[c],low[c]):fv_pressure_trace(p[c],low[c],gx[c],gy[c],gz[c],lx[c],ly[c],lz[c],px[i],py[i],pz[i]);
      kpair next=poor[n]?ktwo_sum(p[n],low[n]):fv_pressure_trace(p[n],low[n],gx[n],gy[n],gz[n],lx[n],ly[n],lz[n],nx[i],ny[i],nz[i]);
      value=ktwo_sum(0.5f*own.hi,0.5f*next.hi);
      value=ktwo_sum(value.hi,value.lo+0.5f*own.lo+0.5f*next.lo);
    }
  } else {
    int b=i-internal;
    if(kind[b]==2 || kind[b]==4) value=ktwo_sum(boundary[b],boundary_low[b]);
    else if(!standard && !poor[c]) value=fv_pressure_trace(p[c],low[c],gx[c],gy[c],gz[c],lx[c],ly[c],lz[c],px[i],py[i],pz[i]);
  }
  face[i]=value.hi;face_low[i]=value.lo;
}

// SIMPLE/SIMPLEC float-pair state (SolverSettings::paired_state). Velocity and
// pressure are hi + lo pairs (prelude kpair_add); the face mass flux is
// F(hi) + F'(lo), with F' the linearized Rhie-Chow flux of the low words.
// Every existing kernel keeps reading hi alone. These kernels add lo where a
// float32 state otherwise sets the residual floor: the momentum residual,
// the face flux, the continuity sum and the corrections.

// Momentum residual of hi + lo in delta form, compensated: rhs - A (hi + lo)
// - V grad(p_lo). Assembly's relaxation term used old = hi; adding back
// (relaxed - physical) lo makes it relax the full pair. Optional assembly low
// words recover the physical row before positive upward matrix rounding.
KUTTA_KERNEL void fv_paired_momentum_residual(int i,KIn(int) rows,KIn(int) cols,
    KIn(float) values,KIn(int) diag,KIn(float) physical,KIn(float) physical_low,KIn(float) hi,KIn(float) lo,
    KIn(float) pressure_gradient_lo,KIn(float) volume,KIn(float) rhs,KIn(float) rhs_low,KOut(float) r) {
  float sum=rhs[i],error=rhs_low[i];
  for(int j=rows[i];j<rows[i+1];++j) {
    float a=values[j],x=hi[cols[j]],product=a*x;
    kpair s=ktwo_sum(sum,-product);
    sum=s.hi;error+=s.lo-kfma(a,x,-product)-a*lo[cols[j]];
  }
  error+=(values[diag[i]]-physical[i])*lo[i]-volume[i]*pressure_gradient_lo[i];
  error-=physical_low[i]*hi[i]+physical_low[i]*lo[i];
  r[i]=sum+error;
}

// Evaluate the same physical equation without cancelling two large
// pseudo-time products of lo. Carry every high/low product through TwoSum.
KUTTA_KERNEL void fv_compensated_momentum_residual(int i,KIn(int) rows,KIn(int) cols,
    KIn(float) values,KIn(int) diag,KIn(float) physical,KIn(float) physical_low,KIn(float) hi,KIn(float) lo,
    KIn(float) pressure_gradient_lo,KIn(float) volume,KIn(float) rhs,KIn(float) rhs_low,KOut(float) r) {
  kpair value=ktwo_sum(rhs[i],rhs_low[i]);
  for(int j=rows[i];j<rows[i+1];++j) {
    value=kpair_product_add(value,-values[j],hi[cols[j]]);
    float low_coefficient=j==diag[i]?physical[i]:values[j];
    value=kpair_product_add(value,-low_coefficient,lo[cols[j]]);
  }
  value=kpair_product_add(value,-physical_low[i],hi[i]);
  value=kpair_product_add(value,-physical_low[i],lo[i]);
  value=kpair_product_add(value,-volume[i],pressure_gradient_lo[i]);
  r[i]=value.hi+value.lo;
}

KUTTA_KERNEL void fv_paired_accumulate(int i,KIn(float) delta,KOut(float) hi,KOut(float) lo,float slack) {
  kpair x=kpair_add(hi[i],lo[i],delta[i],slack);hi[i]=x.hi;lo[i]=x.lo;
}

// Recover the force lost by the high-word Green-Gauss reduction as well
// as the force of low face pressures. Subtract the same rounded high
// gradient used by fv_momentum_source, so the residual sees their full sum.
KUTTA_KERNEL void fv_compensated_pressure_force(int i,KIn(int) offsets,KIn(int) faces,KIn(float) sign,
    KIn(float) sx,KIn(float) sy,KIn(float) sz,KIn(float) volume,KIn(float) pressure,KIn(float) low,
    KOut(float) x,KOut(float) y,KOut(float) z) {
  float a=0.0f,b=0.0f,c=0.0f;
  kpair ax={0.0f,0.0f},ay={0.0f,0.0f},az={0.0f,0.0f};
  for(int j=offsets[i];j<offsets[i+1];++j) {
    int f=faces[j];float p=sign[j]*pressure[f],q=sign[j]*low[f];
    a+=p*sx[f];b+=p*sy[f];c+=p*sz[f];
    ax=kpair_product_add(ax,p,sx[f]);ax=kpair_product_add(ax,q,sx[f]);
    ay=kpair_product_add(ay,p,sy[f]);ay=kpair_product_add(ay,q,sy[f]);
    az=kpair_product_add(az,p,sz[f]);az=kpair_product_add(az,q,sz[f]);
  }
  ax=kpair_product_add(ax,-a/volume[i],volume[i]);
  ay=kpair_product_add(ay,-b/volume[i],volume[i]);
  az=kpair_product_add(az,-c/volume[i],volume[i]);
  x[i]=(ax.hi+ax.lo)/volume[i];y[i]=(ay.hi+ay.lo)/volume[i];z[i]=(az.hi+az.lo)/volume[i];
}

// Linear boundary observations of the low velocity. Use the high flux for
// inflow/outflow selection and only its low word for normal-backflow speed.
KUTTA_KERNEL void fv_paired_velocity_boundary(int i,KIn(int) owner,KIn(int) kind,
    KIn(float) sx,KIn(float) sy,KIn(float) sz,KIn(float) u,KIn(float) v,KIn(float) w,
    KIn(float) flux,KIn(float) flux_low,KOut(float) bu,KOut(float) bv,KOut(float) bw,
    KOut(int) scalar_kind,int internal,float rho) {
  int f=i+internal,c=owner[f],type=kind[i];
  float x=0.0f,y=0.0f,z=0.0f,a2=sx[f]*sx[f]+sy[f]*sy[f]+sz[f]*sz[f];
  scalar_kind[i]=0;
  if(((type==2 || type==4) && flux[f]>=0.0f) || type==3) {
    x=u[c];y=v[c];z=w[c];scalar_kind[i]=2;
  }
  if(type==4 && flux[f]<0.0f) {
    float factor=flux_low[f]/(rho*a2);x=factor*sx[f];y=factor*sy[f];z=factor*sz[f];
  }
  if(type==3) {
    float normal=(x*sx[f]+y*sy[f]+z*sz[f])/a2;
    x-=normal*sx[f];y-=normal*sy[f];z-=normal*sz[f];
  }
  bu[i]=x;bv[i]=y;bw[i]=z;
}

KUTTA_KERNEL void fv_paired_advection(int i,KIn(int) owner,KIn(int) neighbour,KIn(int) kind,
    KIn(float) px,KIn(float) py,KIn(float) pz,KIn(float) nx,KIn(float) ny,KIn(float) nz,
    KIn(float) field,KIn(float) boundary,KIn(float) gx,KIn(float) gy,KIn(float) gz,KIn(float) limiter,
    KIn(float) flux,KIn(float) flux_low,KOut(float) correction,int internal,int second_order) {
  int c=owner[i],up=c;
  float value=field[c],x=px[i],y=py[i],z=pz[i];
  int reconstruct=second_order;
  if(i<internal) {
    if(flux[i]<0.0f) {up=neighbour[i];x=nx[i];y=ny[i];z=nz[i];}
    value=field[up];
  } else {
    if(flux[i]<0.0f) value=boundary[i-internal];
    if(flux[i]<=0.0f || kind[i-internal]==2) reconstruct=0;
  }
  if(reconstruct) value+=limiter[up]*(gx[up]*x+gy[up]*y+gz[up]*z);
  correction[i]=flux_low[i]*value;
}

KUTTA_KERNEL void fv_paired_face_residual(int i,KIn(int) offsets,KIn(int) faces,KIn(float) sign,
    KIn(float) secondary,KIn(float) advection,KIn(float) stress,KIn(float) boundary_rhs,
    KOut(float) residual,int turbulent) {
  float sum=residual[i],error=boundary_rhs[i];
  for(int j=offsets[i];j<offsets[i+1];++j) {
    int f=faces[j];float value=sign[j]*(secondary[f]-advection[f]+(turbulent?stress[f]:0.0f));
    kpair s=ktwo_sum(sum,value);sum=s.hi;error+=s.lo;
  }
  residual[i]=sum+error;
}

// F'(lo): fv_momentum_faces and fv_rhie_chow applied to the low words. The
// legacy face velocity omits low-word gradient terms and the pressure jump
// uses the shared centroid-link form. Those omissions are O(lo), not an exact
// linearization; compensated_momentum uses the high-word face operators.
//
// flux holds the new F(hi) and returns the pair's hi. Rebuilding hi from
// scratch every iteration would not settle: the convective flux sets the
// momentum diagonal, which sets the face weights and F(hi), and rounding
// leaves that loop on a one-ulp cycle. hi therefore stays at the previous
// flux while the pair is within slack of it.
KUTTA_KERNEL void fv_paired_flux(int i,KIn(int) owner,KIn(int) neighbour,KIn(int) kind,
    KIn(float) sx,KIn(float) sy,KIn(float) sz,
    KIn(float) px,KIn(float) py,KIn(float) pz,KIn(float) nx,KIn(float) ny,KIn(float) nz,
    KIn(float) weights,KIn(float) diagonal,KIn(float) volume,KIn(float) diffusion,
    KIn(float) u,KIn(float) v,KIn(float) w,KIn(float) p,
    KIn(float) gx,KIn(float) gy,KIn(float) gz,KIn(float) previous,
    KOut(float) flux,KOut(float) flux_lo,int internal,float rho,float slack_per_area) {
  int c=owner[i];
  float rc=volume[c]/kmax(diagonal[c],1.0e-20f),low=0.0f;
  if(i<internal) {
    int n=neighbour[i];float weight=weights[c]/(weights[c]+weights[n]);
    float x=weight*u[c]+(1.0f-weight)*u[n];
    float y=weight*v[c]+(1.0f-weight)*v[n];
    float z=weight*w[c]+(1.0f-weight)*w[n];
    rc=0.5f*(rc+volume[n]/kmax(diagonal[n],1.0e-20f));
    float dx=px[i]-nx[i],dy=py[i]-ny[i],dz=pz[i]-nz[i];
    float jump=p[c]-p[n]+0.5f*((gx[c]+gx[n])*dx+(gy[c]+gy[n])*dy+(gz[c]+gz[n])*dz);
    low=rho*(x*sx[i]+y*sy[i]+z*sz[i])+rho*rc*diffusion[i]*jump;
  } else if(kind[i-internal]==2 || kind[i-internal]==4) {
    // Pressure boundaries prescribe hi; their low word is zero.
    float jump=p[c]+gx[c]*px[i]+gy[c]*py[i]+gz[c]*pz[i];
    low=rho*(u[c]*sx[i]+v[c]*sy[i]+w[c]*sz[i])+rho*rc*diffusion[i]*jump;
  }
  float slack=slack_per_area*ksqrt(sx[i]*sx[i]+sy[i]*sy[i]+sz[i]*sz[i]);
  kpair x=ktwo_sum(flux[i],-previous[i]);
  float offset=x.hi+(x.lo+low);
  if(kabs(offset)<=slack) {flux[i]=previous[i];flux_lo[i]=offset;}
  else flux_lo[i]=low;
}

KUTTA_KERNEL void fv_paired_flux_offset(int i,KIn(float) sx,KIn(float) sy,KIn(float) sz,
    KIn(float) previous,KOut(float) flux,KOut(float) low,float slack_per_area) {
  float slack=slack_per_area*ksqrt(sx[i]*sx[i]+sy[i]*sy[i]+sz[i]*sz[i]);
  // Transport coefficients use the high word. Outside the hysteresis band,
  // normalize the total flux so products with low velocity retain their cross term.
  kpair total=ktwo_sum(flux[i],low[i]);
  kpair x=ktwo_sum(total.hi,-previous[i]);
  float offset=x.hi+(x.lo+total.lo);
  int same=(previous[i]>0.0f && total.hi>0.0f) || (previous[i]<0.0f && total.hi<0.0f) || (previous[i]==0.0f && total.hi==0.0f);
  if(same && kabs(offset)<=slack) {flux[i]=previous[i];low[i]=offset;}
  else {flux[i]=total.hi;low[i]=total.lo;}
}

// Signed and absolute net outflow of F(hi) + F'(lo). The hi faces are summed
// with TwoSum, so the result is accurate even though it is a small difference
// of large fluxes.
KUTTA_KERNEL void fv_paired_continuity(int i,KIn(int) offsets,KIn(int) faces,
    KIn(float) sign,KIn(float) flux,KIn(float) flux_lo,KOut(float) imbalance,KOut(float) magnitude) {
  float sum=0.0f,error=0.0f;
  for(int j=offsets[i];j<offsets[i+1];++j) {
    int f=faces[j];kpair s=ktwo_sum(sum,sign[j]*flux[f]);
    sum=s.hi;error+=s.lo+sign[j]*flux_lo[f];
  }
  float mass=sum+error;
  imbalance[i]=mass;magnitude[i]=kabs(mass);
}

// fv_pressure_assemble's right-hand side, from the compensated imbalance.
KUTTA_KERNEL void fv_paired_pressure_rhs(int i,KIn(float) imbalance,KOut(float) rhs,int reference) {
  rhs[i]=i==reference?0.0f:-imbalance[i];
}

// slack_per_area: the flux slack per unit face area.
KUTTA_KERNEL void fv_paired_correct_flux(int i,KIn(int) owner,KIn(int) neighbour,
    KIn(float) sx,KIn(float) sy,KIn(float) sz,KIn(float) correction,KIn(float) coefficient,
    KOut(float) flux,KOut(float) flux_lo,int internal,float slack_per_area) {
  int c=owner[i];float pc=correction[c];
  float pn=i<internal?correction[neighbour[i]]:0.0f;
  float slack=slack_per_area*ksqrt(sx[i]*sx[i]+sy[i]*sy[i]+sz[i]*sz[i]);
  kpair x=kpair_add(flux[i],flux_lo[i],coefficient[i]*(pc-pn),slack);flux[i]=x.hi;flux_lo[i]=x.lo;
}

// fv_correct_velocity on the pairs.
KUTTA_KERNEL void fv_paired_correct_velocity(int i,KIn(float) volume,KIn(float) inverse,KIn(float) diagonal,
    KIn(float) gx,KIn(float) gy,KIn(float) gz,KIn(float) correction,
    KOut(float) u,KOut(float) v,KOut(float) w,KOut(float) p,
    KOut(float) u_lo,KOut(float) v_lo,KOut(float) w_lo,KOut(float) p_lo,float alpha,
    float velocity_slack,float pressure_slack) {
  float d=volume[i]*inverse[i];
  kpair x=kpair_add(u[i],u_lo[i],-d*gx[i],velocity_slack);u[i]=x.hi;u_lo[i]=x.lo;
  x=kpair_add(v[i],v_lo[i],-d*gy[i],velocity_slack);v[i]=x.hi;v_lo[i]=x.lo;
  x=kpair_add(w[i],w_lo[i],-d*gz[i],velocity_slack);w[i]=x.hi;w_lo[i]=x.lo;
  float pressure_alpha=kmin(alpha,0.7f*diagonal[i]*inverse[i]);
  x=kpair_add(p[i],p_lo[i],pressure_alpha*correction[i],pressure_slack);p[i]=x.hi;p_lo[i]=x.lo;
}

KUTTA_KERNEL void fv_paired_value(int i,KIn(float) hi,KIn(float) lo,KOut(float) value) {
  value[i]=hi[i]+lo[i];
}
KUTTA_KERNEL void fv_paired_vector_values(int i,KIn(float) x,KIn(float) y,KIn(float) z,
    KIn(float) lx,KIn(float) ly,KIn(float) lz,KOut(float) u,KOut(float) v,KOut(float) w) {
  u[i]=x[i]+lx[i];v[i]=y[i]+ly[i];w[i]=z[i]+lz[i];
}
KUTTA_KERNEL void fv_compensated_correct_flux(int i,KIn(int) owner,KIn(int) neighbour,
    KIn(float) sx,KIn(float) sy,KIn(float) sz,KIn(float) correction,KIn(float) coefficient,
    KOut(float) flux,KOut(float) flux_lo,int internal,float slack_per_area) {
  int c=owner[i];float pn=i<internal?correction[neighbour[i]]:0.0f;
  float slack=slack_per_area*ksqrt(sx[i]*sx[i]+sy[i]*sy[i]+sz[i]*sz[i]);
  kpair x=kpair_add_flux(flux[i],flux_lo[i],coefficient[i]*(correction[c]-pn),slack);
  flux[i]=x.hi;flux_lo[i]=x.lo;
}

KUTTA_KERNEL void fv_compensated_correct_velocity(int i,KIn(float) volume,KIn(float) inverse,KIn(float) diagonal,
    KIn(float) gx,KIn(float) gy,KIn(float) gz,KIn(float) correction,
    KOut(float) u,KOut(float) v,KOut(float) w,KOut(float) p,
    KOut(float) u_lo,KOut(float) v_lo,KOut(float) w_lo,KOut(float) p_lo,float alpha,
    float velocity_slack,float pressure_slack) {
  float d=volume[i]*inverse[i];
  kpair x=kpair_add(u[i],u_lo[i],-d*gx[i],velocity_slack);u[i]=x.hi;u_lo[i]=x.lo;
  x=kpair_add(v[i],v_lo[i],-d*gy[i],velocity_slack);v[i]=x.hi;v_lo[i]=x.lo;
  x=kpair_add(w[i],w_lo[i],-d*gz[i],velocity_slack);w[i]=x.hi;w_lo[i]=x.lo;
  float pressure_alpha=kmin(alpha,0.7f*diagonal[i]*inverse[i]);
  x=kpair_add_relative(p[i],p_lo[i],pressure_alpha*correction[i],pressure_slack);p[i]=x.hi;p_lo[i]=x.lo;
}
