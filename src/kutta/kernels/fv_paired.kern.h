#include "kutta/kernels/prelude.h"

KUTTA_FN kpair fv_pair_product_add(kpair a,float x,float y) {
  float product=x*y;
  kpair sum=ktwo_sum(a.hi,product);
  return ktwo_sum(sum.hi,(a.lo+sum.lo)+kfma(x,y,-product));
}

// Same bounded side reconstruction as fv_pressure_faces, without first
// rounding two O(p) traces and then subtracting away their small difference.
KUTTA_KERNEL void fv_compensated_pressure_jump(int i,KIn(int) owner,KIn(int) neighbour,
    KIn(float) px,KIn(float) py,KIn(float) pz,KIn(float) nx,KIn(float) ny,KIn(float) nz,
    KIn(float) p,KIn(float) gx,KIn(float) gy,KIn(float) gz,KIn(float) boundary,
    KOut(float) jump,int internal) {
  int c=owner[i];
  kpair full=ktwo_sum(p[c],i<internal?-p[neighbour[i]]:-boundary[i-internal]);
  full=fv_pair_product_add(full,gx[c],px[i]);
  full=fv_pair_product_add(full,gy[c],py[i]);
  full=fv_pair_product_add(full,gz[c],pz[i]);
  if(i>=internal) {jump[i]=full.hi+full.lo;return;}
  int n=neighbour[i];
  full=fv_pair_product_add(full,-gx[n],nx[i]);
  full=fv_pair_product_add(full,-gy[n],ny[i]);
  full=fv_pair_product_add(full,-gz[n],nz[i]);
  float side=full.hi+full.lo;
  float dx=px[i]-nx[i],dy=py[i]-ny[i],dz=pz[i]-nz[i];
  float sx=0.5f*(px[i]+nx[i]),sy=0.5f*(py[i]+ny[i]),sz=0.5f*(pz[i]+nz[i]);
  float d2=dx*dx+dy*dy+dz*dz,s2=sx*sx+sy*sy+sz*sz;
  float trust=kmax(0.0f,kmin(1.0f,(4.0f*d2-s2)/(3.0f*d2)));
  if(trust==1.0f) {jump[i]=side;return;}
  kpair shared=ktwo_sum(p[c],-p[n]);
  shared=fv_pair_product_add(shared,gx[c],0.5f*px[i]);shared=fv_pair_product_add(shared,gx[n],0.5f*px[i]);
  shared=fv_pair_product_add(shared,gx[c],-0.5f*nx[i]);shared=fv_pair_product_add(shared,gx[n],-0.5f*nx[i]);
  shared=fv_pair_product_add(shared,gy[c],0.5f*py[i]);shared=fv_pair_product_add(shared,gy[n],0.5f*py[i]);
  shared=fv_pair_product_add(shared,gy[c],-0.5f*ny[i]);shared=fv_pair_product_add(shared,gy[n],-0.5f*ny[i]);
  shared=fv_pair_product_add(shared,gz[c],0.5f*pz[i]);shared=fv_pair_product_add(shared,gz[n],0.5f*pz[i]);
  shared=fv_pair_product_add(shared,gz[c],-0.5f*nz[i]);shared=fv_pair_product_add(shared,gz[n],-0.5f*nz[i]);
  float link=shared.hi+shared.lo;
  jump[i]=link+trust*(side-link);
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
    value=fv_pair_product_add(value,-values[j],hi[cols[j]]);
    float low_coefficient=j==diag[i]?physical[i]:values[j];
    value=fv_pair_product_add(value,-low_coefficient,lo[cols[j]]);
  }
  value=fv_pair_product_add(value,-physical_low[i],hi[i]);
  value=fv_pair_product_add(value,-physical_low[i],lo[i]);
  value=fv_pair_product_add(value,-volume[i],pressure_gradient_lo[i]);
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
    ax=fv_pair_product_add(ax,p,sx[f]);ax=fv_pair_product_add(ax,q,sx[f]);
    ay=fv_pair_product_add(ay,p,sy[f]);ay=fv_pair_product_add(ay,q,sy[f]);
    az=fv_pair_product_add(az,p,sz[f]);az=fv_pair_product_add(az,q,sz[f]);
  }
  ax=fv_pair_product_add(ax,-a/volume[i],volume[i]);
  ay=fv_pair_product_add(ay,-b/volume[i],volume[i]);
  az=fv_pair_product_add(az,-c/volume[i],volume[i]);
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
  if(kabs(offset)<=slack) {flux[i]=previous[i];low[i]=offset;}
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
