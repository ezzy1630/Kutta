#include "kutta/kernels/prelude.h"

KUTTA_KERNEL void fv_projected_momentum_faces(int i,KIn(int) offsets,KIn(int) points,KIn(int) owner,KIn(int) neighbour,KIn(int) mask,
    KIn(float) face_weights,KIn(float) link_weights,KIn(float) field,KIn(float) boundary,KIn(float) diagonal,KIn(float) cell_limiter,
    KOut(float) face,int limit) {
  if(mask[i]==0) return;
  int p=owner[i],n=neighbour[i];float w=diagonal[p]/(diagonal[p]+diagonal[n]);
  float cp=field[p],anchor=0.5f*(cp+field[n]),lo=kmin(cp,field[n]),hi=kmax(cp,field[n]);
  kupper_sum f={0.0f,0.0f},link={0.0f,0.0f};
  for(int j=offsets[i];j<offsets[i+1];++j) {
    int point=points[j];float value=point<0?boundary[-point-1]:field[point];lo=kmin(lo,value);hi=kmax(hi,value);
    f=kupper_add(f,face_weights[j]*(value-anchor));link=kupper_add(link,link_weights[j]*(value-anchor));
  }
  float value=w*cp+(1.0f-w)*field[n]+kupper_value(f)-(1.0f-w)*kupper_value(link);
  if(limit!=0) {
    float theta=kmin(cell_limiter[p],cell_limiter[n]);
    face[i]+=theta*(kmin(hi,kmax(lo,value))-face[i]);
  } else face[i]=value;
}

// Increment of the bounded predictor, with bounds of the complete pair.
// Clipping the low field independently is not its linear response: a low
// extrapolation can exceed its own range while the total stays well inside.
KUTTA_KERNEL void fv_projected_momentum_increment(int i,KIn(int) offsets,KIn(int) points,KIn(int) owner,KIn(int) neighbour,KIn(int) mask,
    KIn(float) face_weights,KIn(float) link_weights,KIn(float) hi,KIn(float) lo,KIn(float) boundary_hi,KIn(float) boundary_lo,
    KIn(float) diagonal,KIn(float) cell_limiter,KOut(float) face) {
  if(mask[i]==0) return;
  int p=owner[i],n=neighbour[i];float w=diagonal[p]/(diagonal[p]+diagonal[n]);
  float ah=0.5f*(hi[p]+hi[n]),al=0.5f*(lo[p]+lo[n]);
  float high_min=kmin(hi[p],hi[n]),high_max=kmax(hi[p],hi[n]);
  kpair minimum={hi[p],lo[p]},maximum=minimum,other={hi[n],lo[n]};
  if(kpair_delta(other,minimum)<0.0f) minimum=other;
  if(kpair_delta(other,maximum)>0.0f) maximum=other;
  kupper_sum fh={0.0f,0.0f},fl={0.0f,0.0f},lh={0.0f,0.0f},ll={0.0f,0.0f};
  for(int j=offsets[i];j<offsets[i+1];++j) {
    int point=points[j];kpair value={point<0?boundary_hi[-point-1]:hi[point],point<0?boundary_lo[-point-1]:lo[point]};
    high_min=kmin(high_min,value.hi);high_max=kmax(high_max,value.hi);
    if(kpair_delta(value,minimum)<0.0f) minimum=value;
    if(kpair_delta(value,maximum)>0.0f) maximum=value;
    fh=kupper_add(fh,face_weights[j]*(value.hi-ah));fl=kupper_add(fl,face_weights[j]*(value.lo-al));
    lh=kupper_add(lh,link_weights[j]*(value.hi-ah));ll=kupper_add(ll,link_weights[j]*(value.lo-al));
  }
  float vh=w*hi[p]+(1.0f-w)*hi[n]+kupper_value(fh)-(1.0f-w)*kupper_value(lh);
  float vl=w*lo[p]+(1.0f-w)*lo[n]+kupper_value(fl)-(1.0f-w)*kupper_value(ll);
  kpair value=ktwo_sum(vh,vl);
  if(kpair_delta(value,minimum)<0.0f) value=minimum;
  if(kpair_delta(value,maximum)>0.0f) value=maximum;
  kpair clipped_high={kmin(high_max,kmax(high_min,vh)),0.0f};
  float delta=kpair_delta(value,clipped_high),theta=kmin(cell_limiter[p],cell_limiter[n]);
  face[i]+=theta*(delta-face[i]);
}

KUTTA_KERNEL void fv_wall_projected_viscous_force(int i,KIn(float) selected,KIn(int) mask,
    KIn(float) sx,KIn(float) sy,KIn(float) sz,
    KIn(float) ux,KIn(float) uy,KIn(float) uz,KIn(float) vx,KIn(float) vy,KIn(float) vz,
    KIn(float) wx,KIn(float) wy,KIn(float) wz,KOut(float) x,KOut(float) y,KOut(float) z,
    int internal,float viscosity) {
  int f=i+internal;if(mask[f]==0) return;float scale=-viscosity*selected[i];
  x[i]=scale*(2.0f*ux[f]*sx[f]+(uy[f]+vx[f])*sy[f]+(uz[f]+wx[f])*sz[f]);
  y[i]=scale*((vx[f]+uy[f])*sx[f]+2.0f*vy[f]*sy[f]+(vz[f]+wy[f])*sz[f]);
  z[i]=scale*((wx[f]+uz[f])*sx[f]+(wy[f]+vz[f])*sy[f]+2.0f*wz[f]*sz[f]);
}

KUTTA_KERNEL void fv_projected_face_gradient(int i,KIn(int) offsets,KIn(int) points,KIn(int) owner,KIn(int) neighbour,
    KIn(float) wx,KIn(float) wy,KIn(float) wz,KIn(int) mask,KIn(float) field,KIn(float) boundary,
    KOut(float) x,KOut(float) y,KOut(float) z,int internal) {
  if(mask[i]==0) return;
  kupper_sum a={0.0f,0.0f},b={0.0f,0.0f},c={0.0f,0.0f};float anchor=i<internal?0.5f*(field[owner[i]]+field[neighbour[i]]):field[owner[i]];
  for(int j=offsets[i];j<offsets[i+1];++j) {
    int p=points[j];float value=(p<0?boundary[-p-1]:field[p])-anchor;
    a=kupper_add(a,wx[j]*value);b=kupper_add(b,wy[j]*value);c=kupper_add(c,wz[j]*value);
  }
  x[i]=kupper_value(a);y[i]=kupper_value(b);z[i]=kupper_value(c);
}

KUTTA_KERNEL void fv_projected_diffusion_correction(int i,KIn(int) offsets,KIn(int) points,KIn(int) owner,KIn(int) neighbour,
    KIn(float) weights,KIn(float) field,KIn(float) boundary,KOut(float) correction,int internal) {
  kupper_sum sum={0.0f,0.0f};float anchor=i<internal?0.5f*(field[owner[i]]+field[neighbour[i]]):field[owner[i]];
  for(int j=offsets[i];j<offsets[i+1];++j) {
    int p=points[j];float value=(p<0?boundary[-p-1]:field[p])-anchor;
    sum=kupper_add(sum,weights[j]*value);
  }
  correction[i]=kupper_value(sum);
}

KUTTA_KERNEL void fv_flux_consistent_gradient(int i,KIn(int) offsets,KIn(int) points,KIn(int) owner,
    KIn(int) neighbour,KIn(int) mask,KIn(float) weights,KIn(float) field,KIn(float) boundary,
    KIn(float) sx,KIn(float) sy,KIn(float) sz,KIn(float) diffusion,
    KOut(float) x,KOut(float) y,KOut(float) z,int internal) {
  if(mask[i]==0) return;
  int p=owner[i];float anchor=i<internal?0.5f*(field[p]+field[neighbour[i]]):field[p];kupper_sum sum={0.0f,0.0f};
  for(int j=offsets[i];j<offsets[i+1];++j) {
    int point=points[j];float delta=(point<0?boundary[-point-1]:field[point])-anchor;
    sum=kupper_add(sum,weights[j]*delta);
  }
  float delta=(i<internal?field[neighbour[i]]:boundary[i-internal])-field[p];
  float flux=diffusion[i]*delta+kupper_value(sum);
  float correction=(flux-x[i]*sx[i]-y[i]*sy[i]-z[i]*sz[i])/(sx[i]*sx[i]+sy[i]*sy[i]+sz[i]*sz[i]);
  x[i]+=correction*sx[i];y[i]+=correction*sy[i];z[i]+=correction*sz[i];
}

KUTTA_KERNEL void fv_projected_transport_faces(int i, KIn(int) owner, KIn(int) neighbour,
    KIn(float) weight, KIn(float) px, KIn(float) py, KIn(float) pz,
    KIn(float) nx, KIn(float) ny, KIn(float) nz,
    KIn(float) tx, KIn(float) ty, KIn(float) tz,
    KIn(float) gx, KIn(float) gy, KIn(float) gz, KIn(float) limiter,
    KIn(float) flux, KIn(float) mu, KIn(int) kind,KIn(int) reconstructed_mask,
    KIn(float) reconstructed_correction,KOut(float) correction,
    int internal, int second_order) {
  int p=owner[i];
  float x=gx[p],y=gy[p],z=gz[p];
  float conv=0.0f;
  if(i<internal) {
    int n=neighbour[i]; float w=weight[i];
    x=w*x+(1.0f-w)*gx[n]; y=w*y+(1.0f-w)*gy[n]; z=w*z+(1.0f-w)*gz[n];
    int up=flux[i]>=0.0f?p:n;
    float rx=flux[i]>=0.0f?px[i]:nx[i],ry=flux[i]>=0.0f?py[i]:ny[i],rz=flux[i]>=0.0f?pz[i]:nz[i];
    conv=limiter[up]*(gx[up]*rx+gy[up]*ry+gz[up]*rz)*flux[i];
  } else {
    if(kind[i-internal]!=0) x=y=z=0.0f;
    if(flux[i]>0.0f && kind[i-internal]!=2) conv=limiter[p]*(gx[p]*px[i]+gy[p]*py[i]+gz[p]*pz[i])*flux[i];
  }
  // Like the legacy kernel, only Dirichlet boundary faces carry diffusion.
  int dirichlet=i<internal || kind[i-internal]==0;
  float diffusion=reconstructed_mask[i]!=0 && dirichlet!=0?reconstructed_correction[i]:x*tx[i]+y*ty[i]+z*tz[i];
  correction[i]=mu[i]*diffusion-(second_order!=0?conv:0.0f);
}

