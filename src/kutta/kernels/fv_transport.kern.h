#include "kutta/kernels/prelude.h"

KUTTA_KERNEL void fv_transport_faces(int i, KIn(int) owner, KIn(int) neighbour,
    KIn(float) weight, KIn(float) px, KIn(float) py, KIn(float) pz,
    KIn(float) nx, KIn(float) ny, KIn(float) nz,
    KIn(float) tx, KIn(float) ty, KIn(float) tz,
    KIn(float) gx, KIn(float) gy, KIn(float) gz, KIn(float) limiter,
    KIn(float) flux, KIn(float) mu, KIn(int) kind, KOut(float) correction,
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
  correction[i]=mu[i]*(x*tx[i]+y*ty[i]+z*tz[i])-(second_order!=0?conv:0.0f);
}

KUTTA_KERNEL void fv_transport_assemble(int i, KIn(int) offsets, KIn(int) faces,
    KIn(float) sign, KIn(int) adjacent, KIn(int) entry, KIn(int) rows, KIn(int) diag,
    KIn(float) diffusion, KIn(float) flux, KIn(float) mu,
    KIn(int) kind, KIn(float) boundary, KIn(float) correction,
    KIn(float) source, KIn(float) volume, KIn(float) old,
    KIn(float) sink,KIn(float) extra_diagonal,KIn(float) extra_rhs,
    KOut(float) values, KOut(float) rhs,KOut(float) rhs_low, KOut(float) diagonal,KOut(float) diagonal_low,
    KOut(float) consistent_inverse, int internal, float rho, float dt, float relax,int compensated) {
  float sink_product=sink[i]*volume[i];
  kupper_sum ap={sink_product,0.0f};ap=kupper_add(ap,extra_diagonal[i]);
  float b=source[i]*volume[i]+extra_rhs[i];
  float error=0.0f;
  if(compensated) {
    float product=source[i]*volume[i];
    // Keep the high word consistent with its carry; a fused initial sum
    // would already include the multiplication error retained below.
    kpair initial=ktwo_sum(product,extra_rhs[i]);
    b=initial.hi;error=initial.lo+kfma(source[i],volume[i],-product);
  }
  for(int j=offsets[i];j<offsets[i+1];++j) if(adjacent[j]>=0) values[entry[j]]=0.0f;
  for(int j=offsets[i];j<offsets[i+1];++j) {
    int f=faces[j],n=adjacent[j];
    float F=sign[j]*flux[f],D=mu[f]*diffusion[f];
    // aP = sum(|stored aN|) + net outward mass flux + boundary anchors/sinks.
    // A negative net flux is a real transport term: do not clamp it to force
    // diagonal dominance during a non-solenoidal SIMPLE transient.
    ap=kupper_add(ap,F);
    if(n>=0) {
      float a=-D+kmin(F,0.0f);
      values[entry[j]]+=a;
    } else {
      int bf=f-internal;
      if(kind[bf]!=0) D=0.0f;
      // Extrapolated outflow has F*phi_P on the implicit diagonal; prescribed
      // backflow is on the RHS and uses the same upwind sign as an inlet.
      ap=kupper_add(ap,D);ap=kupper_add(ap,-kmin(F,0.0f));
      if(compensated) {
        kpair factor=ktwo_sum(D,-kmin(F,0.0f));
        float product=factor.hi*boundary[bf];
        kpair sum=ktwo_sum(b,product);
        error+=sum.lo+kfma(factor.hi,boundary[bf],-product)+factor.lo*boundary[bf];
        b=sum.hi;
      } else b+=(D-kmin(F,0.0f))*boundary[bf];
    }
    float secondary=sign[j]*correction[f];
    if(compensated) {kpair sum=ktwo_sum(b,secondary);error+=sum.lo;b=sum.hi;}
    else b+=secondary;
  }
  for(int j=rows[i];j<rows[i+1];++j) if(j!=diag[i]) ap=kupper_add(ap,-values[j]);
  float physical=kupper_value(ap);
  diagonal[i]=physical;
  // The upward-rounded matrix remains a positive linear operator. Its small
  // row surplus is not physical reaction: retain it for the paired residual.
  if(compensated) {
    kpair difference=ktwo_sum(ap.sum,-physical);
    diagonal_low[i]=difference.hi+(difference.lo+ap.error+kfma(sink[i],volume[i],-sink_product));
  }
  float time=dt>0.0f?rho*volume[i]/dt:0.0f;
  float relaxed=kadd_up(physical,time)/relax;
  if(compensated) {
    kpair factor=ktwo_sum(relaxed,-physical);
    float product=factor.hi*old[i];
    kpair sum=ktwo_sum(b,product);
    error+=sum.lo+kfma(factor.hi,old[i],-product)+factor.lo*old[i];
    b=sum.hi;
  } else b+=(relaxed-physical)*old[i];
  values[diag[i]]=relaxed; rhs[i]=b;
  if(compensated) rhs_low[i]=error;
  // A non-solenoidal startup flux can make the SIMPLEC denominator negative
  // despite relaxation/pseudo time. Use SIMPLE's positive diagonal then;
  // clamping the negative value to 1e-20 creates a spurious 1e20 inverse.
  kupper_sum denominator={relaxed,0.0f};
  for(int j=rows[i];j<rows[i+1];++j) if(j!=diag[i]) denominator=kupper_add(denominator,values[j]);
  float consistent=kupper_value(denominator);
  consistent_inverse[i]=1.0f/kmax(consistent>0.0f?consistent:relaxed,1.0e-20f);
}

// Each matrix is rebuilt on the next assembly; no persistent constraints.
KUTTA_KERNEL void fv_fixed_rows(int i,KIn(int) rows,KIn(int) diag,KIn(int) mask,
    KIn(float) fixed,KOut(float) values,KOut(float) rhs) {
  if(mask[i]==0) return;
  // Preserve the transport row's units. An unscaled identity row for omega
  // can outweigh all interior residuals by orders of magnitude, causing an
  // inexact Krylov solve to stop before resolving interior transport.
  float scale=values[diag[i]];
  for(int j=rows[i];j<rows[i+1];++j) values[j]=0.0f;
  values[diag[i]]=scale;rhs[i]=scale*fixed[i];
}

KUTTA_KERNEL void fv_floor(int i,KOut(float) field,float floor) {
  field[i]=kmax(field[i],floor);
}

KUTTA_KERNEL void fv_scalar_monitor(int i,KIn(float) field,KIn(float) diagonal,
    KIn(float) residual,KIn(int) excluded,KIn(float) volume,KOut(float) absolute,KOut(float) scale,
    KOut(float) fluent_scale,float rho,float dt) {
  absolute[i]=excluded[i]!=0?0.0f:kabs(residual[i]);
  scale[i]=excluded[i]!=0?0.0f:kabs(diagonal[i]*field[i]);
  float time=dt>0.0f?rho*volume[i]/dt:0.0f;
  fluent_scale[i]=excluded[i]!=0?0.0f:kabs((diagonal[i]+time)*field[i]);
}
