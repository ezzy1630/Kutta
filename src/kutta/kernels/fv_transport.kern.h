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
    KIn(float) sign, KIn(int) adjacent, KIn(int) entry, KIn(int) diag,
    KIn(float) diffusion, KIn(float) flux, KIn(float) mu,
    KIn(int) kind, KIn(float) boundary, KIn(float) correction,
    KIn(float) source, KIn(float) volume, KIn(float) old,
    KIn(float) sink,KIn(float) extra_diagonal,KIn(float) extra_rhs,
    KOut(float) values, KOut(float) rhs, KOut(float) diagonal,
    KOut(float) consistent_inverse, int internal, float rho, float dt, float relax) {
  float ap=sink[i]*volume[i]+extra_diagonal[i],offsum=0.0f,b=source[i]*volume[i]+extra_rhs[i];
  for(int j=offsets[i];j<offsets[i+1];++j) if(adjacent[j]>=0) values[entry[j]]=0.0f;
  for(int j=offsets[i];j<offsets[i+1];++j) {
    int f=faces[j],n=adjacent[j];
    float F=sign[j]*flux[f],D=mu[f]*diffusion[f];
    if(n>=0) {
      float a=-D+kmin(F,0.0f);
      values[entry[j]]+=a; offsum+=a;
      ap+=D+kmax(F,0.0f);
    } else {
      int bf=f-internal;
      if(kind[bf]!=0) D=0.0f;
      // Extrapolated outflow has F*phi_P on the implicit diagonal; prescribed
      // backflow is on the RHS and uses the same upwind sign as an inlet.
      ap+=D+kmax(F,0.0f);
      b+=(D-kmin(F,0.0f))*boundary[bf];
    }
    b+=sign[j]*correction[f];
  }
  diagonal[i]=ap;
  float time=dt>0.0f?rho*volume[i]/dt:0.0f;
  float relaxed=(ap+time)/relax;
  b+=(relaxed-ap)*old[i];
  values[diag[i]]=relaxed; rhs[i]=b;
  // SIMPLEC subtracts the neighbour sum from the relaxed diagonal. Pseudo
  // time/under-relaxation keep this positive for convection-dominated cells.
  consistent_inverse[i]=1.0f/kmax(relaxed+offsum,1.0e-20f);
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
