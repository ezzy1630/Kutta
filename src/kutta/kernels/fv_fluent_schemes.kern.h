#include "kutta/kernels/prelude.h"

// Fluent pressure-velocity options outside the default Kutta path. Sources are
// the Fluent Theory Guide sections cited at each kernel.

// SIMPLE (TG "SIMPLE vs. SIMPLEC"): the correction coefficient is the inverse
// of the relaxed momentum diagonal, d = V/a_P, rather than SIMPLEC's
// V/(a_P - sum a_nb). Call immediately after a momentum assembly; the momentum
// components share one CSR value buffer.
KUTTA_KERNEL void fv_simple_inverse(int i,KIn(int) diag,KIn(float) values,KOut(float) inverse) {
  inverse[i]=1.0f/kmax(values[diag[i]],1.0e-20f);
}

// Standard pressure interpolation (TG 23-42): interior faces weight the cell
// pressures by the inverse momentum coefficients,
//   p_f = (p0/a0 + p1/a1)/(1/a0 + 1/a1).
// TG: "Ansys Fluent assumes that the normal pressure gradient at the wall is
// zero", so wall faces (kind 0) take the cell value; symmetry (kind 3) has an
// exactly zero normal gradient. The TG is silent on velocity and mass-flow
// inlets (kind 1). Kutta also uses the cell value there, an inference from
// measurement: on elbow_param it agrees with Fluent's converged fields
// (velocity L2 0.93%), while extrapolating the inlet value diverges
// (validation/results/case_coverage/REPORT.md#standard-pressure).
// Prescribed-pressure faces (kinds 2 and 4) keep their boundary pressure.
// Only the face pressure is replaced: the Rhie-Chow pressure difference is the
// reconstructed one of TG 23-45 for every interpolation scheme.
KUTTA_KERNEL void fv_standard_pressure_faces(int i,KIn(int) owner,KIn(int) neighbour,KIn(int) kind,
    KIn(float) p,KIn(float) diagonal,KOut(float) face,KOut(float) boundary,int internal) {
  int c=owner[i];
  if(i<internal) {
    int n=neighbour[i];
    float wc=1.0f/kmax(diagonal[c],1.0e-20f),wn=1.0f/kmax(diagonal[n],1.0e-20f);
    face[i]=(p[c]*wc+p[n]*wn)/(wc+wn);
    return;
  }
  int b=i-internal;
  if(kind[b]==2 || kind[b]==4) return;
  face[i]=p[c];boundary[b]=p[c];
}

// Unprescribed boundary values seen by the pressure gradient: the cell
// pressure Fluent stores there. With the standard scheme's normal-foot LSQ
// rows this is the zero normal derivative above; second-order pressure
// observes it at the face centroid where its fit keeps such a row, while its
// momentum boundary pressure stays extrapolated. Prescribed faces keep
// fv_pressure_boundary's values.
KUTTA_KERNEL void fv_cell_pressure_boundary(int i,KIn(int) owner,KIn(int) kind,KIn(float) p,
    KOut(float) boundary,int internal) {
  if(kind[i]!=2 && kind[i]!=4) boundary[i]=p[owner[i+internal]];
}

// Wall force per face to the stress on the fluid, matching fv_sst_force's
// convention F = -tau*|A|.
KUTTA_KERNEL void fv_force_to_shear(int i,KIn(float) mask,KIn(float) ax,KIn(float) ay,KIn(float) az,
    KIn(float) fx,KIn(float) fy,KIn(float) fz,KOut(float) sx,KOut(float) sy,KOut(float) sz,int internal) {
  int f=i+internal;float a=ksqrt(ax[f]*ax[f]+ay[f]*ay[f]+az[f]*az[f]);
  float scale=mask[i]>0.0f?-1.0f/a:0.0f;
  sx[i]=scale*fx[i];sy[i]=scale*fy[i];sz[i]=scale*fz[i];
}

// Fluent's poor mesh numerics (User's Guide 36.22.5) give cells of low
// orthogonal quality a locally low-order discretization. Here a flagged
// cell's side of each second-order face pressure (TG 23-43) is its cell
// pressure instead of p + grad(p).r. A first wall-layer cell sheared across
// a sharp convex edge has its centroid tens of cell heights from its face
// centroids, and the extrapolated difference across its two large faces fed
// back through its own pressure force until it diverged
// (validation/results/external_3d_robustness/REPORT.md). Prescribed-pressure
// faces keep their value; the Rhie-Chow jump keeps the full reconstruction.
KUTTA_KERNEL void fv_poor_cell_pressure_faces(int i,KIn(int) owner,KIn(int) neighbour,KIn(int) kind,
    KIn(int) poor,KIn(float) px,KIn(float) py,KIn(float) pz,KIn(float) nx,KIn(float) ny,KIn(float) nz,
    KIn(float) p,KIn(float) gx,KIn(float) gy,KIn(float) gz,KOut(float) face,KOut(float) boundary,int internal) {
  int c=owner[i];
  if(i>=internal) {
    int b=i-internal;
    if(poor[c]==0 || kind[b]==2 || kind[b]==4) return;
    face[i]=p[c];boundary[b]=p[c];
    return;
  }
  int n=neighbour[i];
  if(poor[c]==0 && poor[n]==0) return;
  float pp=poor[c]!=0?p[c]:p[c]+gx[c]*px[i]+gy[c]*py[i]+gz[c]*pz[i];
  float pn=poor[n]!=0?p[n]:p[n]+gx[n]*nx[i]+gy[n]*ny[i]+gz[n]*nz[i];
  face[i]=0.5f*(pp+pn);
}
