#include "kutta/kernels/prelude.h"

// Boundary kind: wall=0, inlet=1, prescribed-backflow outlet=2, symmetry=3,
// pressure outlet with normal backflow=4.
KUTTA_KERNEL void fv_velocity_boundary(int i,KIn(int) owner,KIn(int) kind,
    KIn(float) sx,KIn(float) sy,KIn(float) sz,
    KIn(float) u,KIn(float) v,KIn(float) w,
    KIn(float) prescribed_u,KIn(float) prescribed_v,KIn(float) prescribed_w,
    KIn(float) flux,KOut(float) bu,KOut(float) bv,KOut(float) bw,
    KOut(int) scalar_kind,int internal,float rho) {
  int f=i+internal,p=owner[f],type=kind[i];
  float x=prescribed_u[i],y=prescribed_v[i],z=prescribed_w[i];
  float a2=sx[f]*sx[f]+sy[f]*sy[f]+sz[f]*sz[f];
  scalar_kind[i]=0;
  if((type==2 || type==4) && flux[f]>=0.0f) {
    x=u[p];y=v[p];z=w[p];
    scalar_kind[i]=2; // cell-value outflow: no diffusion or SOU correction
  }
  if(type==4 && flux[f]<0.0f) {
    float area=ksqrt(a2),speed=flux[f]/(rho*area);
    x=speed*(sx[f]/area);y=speed*(sy[f]/area);z=speed*(sz[f]/area);
  }
  if(type==3) {
    x=u[p];y=v[p];z=w[p];
    scalar_kind[i]=2; // symmetry has no net viscous flux
  }
  if(type==0 || type==3) {
    float normal=(x*sx[f]+y*sy[f]+z*sz[f])/a2;
    x-=normal*sx[f];y-=normal*sy[f];z-=normal*sz[f];
  }
  bu[i]=x;bv[i]=y;bw[i]=z;
}

KUTTA_KERNEL void fv_pressure_boundary(int i,KIn(int) owner,KIn(int) kind,
    KIn(float) rx,KIn(float) ry,KIn(float) rz,KIn(float) p,
    KIn(float) gx,KIn(float) gy,KIn(float) gz,KIn(float) prescribed,
    KOut(float) boundary,int internal) {
  int f=i+internal,c=owner[f];
  float x=rx[f],y=ry[f],z=rz[f];
  boundary[i]=(kind[i]==2 || kind[i]==4)?prescribed[i]:p[c]+gx[c]*x+gy[c]*y+gz[c]*z;
}

KUTTA_KERNEL void fv_pressure_faces(int i,KIn(int) owner,KIn(int) neighbour,
    KIn(float) px,KIn(float) py,KIn(float) pz,
    KIn(float) nx,KIn(float) ny,KIn(float) nz,KIn(float) p,
    KIn(float) gx,KIn(float) gy,KIn(float) gz,KIn(float) boundary,
    KOut(float) face,KOut(float) jump,int internal) {
  int c=owner[i];
  float pp=p[c]+gx[c]*px[i]+gy[c]*py[i]+gz[c]*pz[i];
  if(i<internal) {
    int n=neighbour[i];
    float pn=p[n]+gx[n]*nx[i]+gy[n]*ny[i]+gz[n]*nz[i];
    face[i]=0.5f*(pp+pn);jump[i]=pp-pn;
  } else {face[i]=boundary[i-internal];jump[i]=pp-boundary[i-internal];}
}

KUTTA_KERNEL void fv_momentum_source(int i,KIn(int) offsets,KIn(int) faces,
    KIn(float) sign,KIn(float) sx,KIn(float) sy,KIn(float) sz,
    KIn(float) face_pressure,KIn(float) volume,
    KIn(float) body_u,KIn(float) body_v,KIn(float) body_w,
    KOut(float) su,KOut(float) sv,KOut(float) sw) {
  float x=0.0f,y=0.0f,z=0.0f;
  for(int j=offsets[i];j<offsets[i+1];++j) {
    int f=faces[j];float p=sign[j]*face_pressure[f];
    x+=p*sx[f];y+=p*sy[f];z+=p*sz[f];
  }
  su[i]=body_u[i]-x/volume[i];sv[i]=body_v[i]-y/volume[i];sw[i]=body_w[i]-z/volume[i];
}

// Correct the momentum-weighted interpolation point to the actual face.
// Sharing the gradient preserves affine fields with unequal diagonals and
// skewness, without the side mean's wall-layer curvature term (VERIFICATION.md).
KUTTA_KERNEL void fv_momentum_faces(int i,KIn(int) owner,KIn(int) neighbour,KIn(int) kind,
    KIn(float) px,KIn(float) py,KIn(float) pz,KIn(float) nx,KIn(float) ny,KIn(float) nz,
    KIn(float) u,KIn(float) v,KIn(float) w,
    KIn(float) ux,KIn(float) uy,KIn(float) uz,KIn(float) vx,KIn(float) vy,KIn(float) vz,
    KIn(float) wx,KIn(float) wy,KIn(float) wz,KIn(float) diagonal,
    KIn(float) bu,KIn(float) bv,KIn(float) bw,
    KOut(float) fu,KOut(float) fv,KOut(float) fw,int internal) {
  if(i>=internal) {
    int b=i-internal,p=owner[i];
    if(kind[b]==2 || kind[b]==4) {
      // Pressure fixes outlet flux through momentum, not through a prescribed
      // backflow velocity. Backflow data belongs to advected momentum only.
      fu[i]=u[p];fv[i]=v[p];fw[i]=w[p];
    } else {fu[i]=bu[b];fv[i]=bv[b];fw[i]=bw[b];}
    return;
  }
  int p=owner[i],n=neighbour[i];float weight=diagonal[p]/(diagonal[p]+diagonal[n]);
  float rx=weight*px[i]+(1.0f-weight)*nx[i];
  float ry=weight*py[i]+(1.0f-weight)*ny[i];
  float rz=weight*pz[i]+(1.0f-weight)*nz[i];
  fu[i]=weight*u[p]+(1.0f-weight)*u[n]
    +0.5f*((ux[p]+ux[n])*rx+(uy[p]+uy[n])*ry+(uz[p]+uz[n])*rz);
  fv[i]=weight*v[p]+(1.0f-weight)*v[n]
    +0.5f*((vx[p]+vx[n])*rx+(vy[p]+vy[n])*ry+(vz[p]+vz[n])*rz);
  fw[i]=weight*w[p]+(1.0f-weight)*w[n]
    +0.5f*((wx[p]+wx[n])*rx+(wy[p]+wy[n])*ry+(wz[p]+wz[n])*rz);
}

// The physical (unrelaxed, no pseudo-time) momentum diagonal is used in the
// Rhie–Chow term. Thus its converged pressure smoothing is independent of
// relaxation and pseudo time (Majumdar consistency). SIMPLEC uses its own
// consistent inverse for the pressure-correction operator.
KUTTA_KERNEL void fv_rhie_chow(int i,KIn(int) owner,KIn(int) neighbour,
    KIn(int) kind,KIn(float) sx,KIn(float) sy,KIn(float) sz,
    KIn(float) diffusion,KIn(float) volume,KIn(float) diagonal,
    KIn(float) inverse,KIn(float) u,KIn(float) v,KIn(float) w,
    KIn(float) pressure_jump,
    KOut(float) flux,KOut(float) coefficient,int internal,float rho) {
  int c=owner[i];
  float a=diagonal[c],x=u[i],y=v[i],z=w[i];
  float rc=volume[c]/kmax(a,1.0e-20f),dc=volume[c]*inverse[c];
  float jump=pressure_jump[i];
  if(i<internal) {
    int n=neighbour[i];float an=diagonal[n];
    rc=0.5f*(rc+volume[n]/kmax(an,1.0e-20f));
    dc=0.5f*(dc+volume[n]*inverse[n]);
  } else {
    int b=i-internal;
    if(kind[b]==0 || kind[b]==3) {flux[i]=0.0f;coefficient[i]=0.0f;return;}
    if(kind[b]==1) {
      flux[i]=rho*(x*sx[i]+y*sy[i]+z*sz[i]);coefficient[i]=0.0f;return;
    }

  }
  flux[i]=rho*(x*sx[i]+y*sy[i]+z*sz[i])+rho*rc*diffusion[i]*jump;
  coefficient[i]=rho*dc*diffusion[i];
}

KUTTA_KERNEL void fv_pressure_assemble(int i,KIn(int) offsets,KIn(int) faces,
    KIn(float) sign,KIn(int) adjacent,KIn(int) entry,KIn(int) diag,
    KIn(float) coefficient,KIn(float) flux,KOut(float) values,KOut(float) rhs,
    int reference) {
  float ap=0.0f,b=0.0f;
  for(int j=offsets[i];j<offsets[i+1];++j) if(adjacent[j]>=0) values[entry[j]]=0.0f;
  for(int j=offsets[i];j<offsets[i+1];++j) {
    int f=faces[j],n=adjacent[j];float d=coefficient[f];
    b-=sign[j]*flux[f];
    ap+=d;
    // Eliminate the reference column as well as its row, keeping SPD for CG.
    if(n>=0 && n!=reference && i!=reference) values[entry[j]]-=d;
  }
  values[diag[i]]=i==reference?1.0f:ap;rhs[i]=i==reference?0.0f:b;
}

KUTTA_KERNEL void fv_correct_flux(int i,KIn(int) owner,KIn(int) neighbour,
    KIn(float) correction,KIn(float) coefficient,KOut(float) flux,int internal) {
  int c=owner[i];float pc=correction[c];
  float pn=i<internal?correction[neighbour[i]]:0.0f;
  flux[i]+=coefficient[i]*(pc-pn);
}

KUTTA_KERNEL void fv_correct_velocity(int i,KIn(float) volume,KIn(float) inverse,KIn(float) diagonal,
    KIn(float) gx,KIn(float) gy,KIn(float) gz,KIn(float) correction,
    KOut(float) u,KOut(float) v,KOut(float) w,KOut(float) p,float alpha) {
  float d=volume[i]*inverse[i];
  u[i]-=d*gx[i];v[i]-=d*gy[i];w[i]-=d*gz[i];
  // Majumdar's physical coefficient can exceed the SIMPLEC correction
  // coefficient when pseudo time dominates. Bound the pressure update by
  // their local ratio to avoid overshooting that pressure-smoothing mode.
  float pressure_alpha=kmin(alpha,0.7f*diagonal[i]*inverse[i]);
  p[i]+=pressure_alpha*correction[i];
}

KUTTA_KERNEL void fv_monitor_cell(int i,KIn(int) offsets,KIn(int) faces,
    KIn(float) sign,KIn(float) flux,KIn(float) u,KIn(float) v,KIn(float) w,
    KIn(float) diagonal,KIn(float) residual_u,KIn(float) residual_v,KIn(float) residual_w,
    KOut(float) continuity,KOut(float) scale,KOut(float) speed,
    KOut(float) ru,KOut(float) rv,KOut(float) rw) {
  float mass=0.0f;
  for(int j=offsets[i];j<offsets[i+1];++j) mass+=sign[j]*flux[faces[j]];
  float velocity=ksqrt(u[i]*u[i]+v[i]*v[i]+w[i]*w[i]);
  continuity[i]=kabs(mass);scale[i]=kabs(diagonal[i])*velocity;speed[i]=velocity;
  ru[i]=kabs(residual_u[i]);rv[i]=kabs(residual_v[i]);rw[i]=kabs(residual_w[i]);
}

KUTTA_KERNEL void fv_continuity(int i,KIn(int) offsets,KIn(int) faces,
    KIn(float) sign,KIn(float) flux,KOut(float) result) {
  float mass=0.0f;
  for(int j=offsets[i];j<offsets[i+1];++j) mass+=sign[j]*flux[faces[j]];
  result[i]=kabs(mass);
}

KUTTA_KERNEL void fv_wall_pressure_force(int i,KIn(float) selected,KIn(float) sx,
    KIn(float) sy,KIn(float) sz,KIn(float) face_pressure,
    KOut(float) x,KOut(float) y,KOut(float) z,int internal,float reference) {
  int f=i+internal;float p=selected[i]*(face_pressure[f]-reference);
  x[i]=p*sx[f];y[i]=p*sy[f];z[i]=p*sz[f];
}

KUTTA_KERNEL void fv_wall_viscous_force(int i,KIn(int) owner,KIn(float) selected,
    KIn(float) sx,KIn(float) sy,KIn(float) sz,KIn(float) rx,KIn(float) ry,KIn(float) rz,
    KIn(float) u,KIn(float) v,KIn(float) w,
    KIn(float) ux,KIn(float) uy,KIn(float) uz,
    KIn(float) vx,KIn(float) vy,KIn(float) vz,
    KIn(float) wx,KIn(float) wy,KIn(float) wz,
    KIn(float) bu,KIn(float) bv,KIn(float) bw,
    KOut(float) x,KOut(float) y,KOut(float) z,int internal,float viscosity) {
  int f=i+internal,c=owner[f];
  float area=ksqrt(sx[f]*sx[f]+sy[f]*sy[f]+sz[f]*sz[f]);
  float nx=sx[f]/area,ny=sy[f]/area,nz=sz[f]/area;
  float normal=rx[f]*nx+ry[f]*ny+rz[f]*nz;
  float du=(bu[i]-u[c]-ux[c]*rx[f]-uy[c]*ry[f]-uz[c]*rz[f])/normal;
  float dv=(bv[i]-v[c]-vx[c]*rx[f]-vy[c]*ry[f]-vz[c]*rz[f])/normal;
  float dw=(bw[i]-w[c]-wx[c]*rx[f]-wy[c]*ry[f]-wz[c]*rz[f])/normal;
  float a=ux[c]+du*nx,b=uy[c]+du*ny,d=uz[c]+du*nz;
  float e=vx[c]+dv*nx,g=vy[c]+dv*ny,h=vz[c]+dv*nz;
  float j=wx[c]+dw*nx,k=wy[c]+dw*ny,l=wz[c]+dw*nz;
  float scale=-viscosity*selected[i];
  x[i]=scale*(2.0f*a*sx[f]+(b+e)*sy[f]+(d+j)*sz[f]);
  y[i]=scale*((e+b)*sx[f]+2.0f*g*sy[f]+(h+k)*sz[f]);
  z[i]=scale*((j+d)*sx[f]+(k+h)*sy[f]+2.0f*l*sz[f]);
}

// Fluent's coupled pseudo-time momentum normalization. This is independent
// of the physical diagonal used by Rhie-Chow and the legacy steady monitor.
KUTTA_KERNEL void fv_fluent_momentum_scale(int i,KIn(float) ax,KIn(float) ay,KIn(float) az,
    KIn(float) volume,KIn(float) speed,KOut(float) x,KOut(float) y,KOut(float) z,float rho,float dt) {
  float time=dt>0.0f?rho*volume[i]/dt:0.0f;
  x[i]=2.0f*(ax[i]+time)*speed[i];y[i]=2.0f*(ay[i]+time)*speed[i];z[i]=2.0f*(az[i]+time)*speed[i];
}
