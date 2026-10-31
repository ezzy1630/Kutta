#include "kutta/kernels/prelude.h"

// Weighted quadratic moments sum_j w_j b(x_j-x_c)(q_j-q_c), basis
// b(r)=(x,y,z,x^2/2,xy,xz,y^2/2,yz,z^2/2,1). The weights are walk counts of
// (I+A)^hops, so one gather per hop builds every cell's moment from its face
// neighbours. Rows are cell-major so a neighbour gather touches one or two
// cache lines: a moment row holds the ten moments and then the cloud's min and
// max, a geometry row the ten moments of the constant 1, a coefficient row the
// three functionals' ten terms.
#define KMOMENT_ROW 12
struct kmoment { float v[10]; };

// Moment about x_d expressed about x_c, with d=x_d-x_c.
KUTTA_FN kmoment kmoment_shift(kmoment m,float dx,float dy,float dz) {
  kmoment o;float c=m.v[9];
  o.v[0]=m.v[0]+dx*c;o.v[1]=m.v[1]+dy*c;o.v[2]=m.v[2]+dz*c;
  o.v[3]=m.v[3]+dx*m.v[0]+0.5f*dx*dx*c;
  o.v[4]=m.v[4]+dy*m.v[0]+dx*m.v[1]+dx*dy*c;
  o.v[5]=m.v[5]+dz*m.v[0]+dx*m.v[2]+dx*dz*c;
  o.v[6]=m.v[6]+dy*m.v[1]+0.5f*dy*dy*c;
  o.v[7]=m.v[7]+dz*m.v[1]+dy*m.v[2]+dy*dz*c;
  o.v[8]=m.v[8]+dz*m.v[2]+0.5f*dz*dz*c;
  o.v[9]=c;
  return o;
}
// Moment of the anchored side, M_s+(q_s-anchor)G_s, moved by (dx,dy,dz).
KUTTA_FN kmoment kmoment_side(KIn(float) moment,KIn(float) geometry,int s,float offset,float dx,float dy,float dz) {
  kmoment m;for(int k=0;k<10;++k) m.v[k]=moment[s*KMOMENT_ROW+k]+offset*geometry[s*10+k];
  return kmoment_shift(m,dx,dy,dz);
}
KUTTA_FN float kmoment_apply(KIn(float) coefficients,int slot,int functional,kmoment m) {
  kupper_sum sum={0.0f,0.0f};
  for(int k=0;k<10;++k) sum=kupper_add(sum,coefficients[slot*30+functional*10+k]*m.v[k]);
  return kupper_value(sum);
}

// Hop 0: the cell itself and its independently prescribed boundary values.
KUTTA_KERNEL void fv_moment_seed(int i,KIn(int) offsets,KIn(int) adjacent,KIn(int) observed,
    KIn(float) rx,KIn(float) ry,KIn(float) rz,KIn(float) field,KIn(float) boundary,
    KOut(float) moment,int internal) {
  kmoment m;for(int k=0;k<10;++k) m.v[k]=0.0f;
  float q=field[i],low=q,high=q;
  for(int j=offsets[i];j<offsets[i+1];++j) {
    int b=adjacent[j];if(b>=0) continue;b=-b-1;if(observed[b]==0) continue;
    int f=b+internal;float x=rx[f],y=ry[f],z=rz[f],d=boundary[b]-q;
    m.v[0]+=x*d;m.v[1]+=y*d;m.v[2]+=z*d;m.v[3]+=0.5f*x*x*d;m.v[4]+=x*y*d;m.v[5]+=x*z*d;
    m.v[6]+=0.5f*y*y*d;m.v[7]+=y*z*d;m.v[8]+=0.5f*z*z*d;m.v[9]+=d;
    low=kmin(low,boundary[b]);high=kmax(high,boundary[b]);
  }
  for(int k=0;k<10;++k) moment[i*KMOMENT_ROW+k]=m.v[k];
  moment[i*KMOMENT_ROW+10]=low;moment[i*KMOMENT_ROW+11]=high;
}

// M_c <- M_c + sum_d shift(M_d+(q_d-q_c)G_d, x_d-x_c), in fixed face order.
KUTTA_KERNEL void fv_moment_hop(int i,KIn(int) offsets,KIn(int) faces,KIn(int) adjacent,KIn(float) sign,
    KIn(float) dx,KIn(float) dy,KIn(float) dz,KIn(float) field,KIn(float) geometry,
    KIn(float) moment,KOut(float) next) {
  kmoment m;for(int k=0;k<10;++k) m.v[k]=moment[i*KMOMENT_ROW+k];
  float q=field[i],low=moment[i*KMOMENT_ROW+10],high=moment[i*KMOMENT_ROW+11];
  for(int j=offsets[i];j<offsets[i+1];++j) {
    int d=adjacent[j];if(d<0) continue;int f=faces[j];float s=sign[j];
    kmoment side=kmoment_side(moment,geometry,d,field[d]-q,s*dx[f],s*dy[f],s*dz[f]);
    for(int k=0;k<10;++k) m.v[k]+=side.v[k];
    low=kmin(low,moment[d*KMOMENT_ROW+10]);high=kmax(high,moment[d*KMOMENT_ROW+11]);
  }
  for(int k=0;k<10;++k) next[i*KMOMENT_ROW+k]=m.v[k];
  next[i*KMOMENT_ROW+10]=low;next[i*KMOMENT_ROW+11]=high;
}

KUTTA_FN kmoment kmoment_face(int i,KIn(int) owner,KIn(int) neighbour,
    KIn(float) px,KIn(float) py,KIn(float) pz,KIn(float) nx,KIn(float) ny,KIn(float) nz,
    KIn(float) field,KIn(float) geometry,KIn(float) moment,int internal,float anchor) {
  int p=owner[i];
  kmoment m=kmoment_side(moment,geometry,p,field[p]-anchor,-px[i],-py[i],-pz[i]);
  if(i<internal) {
    int n=neighbour[i];
    kmoment o=kmoment_side(moment,geometry,n,field[n]-anchor,-nx[i],-ny[i],-nz[i]);
    for(int k=0;k<10;++k) m.v[k]+=o.v[k];
  }
  return m;
}

// Functional 0: the deferred tangential diffusion correction (per unit
// diffusivity). Uncorrected faces write zero.
KUTTA_KERNEL void fv_moment_diffusion(int i,KIn(int) owner,KIn(int) neighbour,KIn(int) slot,
    KIn(float) px,KIn(float) py,KIn(float) pz,KIn(float) nx,KIn(float) ny,KIn(float) nz,
    KIn(float) field,KIn(float) geometry,KIn(float) moment,KIn(float) coefficients,
    KOut(float) correction,int internal) {
  int s=slot[i];if(s<0) {correction[i]=0.0f;return;}
  int p=owner[i];float anchor=i<internal?0.5f*(field[p]+field[neighbour[i]]):field[p];
  kmoment m=kmoment_face(i,owner,neighbour,px,py,pz,nx,ny,nz,field,geometry,moment,internal,anchor);
  correction[i]=kmoment_apply(coefficients,s,0,m);
}

// Functionals 1 and 2: normal-line face value and centroid link, as in
// fv_projected_momentum_faces.
KUTTA_KERNEL void fv_moment_momentum(int i,KIn(int) owner,KIn(int) neighbour,KIn(int) slot,
    KIn(float) px,KIn(float) py,KIn(float) pz,KIn(float) nx,KIn(float) ny,KIn(float) nz,
    KIn(float) field,KIn(float) geometry,KIn(float) moment,KIn(float) coefficients,
    KIn(float) diagonal,KIn(float) cell_limiter,KOut(float) face,int internal,int limit) {
  int s=slot[i];if(s<0) return;
  int p=owner[i],n=neighbour[i];float w=diagonal[p]/(diagonal[p]+diagonal[n]);
  float cp=field[p],anchor=0.5f*(cp+field[n]);
  kmoment m=kmoment_face(i,owner,neighbour,px,py,pz,nx,ny,nz,field,geometry,moment,internal,anchor);
  float value=w*cp+(1.0f-w)*field[n]+kmoment_apply(coefficients,s,1,m)-(1.0f-w)*kmoment_apply(coefficients,s,2,m);
  if(limit!=0) {
    float theta=kmin(cell_limiter[p],cell_limiter[n]);
    float low=kmin(moment[p*KMOMENT_ROW+10],moment[n*KMOMENT_ROW+10]);
    float high=kmax(moment[p*KMOMENT_ROW+11],moment[n*KMOMENT_ROW+11]);
    face[i]+=theta*(kmin(high,kmax(low,value))-face[i]);
  } else face[i]=value;
}

// Resolved walls: the momentum wall flux per unit viscosity is
// diffusion*(u_b-u_P)+C. The legacy force's normal derivative equals that
// flux with C=g_P.t; on projected faces shift it by (C-g_P.t)/|S|.
KUTTA_KERNEL void fv_wall_projected_shift(int i,KIn(int) owner,KIn(int) mask,KIn(float) sx,KIn(float) sy,KIn(float) sz,
    KIn(float) tx,KIn(float) ty,KIn(float) tz,KIn(float) gx,KIn(float) gy,KIn(float) gz,KIn(float) correction,
    KOut(float) shift,int internal) {
  int f=i+internal;
  if(mask[f]==0) {shift[i]=0.0f;return;}
  int p=owner[f];float area=ksqrt(sx[f]*sx[f]+sy[f]*sy[f]+sz[f]*sz[f]);
  shift[i]=(correction[f]-(gx[p]*tx[f]+gy[p]*ty[f]+gz[p]*tz[f]))/area;
}
// Adds the Newtonian traction of normal-derivative shifts (du,dv,dw).
KUTTA_KERNEL void fv_wall_normal_force_shift(int i,KIn(float) selected,KIn(float) sx,KIn(float) sy,KIn(float) sz,
    KIn(float) du,KIn(float) dv,KIn(float) dw,KOut(float) x,KOut(float) y,KOut(float) z,int internal,float viscosity) {
  int f=i+internal;float area=ksqrt(sx[f]*sx[f]+sy[f]*sy[f]+sz[f]*sz[f]);
  float nx=sx[f]/area,ny=sy[f]/area,nz=sz[f]/area,scale=-viscosity*selected[i];
  float a=du[i]*nx,b=du[i]*ny,d=du[i]*nz,e=dv[i]*nx,g=dv[i]*ny,h=dv[i]*nz,j=dw[i]*nx,k=dw[i]*ny,l=dw[i]*nz;
  x[i]+=scale*(2.0f*a*sx[f]+(b+e)*sy[f]+(d+j)*sz[f]);
  y[i]+=scale*((e+b)*sx[f]+2.0f*g*sy[f]+(h+k)*sz[f]);
  z[i]+=scale*((j+d)*sx[f]+(k+h)*sy[f]+2.0f*l*sz[f]);
}
