#include <algorithm>
#include <array>

#include "../backend_helpers.hpp"
#include "../linalg/instrumented_device.hpp"
#include "kernels.hpp"
#include "kutta/fv/transport.hpp"
#include "kutta/io/mesh_io.hpp"
#include "kutta/linalg/amg.hpp"
#include "kutta/mesh/detail.hpp"

using namespace kutta;
namespace {
PolyMesh graded_mesh() {
  auto mesh=make_box_hex(4,8,1,{0,0,0},{1,8,.01});
  std::array<double,9> y{};
  for(int j=0;j<8;++j) y[j+1]=y[j]+std::pow(10.,j-8);
  for(auto& p:mesh.points) p.y=y[std::lround(p.y)];
  validate(mesh);
  CHECK((y[8]-y[7])/(y[1]-y[0])>=1e6);
  CHECK(.25/(y[1]-y[0])>=1e4);
  return mesh;
}

void check_m_matrix(const CsrPattern& p,const std::vector<float>& values,bool positive,bool symmetric=true) {
  const auto n=p.n_rows;
  std::vector<double> dense(std::size_t(n)*n),scale(n);
  for(Index i=0;i<n;++i) {
    const double diagonal=values[p.diag[i]];
    double off=0;
    for(Index j=p.row_ptr[i];j<p.row_ptr[i+1];++j) {
      const auto column=p.col[j];
      dense[std::size_t(i)*n+column]=values[j];
      if(column!=i) {CHECK(values[j]<=0);off-=double(values[j]);}
    }
    CHECK(diagonal>=off);
    REQUIRE(diagonal>0);
    scale[i]=std::sqrt(diagonal);
  }
  if(symmetric) for(Index i=0;i<n;++i) for(Index j=0;j<n;++j) {
    CHECK(dense[std::size_t(i)*n+j]==dense[std::size_t(j)*n+i]);
  }
  for(Index i=0;i<n;++i) for(Index j=0;j<n;++j)
    dense[std::size_t(i)*n+j]/=scale[i]*scale[j];
  if(!positive) return;
  // Independent double Cholesky of D^-1/2 A D^-1/2 (a congruence preserves
  // definiteness while avoiding the raw coefficient range in the check).
  for(Index i=0;i<n;++i) for(Index j=0;j<=i;++j) {
    double value=dense[std::size_t(i)*n+j];
    for(Index k=0;k<j;++k) value-=dense[std::size_t(i)*n+k]*dense[std::size_t(j)*n+k];
    if(i==j) {REQUIRE(value>0);dense[std::size_t(i)*n+j]=std::sqrt(value);}
    else dense[std::size_t(i)*n+j]=value/dense[std::size_t(j)*n+j];
  }
}

TEST_CASE("FV momentum and SST transport retain graded diffusion dominance") {
  auto mesh=graded_mesh();const auto g=compute_geometry(mesh);
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto m=upload_mesh(*dev,mesh,g);auto s=fv::make_stencil(*dev,mesh,g);
    fv::TransportEquation equation(*dev,m,s,mesh,g);
    auto field=dev->alloc<float>(m.n_cells),source=dev->alloc<float>(m.n_cells);
    auto boundary=dev->alloc<float>(m.n_boundary_faces()),mu=dev->alloc<float>(m.n_faces);
    auto kind=dev->alloc<Index>(m.n_boundary_faces());
    auto sink=dev->alloc<float>(m.n_cells),extra=dev->alloc<float>(m.n_cells);
    dev->fill(field,0);dev->fill(source,0);dev->fill(boundary,0);dev->fill(mu,.7f);dev->fill(kind,0);
    dev->fill(sink,.3f);dev->fill(extra,1e-9f);
    for(int mode=0;mode<3;++mode) {
      std::vector<float> flux(m.n_faces,0);
      if(mode==2) for(Index f=0;f<m.n_faces;++f) flux[f]=float(g.face_area[f].x);
      auto mass=dev->upload<float>(flux);
      fv::TransportSettings settings;settings.relaxation=1;
      fv::TransportTerms terms;
      if(mode==1) {terms.sink=&sink;terms.diagonal=&extra;}
      equation.assemble(field,boundary,kind,mass,mu,source,settings,terms);
      check_m_matrix(s.pattern,dev->download(equation.matrix.val),mode!=2,mode!=2);
    }
  }
}

TEST_CASE("pressure M-matrix AMG preserves stored symmetry and dominance through refresh") {
  auto mesh=graded_mesh();const auto g=compute_geometry(mesh);
  for(auto backend:Device::available()) {
    INFO(backend_name(backend));test_linalg::InstrumentedDevice dev(backend);
    auto m=upload_mesh(dev,mesh,g);auto s=fv::make_stencil(dev,mesh,g);auto matrix=upload(dev,s.pattern);
    auto flux=dev.alloc<float>(m.n_faces),rhs=dev.alloc<float>(m.n_cells);dev.fill(flux,0);
    auto coefficients=dev.download(s.diffusion);
    std::unique_ptr<LinearSolver> solver;
    for(int refresh=0;refresh<2;++refresh) {
      if(refresh) for(Index f=0;f<m.n_faces;++f) coefficients[f]*=1.13f+.07f*float(f%7);
      auto d=dev.upload<float>(coefficients);
      kern::fv_pressure_assemble(dev,m.n_cells,m.cell_face_offsets,m.cell_faces,m.cell_face_sign,
        s.adjacent,s.matrix_entry,matrix.row_ptr,matrix.diag,d,flux,matrix.val,rhs,-1);
      const auto values=dev.download(matrix.val);
      if(!solver) {
        CsrMatrixHost host{s.pattern,{values.begin(),values.end()}};
        AmgOptions options;options.structure=OperatorStructure::SymmetricMMatrix;options.coarse_rows=2;
        solver=make_amg_solver(SolverKind::AmgPcg,dev,host,options);
      }
      int coarse_levels=0;
      dev.after_launch=[&](KernelId id,std::size_t n,std::span<const KernelArg> args) {
        if(id!=KernelId::amg_m_matrix_diagonal) return;
        CsrPattern p;p.n_rows=Index(n);
        p.row_ptr=dev.read<Index>(*args[5].buffer);p.diag=dev.read<Index>(*args[6].buffer);
        // Coarse columns are captured from the first subsequent operator use
        // below; here check each stored diagonal against its actual row edges.
        const auto coarse=dev.read<float>(*args[7].buffer);
        for(Index i=0;i<p.n_rows;++i) {
          double off=0;
          for(Index j=p.row_ptr[i];j<p.row_ptr[i+1];++j) if(j!=p.diag[i]) {CHECK(coarse[j]<=0);off-=double(coarse[j]);}
          CHECK(double(coarse[p.diag[i]])>=off);
        }
        ++coarse_levels;
      };
      solver->refresh(dev,matrix);dev.after_launch={};CHECK(coarse_levels>=2);
      int used_levels=0;
      dev.after_launch=[&](KernelId id,std::size_t n,std::span<const KernelArg> args) {
        if(id==KernelId::amg_coarse_solve) n=args[7].get<Index>();
        if(id!=KernelId::amg_smooth && id!=KernelId::amg_coarse_solve) return;
        CsrPattern p;p.n_rows=Index(n);p.row_ptr=dev.read<Index>(*args[0].buffer);
        p.col=dev.read<Index>(*args[1].buffer);p.diag.resize(n);
        for(Index i=0;i<p.n_rows;++i)
          p.diag[i]=Index(std::lower_bound(p.col.begin()+p.row_ptr[i],p.col.begin()+p.row_ptr[i+1],i)-p.col.begin());
        check_m_matrix(p,dev.read<float>(*args[2].buffer),true);++used_levels;
      };
      auto x=dev.alloc<float>(m.n_cells),b=dev.alloc<float>(m.n_cells);
      dev.fill(x,1);spmv(dev,matrix,x,b);dev.fill(x,0);
      solver->solve(dev,matrix,x,b,{0,0,1});dev.after_launch={};CHECK(used_levels>0);
    }
  }
}
}

TEST_CASE("FV graded float pressure is symmetric diagonally dominant and anchored SPD") {
  auto mesh=graded_mesh();const auto g=compute_geometry(mesh);
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto m=upload_mesh(*dev,mesh,g);auto s=fv::make_stencil(*dev,mesh,g);
    auto matrix=upload(*dev,s.pattern);
    auto flux=dev->alloc<float>(m.n_faces),rhs=dev->alloc<float>(m.n_cells);
    dev->fill(flux,0);
    const auto diffusion=dev->download(s.diffusion);
    for(int mode=0;mode<3;++mode) {
      auto coefficient=diffusion;
      if(mode!=0) std::fill(coefficient.begin()+m.n_internal_faces,coefficient.end(),0);
      auto d=dev->upload<float>(coefficient);
      const int reference=mode==2?m.n_cells/2:-1;
      kern::fv_pressure_assemble(*dev,m.n_cells,m.cell_face_offsets,m.cell_faces,m.cell_face_sign,
        s.adjacent,s.matrix_entry,matrix.row_ptr,matrix.diag,d,flux,matrix.val,rhs,reference);
      const auto values=dev->download(matrix.val);
      check_m_matrix(s.pattern,values,mode!=1);
      if(reference>=0) {
        CHECK(values[s.pattern.diag[reference]]==1);
        for(Index i=0;i<m.n_cells;++i) for(Index j=s.pattern.row_ptr[i];j<s.pattern.row_ptr[i+1];++j)
          if(i!=reference && s.pattern.col[j]==reference) CHECK(values[j]==0);
      }
      // Reassembling must replace every entry and reproduce the stored bits.
      kern::fv_pressure_assemble(*dev,m.n_cells,m.cell_face_offsets,m.cell_faces,m.cell_face_sign,
        s.adjacent,s.matrix_entry,matrix.row_ptr,matrix.diag,d,flux,matrix.val,rhs,reference);
      CHECK(dev->download(matrix.val)==values);
    }
  }
}

TEST_CASE("FV pressure diagonal bounds final entries with repeated neighbour faces") {
  KUTTA_FOR_EACH_BACKEND(dev) {
    // Two shared faces and a weak third-cell edge: summing face magnitudes
    // independently does not in general bound the rounded grouped CSR entry.
    auto offsets=dev->upload<Index>(std::vector<Index>{0,3,5,6});
    auto faces=dev->upload<Index>(std::vector<Index>{0,1,2,0,1,2});
    auto sign=dev->upload<float>(std::vector<float>{1,1,1,-1,-1,-1});
    auto adjacent=dev->upload<Index>(std::vector<Index>{1,1,2,0,0,0});
    auto entry=dev->upload<Index>(std::vector<Index>{1,1,2,3,3,5});
    auto rows=dev->upload<Index>(std::vector<Index>{0,3,5,7});
    auto diag=dev->upload<Index>(std::vector<Index>{0,4,6});
    auto coefficient=dev->upload<float>(std::vector<float>{1.0f,std::ldexp(1.f,-24),std::ldexp(1.f,-25)});
    auto flux=dev->alloc<float>(3),values=dev->alloc<float>(7),rhs=dev->alloc<float>(3);
    dev->fill(flux,0);
    kern::fv_pressure_assemble(*dev,3,offsets,faces,sign,adjacent,entry,rows,diag,coefficient,flux,values,rhs,-1);
    const CsrPattern pattern{3,{0,3,5,7},{0,1,2,0,1,0,2},{0,4,6}};
    check_m_matrix(pattern,dev->download(values),false);
  }
}

TEST_CASE("pressure M-matrix coarsening retains explicit constant nullspace projection") {
  auto mesh=make_box_hex(8,8,1,{0,0,0},{1,1,.01});const auto g=compute_geometry(mesh);
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto m=upload_mesh(*dev,mesh,g);auto s=fv::make_stencil(*dev,mesh,g);auto matrix=upload(*dev,s.pattern);
    std::vector<float> coefficients(m.n_faces,0);
    std::fill(coefficients.begin(),coefficients.begin()+m.n_internal_faces,1);
    auto d=dev->upload<float>(coefficients),flux=dev->alloc<float>(m.n_faces),b=dev->alloc<float>(m.n_cells);
    dev->fill(flux,0);
    kern::fv_pressure_assemble(*dev,m.n_cells,m.cell_face_offsets,m.cell_faces,m.cell_face_sign,
      s.adjacent,s.matrix_entry,matrix.row_ptr,matrix.diag,d,flux,matrix.val,b,-1);
    const auto values=dev->download(matrix.val);CsrMatrixHost host{s.pattern,{values.begin(),values.end()}};
    auto solver=make_solver(SolverKind::AmgPcg,*dev,host,
      {.nullspace=NullSpace::Constant,.structure=OperatorStructure::SymmetricMMatrix});
    std::vector<float> exact(m.n_cells);
    for(Index i=0;i<m.n_cells;++i) exact[i]=float(i)-float(m.n_cells-1)/2;
    auto x=dev->upload<float>(exact);spmv(*dev,matrix,x,b);dev->fill(x,0);
    solver->refresh(*dev,matrix);const auto stats=solver->solve(*dev,matrix,x,b,{1e-6,0,200});
    CHECK(stats.converged);CHECK(std::abs(dev->sum(x,x.size()))<1e-5);
    auto r=dev->alloc<float>(m.n_cells);residual(*dev,matrix,x,b,r);
    CHECK(std::sqrt(dev->dot(r,r,r.size()))<=1e-6*stats.initial_residual);
    dev->fill(b,1);CHECK_THROWS(solver->solve(*dev,matrix,x,b,{}));
  }
}
