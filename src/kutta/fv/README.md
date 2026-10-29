# Finite volume

`upload_mesh()` preserves the caller's cell/face numbering. Call
`renumber_cells_rcm()` before computing geometry or fields when loading a case.
The shared `device_mesh.hpp` contract is unchanged. All device geometry is
stored as float differences computed from host double geometry.

Wall distances use a host double BVH over the same arithmetic-center triangle
fans as mesh geometry. The nearest point can lie in a triangle, on an edge or
at a vertex; an infinite wall plane is not substituted for the finite surface.
No wall patches means infinite wall distance. Non-wall normal distances are
zero. Warped polygons follow the mesh's triangle-fan surface convention.

`make_stencil()` builds the sorted cell-adjacency CSR pattern and incidence to
CSR maps once. The LSQ normal matrix uses weights `1/|d|^2`, with neighbour
centroids on internal faces and face centroids on boundaries. Its symmetric
inverse is stored in six SoA planes. `gradient()` takes actual boundary face
values; callers must supply extrapolated values for Neumann conditions.
`limit_gradient()` computes the cell-to-face Barth–Jespersen scalar separately
so diffusion and pressure can use the unlimited gradient.

Transport uses implicit first-order upwind with deferred second-order upwind
and over-relaxed non-orthogonal diffusion. Assembly gathers into each owned
CSR row, including multiple faces between a cell pair. The linear solver
policy is centralized in `transport.hpp`. All operator work remains on Device;
only the linear solver and explicitly requested monitors synchronize with the
host.
