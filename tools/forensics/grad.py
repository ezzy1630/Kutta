"""Cell-centred gradients on the polyhedral mesh.

Least-squares cell-based (Fluent's default family) with selectable distance
weighting, plus cell-based Green-Gauss for a method-sensitivity check.
Boundary faces act as neighbours located at the face centroid carrying the saved
boundary face value (wall velocity = 0 in this case).
"""
import numpy as np


def neighbour_pairs(mesh, nbf):
    """(cell, dx, face_or_cell index, is_boundary) rows for LSQ assembly."""
    c0, c1, cen, fc = mesh["c0"], mesh["c1"], mesh["cen"], mesh["fcent"]
    inn = np.nonzero(c1 >= 0)[0]
    bnd = np.nonzero(c1 < 0)[0]
    assert bnd.max() < nbf
    cell = np.concatenate([c0[inn], c1[inn], c0[bnd]])
    other = np.concatenate([c1[inn], c0[inn], bnd])
    is_b = np.concatenate([np.zeros(2 * inn.size, bool), np.ones(bnd.size, bool)])
    dx = np.where(is_b[:, None], fc[other] - cen[cell], cen[np.where(is_b, 0, other)] - cen[cell])
    return cell, other, is_b, dx


def lsq_gradient(mesh, phi_cell, phi_bface, weight_power=0, pairs=None):
    """phi_cell: (n, ...) ; phi_bface: (nbf, ...) boundary-face values."""
    n = mesh["n_cells"]
    cell, other, is_b, dx = pairs if pairs is not None else neighbour_pairs(mesh, phi_bface.shape[0])
    w = np.linalg.norm(dx, axis=1) ** (-float(weight_power))
    M = np.zeros((n, 3, 3))
    np.add.at(M, cell, w[:, None, None] * dx[:, :, None] * dx[:, None, :])
    Minv = np.linalg.inv(M)
    phi_cell = np.asarray(phi_cell)
    tail = phi_cell.shape[1:]
    nb = np.where(is_b.reshape((-1,) + (1,) * len(tail)), phi_bface[np.where(is_b, other, 0)], phi_cell[np.where(is_b, 0, other)])
    dphi = nb - phi_cell[cell]
    rhs = np.zeros((n, 3) + tail)
    np.add.at(rhs, cell, (w[:, None] * dx).reshape((-1, 3) + (1,) * len(tail)) * dphi[:, None])
    return np.einsum("nij,nj...->ni...", Minv, rhs)


def green_gauss(mesh, phi_cell, phi_bface):
    """Cell-based Green-Gauss with arithmetic-mean interior face values."""
    c0, c1, A, vol, n = mesh["c0"], mesh["c1"], mesh["farea"], mesh["vol"], mesh["n_cells"]
    inn = c1 >= 0
    phi_cell = np.asarray(phi_cell)
    tail = phi_cell.shape[1:]
    nbf = phi_bface.shape[0]
    pf = np.zeros((c0.size,) + tail)
    pf[inn] = 0.5 * (phi_cell[c0[inn]] + phi_cell[c1[inn]])
    pf[:nbf] = phi_bface
    # Outward area for c0 is -A (RH normal points into c0).
    contrib = A.reshape((-1, 3) + (1,) * len(tail)) * pf[:, None]
    g = np.zeros((n, 3) + tail)
    np.add.at(g, c0, -contrib)
    np.add.at(g, c1[inn], contrib[inn])
    return g / vol.reshape((-1,) + (1,) * (len(tail) + 1))


def strain_vorticity(gradU):
    """gradU[n, i, j] = d u_j / d x_i.  Returns S = sqrt(2 SijSij), Omega = sqrt(2 WijWij)."""
    L = np.transpose(gradU, (0, 2, 1))  # L[n, j, i] = du_j/dx_i  -> velocity-gradient tensor
    Sij = 0.5 * (L + np.transpose(L, (0, 2, 1)))
    Wij = 0.5 * (L - np.transpose(L, (0, 2, 1)))
    S = np.sqrt(2.0 * np.einsum("nij,nij->n", Sij, Sij))
    W = np.sqrt(2.0 * np.einsum("nij,nij->n", Wij, Wij))
    return S, W, L
