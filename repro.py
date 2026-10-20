"""High-order gmsh tetrahedra on a curved boundary: the element maps are not
smooth enough from order 3.

For a curved isoparametric element of order k, optimal convergence needs the
element map T to satisfy |D^m T| <= C h^m for m <= k (Ciarlet & Raviart, 1972).
This script meshes the unit sphere at several sizes, calls setOrder(k), and
measures how the k-th derivative of T scales with h on the tetrahedra that have
a high-order node on the sphere, and on the boundary triangles alone. It then
reports, at the finest size, which of the other high-order nodes of those
tetrahedra gmsh left at their straight-sided (affine) position.

With --control, the high-order nodes are placed on straight edges
(Mesh.SecondOrderLinear = 1) and every node is then moved by the smooth map
x -> x + 0.1 sin(2 x). Those element maps interpolate one smooth map, so
|D^k T| must scale like h^k: this checks the measurement itself.

Usage:
    python repro.py [--orders 2 3 4] [--sizes 0.4 0.2 0.1]
                    [--optimize 0 [1 ...]] [--control]
"""

from __future__ import annotations

import argparse
import itertools
import math

import gmsh
import numpy as np


def mesh_sphere(h: float, order: int, optimize: int, linear: bool) -> None:
    """Meshes the unit ball with size h and elevates the mesh to `order`."""
    gmsh.clear()
    gmsh.model.add("sphere")
    gmsh.model.occ.addSphere(0, 0, 0, 1)
    gmsh.model.occ.synchronize()
    gmsh.option.setNumber("Mesh.MeshSizeMin", h)
    gmsh.option.setNumber("Mesh.MeshSizeMax", h)
    gmsh.option.setNumber("Mesh.HighOrderOptimize", optimize)
    gmsh.option.setNumber("Mesh.SecondOrderLinear", int(linear))
    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.setOrder(order)


def elements(name: str, order: int) -> tuple[np.ndarray, np.ndarray, int]:
    """Node coordinates of every element of a family ("Tetrahedron" or
    "Triangle"), (ne, nn, 3), their node tags, (ne, nn), and the element type."""
    etype = gmsh.model.mesh.getElementType(name, order)
    _, node_tags = gmsh.model.mesh.getElementsByType(etype)
    tags, coords, _ = gmsh.model.mesh.getNodes()
    xyz = np.zeros((int(tags.max()) + 1, 3))
    xyz[tags.astype(int)] = coords.reshape(-1, 3)
    nn = gmsh.model.mesh.getElementProperties(etype)[3]
    conn = node_tags.astype(int).reshape(-1, nn)
    return xyz[conn], conn, etype


def on_sphere(conn: np.ndarray) -> np.ndarray:
    """Whether each node of each element lies on the sphere, (ne, nn)."""
    bnd = []
    for dim, tag in gmsh.model.getEntities(2):
        tags, _, _ = gmsh.model.mesh.getNodes(dim, tag, includeBoundary=True)
        bnd.append(tags.astype(int))
    return np.isin(conn, np.concatenate(bnd))


def max_derivative(X: np.ndarray, etype: int, m: int, dim: int = 3) -> np.ndarray:
    """Largest component of D^m T in each element, for maps of degree m from a
    reference element of dimension `dim`.

    D^m T is constant for a map of degree m, and a forward difference of order
    m is exact for it whatever the step.
    """
    p0, s = np.r_[np.full(dim, 0.1), np.zeros(3 - dim)], 0.1
    out = np.zeros(len(X))
    for dirs in itertools.combinations_with_replacement(range(dim), m):
        pts, signs = [], []
        for steps in itertools.product((0, 1), repeat=m):
            p = p0.copy()
            for d, on in zip(dirs, steps):
                p[d] += s * on
            pts.append(p)
            signs.append((-1) ** (m - sum(steps)))
        pts = np.array(pts)
        _, bf, _ = gmsh.model.mesh.getBasisFunctions(etype, pts.ravel(), "Lagrange")
        T = np.einsum("pn,enc->epc", np.reshape(bf, (len(pts), -1)), X)
        D = np.einsum("p,epc->ec", np.array(signs, float), T) / s**m
        out = np.maximum(out, np.abs(D).max(1))
    return out


def node_report(X: np.ndarray, etype: int, curved: np.ndarray, bnd: np.ndarray):
    """Share of straight-placed nodes among the high-order nodes of the curved
    elements that are not on the sphere, by the entity they lie on.

    Yields (description, number of nodes, share at the affine position).
    """
    ref = np.reshape(gmsh.model.mesh.getElementProperties(etype)[4], (-1, 3))
    affine = X[:, None, 0] + np.einsum("nd,edc->enc", ref, X[:, 1:4] - X[:, None, 0])
    size = np.linalg.norm(X[:, 1] - X[:, 0], axis=1)[:, None]
    straight = np.linalg.norm(X - affine, axis=2) <= 1e-10 * size
    bary = np.c_[1 - ref.sum(1), ref] > 1e-9  # vertices of the node's entity
    vert_bnd = bnd[:, :4]
    kinds = {2: "edge", 3: "face", 4: "volume"}
    for n_vert in (2, 3, 4):
        cols = np.where(bary.sum(1) == n_vert)[0]
        if not len(cols):
            continue
        # number of the entity's vertices that lie on the sphere
        on = np.stack([vert_bnd[:, bary[c]].sum(1) for c in cols], 1)
        mask = curved[:, None] & ~bnd[:, cols]
        for b in range(n_vert + 1):
            sel = mask & (on == b)
            if sel.any():
                what = f"{kinds[n_vert]} nodes, {b} of the {kinds[n_vert]}'s"
                what += f" {n_vert} vertices on the sphere"
                yield what, int(sel.sum()), float(straight[:, cols][sel].mean())


def table(k: int, sizes: list[float], optimize: int, control: bool) -> None:
    """Prints the scaling of |D^k T| and the node report for one order."""
    print(f"\norder {k}: median of max |D^{k} T| over the curved tets and over")
    print("the boundary triangles")
    print(f"{'h':>6} {'tets':>8} {'curved':>7} {'|D^k T|':>10} {'rate':>6}"
          f" {'triangles':>10} {'|D^k T|':>10} {'rate':>6}")
    prev = None
    for h in sizes:
        mesh_sphere(h, k, optimize, linear=control)
        X, conn, etype = elements("Tetrahedron", k)
        Xt, _, ttype = elements("Triangle", k)
        bnd = on_sphere(conn)
        curved = bnd[:, 4:].any(1)  # a high-order node on the sphere
        if control:
            X = X + 0.1 * np.sin(2 * X)
            Xt = Xt + 0.1 * np.sin(2 * Xt)
        d = (
            float(np.median(max_derivative(X[curved], etype, k))),
            float(np.median(max_derivative(Xt, ttype, k, dim=2))),
        )
        rate = ["", ""]
        if prev is not None:
            rate = [f"{math.log(p / q) / math.log(prev[0] / h):6.2f}"
                    for p, q in zip(prev[1], d)]
        print(f"{h:6.3f} {len(X):8d} {int(curved.sum()):7d} {d[0]:10.3e} {rate[0]:>6}"
              f" {len(Xt):10d} {d[1]:10.3e} {rate[1]:>6}")
        prev = (h, d)
    print(f"expected rate: {k}")
    print(f"high-order nodes of curved tets not on the sphere (h = {h}):")
    for what, n, share in node_report(X, etype, curved, bnd):
        print(f"  {what}: {n}, {share:.0%} at the straight position")


def main() -> None:
    """Runs the measurement and prints one table per order and optimize value."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--orders", type=int, nargs="+", default=[2, 3, 4])
    ap.add_argument("--sizes", type=float, nargs="+", default=[0.4, 0.2, 0.1])
    ap.add_argument(
        "--optimize", type=int, nargs="+", default=[0], help="Mesh.HighOrderOptimize"
    )
    ap.add_argument("--control", action="store_true", help="smooth-map control")
    args = ap.parse_args()

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    for optimize in args.optimize:
        head = f"gmsh {gmsh.__version__}, Mesh.HighOrderOptimize = {optimize}"
        if args.control:
            head += ", CONTROL: straight nodes moved by a smooth map"
        print(head)
        for k in args.orders:
            table(k, args.sizes, optimize, args.control)
        print()
    gmsh.finalize()


if __name__ == "__main__":
    main()
