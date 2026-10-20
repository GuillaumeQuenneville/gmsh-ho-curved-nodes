# gmsh high-order nodes on curved boundaries

Minimal reproduction: with `setOrder(k)` for k ≥ 3, the element maps of the
tetrahedra on a curved boundary have |D^k T| ∝ h² instead of h^k. That breaks
the Ciarlet–Raviart condition |D^m T| ≤ C h^m, so isoparametric finite elements
lose one order of convergence. Order 2 is fine.

With [pixi](https://pixi.sh) (gmsh 4.15.2 from conda-forge, pinned in `pixi.lock`):

```sh
pixi run repro     # unit ball, orders 2-4, h = 0.4, 0.2, 0.1 (a few seconds)
pixi run control   # same meshes, straight nodes moved by a smooth map
pixi run fine      # order 3 down to h = 0.05
pixi run optimize  # order 3 with Mesh.HighOrderOptimize = 0..4: same result
pixi run outputs   # regenerate output.txt and output-control.txt
```

Or with pip: `pip install -r requirements.txt`, then `python repro.py --help`.

## What it measures

- **Mesh:** the OCC unit sphere at mesh sizes h, then `setOrder(k)`.
- **Which elements:** the tetrahedra with a high-order node on the sphere.
- **Metric:** the median over those elements of the largest component of D^k T.
  - T has degree k, so D^k T is constant on each element.
  - It is computed exactly, as a k-th forward difference of the Lagrange basis
    from `getBasisFunctions`.
- **Control:** `--control` places straight-sided nodes (`Mesh.SecondOrderLinear = 1`)
  and moves every node by x → x + 0.1 sin(2x). That gives rates of exactly k, so
  the measurement itself is sound.
- **Node report:** for each element, the script also says which of the other
  high-order nodes gmsh left at their straight-sided position.

## Results (gmsh 4.15.2, Linux)

Full output: [`output.txt`](output.txt) and [`output-control.txt`](output-control.txt).

| k | gmsh nodes, rate of \|D^k T\| | control, rate | expected |
|---|---|---|---|
| 2 | 2.07, 1.97 | 1.93, 1.97 | 2 |
| 3 | 2.28, 2.04 (1.98 at h = 0.05) | 3.19, 3.07 | 3 |
| 4 | 2.15, 1.99 | 4.11, 3.96 | 4 |

Where gmsh places the other high-order nodes of the curved elements:
- **Interior edges:** every node on them is at its straight position.
- **Interior faces with a boundary edge:** their face nodes are moved.
- **Element interiors (order 4):** the volume nodes are moved.

The way those nodes are moved is not consistent with a smooth map, so
|D^k T| stays at h².
