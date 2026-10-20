# gmsh high-order nodes on curved boundaries

Minimal reproduction: with `setOrder(k)` for k ≥ 3, the element maps of the
tetrahedra on a curved boundary have |D^k T| ∝ h² instead of h^k. That breaks
the Ciarlet–Raviart condition |D^m T| ≤ C h^m, so isoparametric finite elements
lose one order of convergence. Order 2 is fine, and so are the boundary
triangles themselves at every order: the defect is in how the boundary
curvature is carried into the tetrahedra. [Background](#background) explains
the terms.

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
- **Which elements:**
  - the tetrahedra with a high-order node on the sphere;
  - separately, the boundary triangles on their own.
- **Metric:** the median over those elements of the largest component of D^k T.
  - T has degree k, so D^k T is constant on each element.
  - It is computed exactly, as a k-th forward difference of the Lagrange basis
    from `getBasisFunctions`.
- **Control:** `--control` places straight-sided nodes (`Mesh.SecondOrderLinear = 1`)
  and moves every node by x → x + 0.1 sin(2x). That gives rates of exactly k, so
  the measurement itself is sound.
- **Node report:** for the curved tetrahedra, the script also says which of the
  other high-order nodes gmsh left at their straight-sided position.

## Results (gmsh 4.15.2, Linux)

Full output: [`output.txt`](output.txt) and [`output-control.txt`](output-control.txt).
Rates of |D^k T| for h = 0.4 → 0.2 → 0.1:

| k | tetrahedra, gmsh nodes | boundary triangles, gmsh nodes | tetrahedra, control | expected |
|---|---|---|---|---|
| 2 | 2.07, 1.97 | 2.03, 2.01 | 1.93, 1.97 | 2 |
| 3 | 2.28, 2.04 (1.98 at h = 0.05) | 3.02, 2.98 | 3.19, 3.07 | 3 |
| 4 | 2.15, 1.99 | 4.07, 3.99 | 4.11, 3.96 | 4 |

Where gmsh places the other high-order nodes of the curved tetrahedra:
- **Interior edges:** every node on them is at its straight position.
- **Interior faces with a boundary edge:** their face nodes are moved.
- **Element interiors (order 4):** the volume nodes are moved.

## Background

### The element map T

Each element K is the image of the reference tetrahedron $\hat K$
(vertices 0, $e_1$, $e_2$, $e_3$) under its element map

$$T(\xi) = \sum_a x_a \hat\varphi_a(\xi),$$

where the $x_a$ are the element's nodes as gmsh gives them and the
$\hat\varphi_a$ are the order-k Lagrange basis functions. So T is a vector-valued
polynomial of degree k in the reference coordinates ξ. Isoparametric finite
elements use the same basis for the solution and for T.

### What D^m T is

$D^m T$ is the m-th derivative of T with respect to ξ: the tensor of all
partial derivatives $\partial^m T_i / \partial\xi_{j_1}\cdots\partial\xi_{j_m}$.
- $DT$ is the Jacobian.
- $D^2T$ measures how far the element is from affine: it is zero on a
  straight-sided element.
- For a map of degree k, $D^k T$ is constant and $D^{k+1}T = 0$.

The script reports |D^m T| as the largest component. Any norm gives the same rates.

### Why |D^m T| should scale like h^m

Suppose the nodes of every element are sampled from one smooth map F, for
example the exact curved geometry parametrized over a straight-sided element
$A(\xi)$ of size h. Then $T \approx F \circ A$, with A affine and $|DA| \sim h$,
and the chain rule gives

$$D^m T \approx D^m F\,[DA, \dots, DA] \sim |D^m F|\; h^m.$$

So each additional derivative of a well-curved element shrinks by one factor of h.
The `--control` run does exactly this, and it measures rates of exactly m.

### Why it matters

Finite element error estimates on curved elements pull the solution back to the
reference element, $\hat u = u \circ T$. Bounding the derivatives of $\hat u$ up
to order k + 1 with the chain rule brings in every $D^m T$ with m ≤ k. The
optimal estimate $|u - I_h u|_{H^1} \le C h^k$ holds when $|D^m T| \le C h^m$ for
all m ≤ k (Ciarlet & Raviart 1972; Lenoir 1986). Lenoir also shows how to build
meshes that satisfy it.

At k = 2 the condition stops at $D^2T \sim h^2$, which any placement of the
boundary nodes on the surface gives. From k = 3 it also needs $D^3 T \sim h^3$.
In our isoparametric Poisson solver on the unit ball, the gmsh nodes give these
H1 rates:
- 2.78 instead of 3 at k = 3;
- 2.75 instead of 4 at k = 4.

### Where the h² comes from

Write $T = A + E$, where A is the affine map through the four vertices and E is
the displacement of the high-order nodes from their straight positions.
- **On a boundary face**, the nodes sit on the sphere. The displacement δ there is
  $O(h^2)$ (curvature times h²), and its derivatives scale correctly. That's the
  boundary-triangle column above.
- **Split δ by degree:** $\delta = \delta_2 + \delta_{\ge 3}$, where $\delta_2$
  is its quadratic part, of size $h^2$, and the rest is $O(h^3)$.
- **Extension into the tetrahedron:** gmsh then moves the face and volume nodes.
  If that extension reproduced polynomials, it would map the quadratic $\delta_2$
  to a quadratic inside the element. Then $D^3 E$ would only see the $O(h^3)$
  part, and $D^3 T \sim h^3$.
- **What the measurement shows:** $D^3 T \sim h^2$ means the quadratic part of
  size $h^2$ comes out of the extension with a cubic (and, at k = 4, quartic)
  component of the same size.

Straight interior edges are not the problem by themselves. The vertices of a
boundary face don't move, so $\delta_2 = \sum_{i<j} c_{ij} \lambda_i \lambda_j$ in
the face's barycentric coordinates. Its natural extension uses the same formula
with the tetrahedron's barycentric coordinates, and that extension is zero on
every interior edge too. What matters is the values given to the face and volume
nodes.

### What fixes it

Place every node that is not on a boundary entity by a polynomial extension of
the boundary displacement, using the edge and face blending functions of
Szabó–Babuška or Lenoir's construction.
- In a sphere-only prototype, this gave $|D^3T| \sim h^{3.06}$.
- It also gave optimal H1 rates: 3.36 at k = 3 and 4.72 at k = 4.

### References

- P. G. Ciarlet, P.-A. Raviart, "Interpolation theory over curved elements, with
  applications to finite element methods", *Comput. Methods Appl. Mech. Engrg.*
  1 (1972) 217–249.
- M. Lenoir, "Optimal isoparametric finite elements and error estimates for
  domains involving curved boundaries", *SIAM J. Numer. Anal.* 23 (1986) 562–580.
- B. Szabó, I. Babuška, *Finite Element Analysis*, Wiley (1991), chapter on
  mapping and blending functions.
