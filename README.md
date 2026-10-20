# gmsh high-order nodes on curved boundaries

Minimal reproduction: with `setOrder(k)` for $`k \ge 3`$, the element maps $`T`$
of the tetrahedra on a curved boundary have $`\lvert D^k T \rvert \propto h^2`$
instead of $`h^k`$. That breaks the Ciarlet–Raviart condition
$`\lvert D^m T \rvert \le C h^m`$ and caps the $`H^1`$ convergence rate of
isoparametric finite elements at about $`h^{5/2}`$ for every $`k \ge 3`$.
Order 2 is fine, and so are the boundary triangles
themselves at every order: the defect is in how the boundary curvature is
carried into the tetrahedra. A degree-preserving replacement for gmsh's
inner-node weights restores the rates ([What fixes it](#what-fixes-it)), and
[Background](#background) explains the terms.

With [pixi](https://pixi.sh) (gmsh 4.15.2 from conda-forge, pinned in `pixi.lock`):

```sh
pixi run repro     # unit ball, orders 2-4, h = 0.4, 0.2, 0.1 (a few seconds)
pixi run control   # same meshes, straight nodes moved by a smooth map
pixi run fine      # order 3 down to h = 0.05
pixi run optimize  # order 3 with Mesh.HighOrderOptimize = 0..4: same result
pixi run extension # interior nodes recomputed by a polynomial extension (the fix)
pixi run outputs   # regenerate the three output files
```

Or with pip: `pip install -r requirements.txt`, then `python repro.py --help`.

## What it measures

- **Mesh:** the OCC unit sphere at mesh sizes $`h`$, then `setOrder(k)`.
- **Which elements:**
  - the tetrahedra with a high-order node on the sphere;
  - separately, the boundary triangles on their own.
- **Metric:** the median over those elements of the largest component of
  $`D^k T`$.
  - $`T`$ has degree $`k`$, so $`D^k T`$ is constant on each element.
  - It is computed exactly, as a $`k`$-th forward difference of the Lagrange
    basis from `getBasisFunctions`.
- **Control:** `--control` places straight-sided nodes
  (`Mesh.SecondOrderLinear = 1`) and moves every node by
  $`x \mapsto x + 0.1 \sin(2x)`$. That gives rates of exactly $`k`$, so the
  measurement itself is sound.
- **Node report:** for the curved tetrahedra, the script also says which of the
  other high-order nodes gmsh left at their straight-sided position.

## Results (gmsh 4.15.2, Linux)

Full output: [`output.txt`](output.txt), [`output-control.txt`](output-control.txt)
and [`output-extension.txt`](output-extension.txt).
Rates of $`\lvert D^k T \rvert`$, each cell giving the rate for $`h`$ from 0.4 → 0.2, then 0.2 → 0.1:

| $`k`$ | tetrahedra (gmsh nodes) | boundary triangles (gmsh nodes) | tetrahedra (control) | tetrahedra (extension) | expected |
|---|---|---|---|---|---|
| 2 | 2.07, 1.97 | 2.03, 2.01 | 1.93, 1.97 | 2.07, 1.97 | 2 |
| 3 | 2.28, 2.04 (0.1 → 0.05: 1.98) | 3.02, 2.98 | 3.19, 3.07 | 3.06, 2.97 | 3 |
| 4 | 2.15, 1.99 | 4.07, 3.99 | 4.11, 3.96 | 4.08, 3.95 | 4 |

Where gmsh places the other high-order nodes of the curved tetrahedra:
- **Interior edges:** every node on them is at its straight position.
- **Interior faces with a boundary edge:** their face nodes are moved.
- **Element interiors (order 4):** the volume nodes are moved.

## Background

### The element map

Each element $`K`$ is the image of a fixed reference tetrahedron $`\hat K`$
under its element map

```math
T(\xi) = \sum_a x_a \, \hat\varphi_a(\xi),
```

where the $`x_a`$ are the element's nodes as gmsh gives them and the
$`\hat\varphi_a`$ are the order-$`k`$ Lagrange basis functions. So $`T`$ is a
vector-valued polynomial of degree $`k`$ in the reference coordinates $`\xi`$.
Isoparametric finite elements use the same basis for the solution and for $`T`$.

### The derivatives of the element map

$`D^m T`$ is the $`m`$-th derivative of $`T`$ with respect to $`\xi`$: the
tensor of all partial derivatives
$`\partial^m T_i / \partial\xi_{j_1} \cdots \partial\xi_{j_m}`$.
- $`DT`$ is the Jacobian.
- $`D^2 T`$ measures how far the element is from affine: it is zero on a
  straight-sided element.
- For a map of degree $`k`$, $`D^k T`$ is constant and $`D^{k+1} T = 0`$.

The script reports $`\lvert D^m T \rvert`$ as the largest component. Any norm
gives the same rates.

### Why the derivatives should scale like powers of h

Suppose the nodes of every element are sampled from one smooth map $`F`$, for
example the exact curved geometry parametrized over a straight-sided element
$`A(\xi)`$ of size $`h`$. Then $`T \approx F \circ A`$, with $`A`$ affine and
$`\lvert DA \rvert \sim h`$, and the chain rule gives

```math
D^m T \approx D^m F\,[DA, \dots, DA] \sim \lvert D^m F \rvert \, h^m .
```

So each additional derivative of a well-curved element shrinks by one factor of
$`h`$. The `--control` run does exactly this, and it measures rates of exactly
$`m`$.

### Why it matters

Finite element error estimates on curved elements pull the solution back to the
reference element, $`\hat u = u \circ T`$. Bounding the derivatives of
$`\hat u`$ up to order $`k+1`$ with the chain rule brings in every $`D^m T`$
with $`m \le k`$. The optimal estimate

```math
\lvert u - I_h u \rvert_{H^1} \le C h^k
```

holds when $`\lvert D^m T \rvert \le C h^m`$ for all $`m \le k`$ (Ciarlet &
Raviart 1972; Lenoir 1986). Lenoir also shows how to build meshes that satisfy it.

At $`k = 2`$ the condition stops at $`D^2 T \sim h^2`$, which any placement of
the boundary nodes on the surface gives. From $`k = 3`$ it also needs
$`D^3 T \sim h^3`$. Without it, the $`H^1`$ rate is capped at about
$`h^{5/2}`$ for every $`k \ge 3`$, so the loss grows with $`k`$. In our
isoparametric Poisson solver on the unit ball, the gmsh nodes give these
$`H^1`$ rates:
- 2.78 instead of 3 at $`k = 3`$;
- 2.75 instead of 4 at $`k = 4`$.

### Where the h² comes from

Write $`T = A + E`$, where $`A`$ is the affine map through the four vertices and
$`E`$ is the displacement of the high-order nodes from their straight positions.
- **On a boundary face**, the nodes sit on the sphere. The displacement
  $`\delta`$ there is $`O(h^2)`$ (curvature times $`h^2`$), and its derivatives
  scale correctly. That's the boundary-triangle column above.
- **Split by degree:** $`\delta = \delta_2 + \delta_{\ge 3}`$, where
  $`\delta_2`$ is its quadratic part, of size $`h^2`$, and the rest is
  $`O(h^3)`$.
- **Extension into the tetrahedron:** gmsh then moves the face and volume nodes.
  If that extension reproduced polynomials, it would map the quadratic
  $`\delta_2`$ to a quadratic inside the element. Then $`D^3 E`$ would only see
  the $`O(h^3)`$ part, and $`D^3 T \sim h^3`$.
- **What the measurement shows:** $`D^3 T \sim h^2`$ means the quadratic part of
  size $`h^2`$ comes out of the extension with a cubic (and, at $`k = 4`$,
  quartic) component of the same size.

Straight interior edges are not the problem by themselves. The vertices of a
boundary face don't move, so $`\delta_2 = \sum_{i<j} c_{ij} \lambda_i \lambda_j`$
in the face's barycentric coordinates $`\lambda_i`$. Its natural extension uses
the same formula with the tetrahedron's barycentric coordinates, and that
extension is zero on every interior edge too. What matters is the values given
to the face and volume nodes.

### What fixes it

gmsh places the interior face and volume nodes with fixed weights from
`getInnerVertexPlacement` (`src/numeric/InnerVertexPlacement.cpp`): a linear
transfinite blend of the element's boundary nodes. That blend does not
reproduce a quadratic edge displacement. For the bump $`\lambda_1 \lambda_2`$ on
one edge of an order-3 triangle, it puts the face node at 0.148 instead of the
polynomial's 0.111.

`--extension` keeps gmsh's vertices, edge nodes and boundary nodes, and
recomputes the interior face and volume nodes from a degree-preserving
extension:
- **Edges:** each edge's displacement is fitted with the Szabó–Babuška forms
  $`\lambda_a \lambda_b (\lambda_b - \lambda_a)^j`$, $`j \le k-2`$.
- **Boundary faces:** what remains is fitted with
  $`\lambda_a \lambda_b \lambda_c \, \lambda_b^i \lambda_c^j`$,
  $`i + j \le k-3`$.

A quadratic displacement then stays quadratic inside the element, and the rates
become $`k`$ (last column above). Like gmsh's rule, the new nodes are a fixed
linear combination of the others for each order. Neighbours get the same nodes
on a shared face, to round-off. The elements stay valid: on the coarsest mesh,
the smallest $`\det J`$ relative to the affine element is 0.893 against gmsh's
0.881 at $`k = 3`$, and 0.861 against 0.858 at $`k = 4`$.

A similar construction in our isoparametric solver gave optimal $`H^1`$ rates:
3.36 at $`k = 3`$ and 4.72 at $`k = 4`$.

### References

- P. G. Ciarlet, P.-A. Raviart, "Interpolation theory over curved elements, with
  applications to finite element methods", *Comput. Methods Appl. Mech. Engrg.*
  1 (1972) 217–249.
- M. Lenoir, "Optimal isoparametric finite elements and error estimates for
  domains involving curved boundaries", *SIAM J. Numer. Anal.* 23 (1986) 562–580.
- B. Szabó, I. Babuška, *Finite Element Analysis*, Wiley (1991), chapter on
  mapping and blending functions.
