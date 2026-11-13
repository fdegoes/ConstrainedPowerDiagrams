# Unified Variational Formulation for Shifted Power Diagrams

This document summarizes the unified objective function and the derivatives required to compute a Shifted Power Diagram $(x', u', w')$ that is invariant to affine transformations of the domain boundary.

---

## I. The Unified Objective Function
The objective $\mathcal{F}(x', u', w')$ is defined using the constant stiffness weights $\omega_{ij} = A_{ij}/l_{ij}$ and interpolation constants $\tau_{ij}$ from the reference diagram.

$$\mathcal{F} = \sum_{i \sim j} \omega_{ij} \langle x'_i, x'_j + u'_j \rangle + \sum_i w'_i \left( \sum_{j \sim i} \omega_{ij} b_{ij}(x', u') \right)$$

Where the geometric power-distance residual $b_{ij}$ is:
$$b_{ij} = \|x'_j + u'_j\|^2 - \|x'_i + u'_i\|^2 + 2\langle z'_{ij}, (x'_i + u'_i) - (x'_j + u'_j) \rangle$$

*   **Boundary Point Interpolation:** $z'_{ij} = (1-\tau_{ij})x'_i + \tau_{ij}x'_j$.
*   **Dual Site:** $Y_i = x'_i + u'_i$.

---

## II. First Derivatives (The Residuals $R$)
The stationary point $\nabla \mathcal{F} = 0$ recovers the following three equilibrium conditions:

1. **Weight Gradient ($\nabla_{w'}$):** 
   $$R_{w_i} = \sum_{j \sim i} \omega_{ij} b_{ij}(x', u') = 0$$
   *Ensures power-distance consistency at the cell boundaries.*

2. **Shift Gradient ($\nabla_{u'}$):** 
   $$R_{u_i} = \sum_{j \sim i} \omega_{ij} x'_j + 2 \sum_{j \sim i} \omega_{ij} (w'_i - w'_j) (z'_{ij} - (x'_i + u'_i)) = 0$$
   *Forces the shift $u'$ to satisfy the harmonic metric correction.*

3. **Site Gradient ($\nabla_{x'}$):** 
   $$R_{x_i} = \sum_{j \sim i} \omega_{ij} (x'_j + u'_j) + 2 \sum_{j \sim i} \omega_{ij} (w'_i - w'_j) (1-\tau_{ij})((x'_i + u'_i) - (x'_j + u'_j)) = 0$$
   *Forces the sites $x'$ to satisfy the harmonic site identity.*

---

## III. Second Derivatives (The Jacobian Blocks)
To solve via Newton-Raphson, use the $3 \times 3$ symmetric block Jacobian (Hessian of $\mathcal{F}$):

$$
\mathbf{J} = \begin{bmatrix}
\mathbf{H}_{xx} & \mathbf{H}_{xu} & \mathbf{J}_{xw} \\
\mathbf{H}_{ux} & \mathbf{H}_{uu} & \mathbf{J}_{uw} \\
\mathbf{J}_{wx} & \mathbf{J}_{wu} & \mathbf{0}
\end{bmatrix}
$$

*   **Linearized Geometry ($\mathbf{J}_{wu}$):** $2 \omega_{ij} (z'_{ij} - (x'_j + u'_j))^T$
*   **Linearized Geometry ($\mathbf{J}_{xw}$):** $2 \omega_{ij} (1-\tau_{ij})((x'_i + u'_i) - (x'_j + u'_j))^T$
*   **Shift Stiffness ($\mathbf{H}_{uu}$):** $-2 \sum_{j \sim i} \omega_{ij} (w'_i - w'_j) \mathbf{I}$
*   **Coupling Block ($\mathbf{H}_{xu}$):** $\omega_{ij} \mathbf{I} - 2 \omega_{ij} (w'_i - w'_j) (1-\tau_{ij}) \mathbf{I}$

---

## IV. Special Case: Affine Transformation ($x' = Ax + t$)
If the boundary is transformed affinely, the stationary point matches the following analytical expressions:


| Parameter | Expression | Description |
| :--- | :--- | :--- |
| **Site** | $x'_i = A P_i + t$ | Affine map of original sites. |
| **Shift** | $u'_i = ((A^T)^{-1} - A) P_i + \text{const}$ | Corrects for non-isometric shear/scale. |
| **Weight** | $w'_i = w_i + P_i^T ((A^T A)^{-1} - I) P_i + \dots$ | Updates for metric change. |

---

## V. Proof of Affine Invariance
**Property:** Under affine transformation $A$, the inner product $\langle A E_{ij}, (A^T)^{-1} \Delta P \rangle = E_{ij}^T \Delta P = 0$ is preserved. 

Because $x'$ is driven by the boundary Laplacian and $u'$ is driven by the metric-harmonic Laplacian, the resulting shifted power diagram's geometry is the **exact affine image** of the original.
