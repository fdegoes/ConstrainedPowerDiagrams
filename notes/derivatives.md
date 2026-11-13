# Derivation of Gradients for the Unified Integrated Ratio-OT Formulation

This document provides the step-by-step derivation of the residuals (first derivatives) for the unified objective function $\mathcal{F}$.

---

## 1. The Objective Function
The energy $\mathcal{F}(x', Y, w')$ is defined as:
$$\mathcal{F} = \alpha \sum_{i \sim j} \omega_{ij} \langle x'_i - x'_j, Y_i - Y_j \rangle + \sum_i \int_{\Omega'_i} \left( \|x - Y_i\|^2 - w'_i \frac{V'_i}{V_i} \right) \, dx$$

Where:
*   $V'_i = \sum_{j \sim i} \omega_{ij} \left( \frac{1}{2}\|Y_i - Y_j\|^2 + w'_i - w'_j \right)$
*   $\Omega'_i$ is the Power Cell defined by the effective weight $W^{eff}_i = w'_i \frac{V'_i}{V_i}$.

---

## 2. Gradient with respect to Weights ($\nabla_{w'_k}$)

Using the **Reynolds Transport Theorem**, the derivative of the integral over the partition $\Omega'$ with respect to the weights vanishes at the boundaries because the power distance is continuous across the edges. Thus, we differentiate only the integrand.

**Step A: Partial Derivative of the Integrand**
The weight $w'_k$ appears in the integral of cell $k$ and in the $V'_i$ terms of all cells $i$ adjacent to $k$.
$$\frac{\partial \mathcal{F}}{\partial w'_k} = \int_{\Omega'_k} \left( -\frac{V'_k}{V_k} \right) dx - \sum_i \int_{\Omega'_i} \frac{w'_i}{V_i} \frac{\partial V'_i}{\partial w'_k} dx$$

**Step B: Differentiating the Linearized Volume $V'$**
*   For $i = k$: $\frac{\partial V'_k}{\partial w'_k} = \sum_{j \sim k} \omega_{kj}$
*   For $i \in \text{neighbors}(k)$: $\frac{\partial V'_i}{\partial w'_k} = -\omega_{ik}$

**Step C: Substitution and Assembly**
$$R_{w_k} = -|\Omega'_k| \frac{V'_k}{V_k} - \left( \frac{w'_k |\Omega'_k|}{V_k} \sum_{j \sim k} \omega_{kj} + \sum_{j \sim k} \frac{w'_j |\Omega'_j|}{V_j} (-\omega_{jk}) \right)$$
Rearranging the terms:
$$R_{w_k} = -|\Omega'_k| \frac{V'_k}{V_k} + \sum_{j \sim k} \omega_{kj} \left( |\Omega'_j| \frac{w'_j}{V_j} - |\Omega'_k| \frac{w'_k}{V_k} \right) = 0$$

---

## 3. Gradient with respect to Dual Sites ($\nabla_{Y_k}$)

**Step A: Dirichlet Term**
$$\frac{\partial}{\partial Y_k} \alpha \sum_{i \sim j} \omega_{ij} \langle x'_i - x'_j, Y_i - Y_j \rangle = \alpha \sum_{j \sim k} \omega_{kj} (x'_k - x'_j)$$

**Step B: Optimal Transport Term**
Differentiating the squared distance integral:
$$\int_{\Omega'_k} \frac{\partial}{\partial Y_k} \|x - Y_k\|^2 dx = \int_{\Omega'_k} -2(x - Y_k) dx = -2 |\Omega'_k| (c_k - Y_k) = 2 |\Omega'_k| (Y_k - c_k)$$

**Step C: Volume Ratio Term**
Differentiating the $V'_i$ dependency:
$$\frac{\partial}{\partial Y_k} \sum_i \int_{\Omega'_i} -w'_i \frac{V'_i}{V_i} dx = - \sum_i \frac{w'_i |\Omega'_i|}{V_i} \frac{\partial V'_i}{\partial Y_k}$$
Using $\frac{\partial V'_i}{\partial Y_k} = \omega_{ik}(Y_k - Y_i)$ and assuming $|\Omega'|$ is locally consistent:
$$\text{Term} = - |\Omega'_k| \sum_{j \sim k} \omega_{kj} \left( \frac{w'_k}{V_k} + \frac{w'_j}{V_j} \right) (Y_k - Y_j)$$

**Step D: Final Assembly for $R_{Y_k}$**
$$R_{Y_k} = \alpha \mathbf{L}x' + 2 |\Omega'_k| (Y_k - c_k) - |\Omega'_k| \sum_{j \sim k} \omega_{kj} \left( \frac{w'_k}{V_k} + \frac{w'_j}{V_j} \right) (Y_k - Y_j) = 0$$

---

## 4. Gradient with respect to Site Positions ($\nabla_{x'_k}$)

**Step A: Dirichlet Term**
$$\frac{\partial}{\partial x'_k} \alpha \sum_{i \sim j} \omega_{ij} \langle x'_i - x'_j, Y_i - Y_j \rangle = \alpha \sum_{j \sim k} \omega_{kj} (Y_k - Y_j)$$

**Step B: Integral Terms**
Since $x'_k$ only appears inside the $Y$ and $W$ terms of the integral (which were already differentiated in the $Y$ derivation), the integral gradients are identical to the $Y$ derivation.

**Step C: Final Assembly for $R_{x_k}$**
$$R_{x_k} = \alpha \mathbf{L}Y + 2 |\Omega'_k| (Y_k - c_k) - |\Omega'_k| \sum_{j \sim k} \omega_{kj} \left( \frac{w'_k}{V_k} + \frac{w'_j}{V_j} \right) (Y_k - Y_j) = 0$$
