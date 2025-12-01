# 📘 Pivotal Based Inference

**Authors:**  
*Kundan Singh¹, Yogesh Mani Tripathi², and Liang Wang³*

---

## 📌 Shape–Scale Family (SSF) Distribution

Following Maswadah (2022), the probability density function (PDF) of the **shape–scale family (SSF)** is:

$$
f(x;\alpha,\beta)
= \alpha\beta\,[g(x)]^{\alpha-1}\, g'(x)\,
\exp\left\{-\beta [g(x)]^\alpha\right\},
\qquad x>0,\ \alpha,\beta>0.
$$

where:

- \(g(x)\) is differentiable and increasing,  
- \(g(0^+)=0\),  
- \(g(x)\to\infty\) as \(x\to\infty\).

The cumulative distribution function (CDF) is:

$$
F(x;\alpha,\beta)=1 - \exp\{-\beta [g(x)]^\alpha\}.
$$

The survival function (SF), hazard rate function (HRF), and mean time to failure (MTTF) are:

$$
S(x;\alpha,\beta)=\exp\{-\beta[g(x)]^\alpha\},
$$

$$
H(x;\alpha,\beta)=\alpha\beta[g(x)]^{\alpha-1}g'(x),
$$

$$
\text{MTTF}(\alpha,\beta)=\int_0^\infty x f(x;\alpha,\beta)\,dx.
$$

➡️ Many well-known lifetime models (Weibull, Pareto, Lomax, Burr-XII, exponential, etc.) are special cases depending on the function \(g(x)\).

---

# 📌 Pivotal Quantities for Generalized Inference

To construct generalized point estimators (GPEs) and generalized confidence intervals (GCIs) under **block progressive censoring**, the following pivotal quantities are used.

---

## **Theorem 1**

For \(i = 1,2,\ldots,m\), define:

$$
P^{\mathcal{X}_i}(\alpha)
= 2 \sum_{j=1}^{s_i-1}
\ln\left[
\frac{
\sum_{r=1}^{s_i-1}(R_{ir}+1)[g(X_{ir})]^\alpha +
\left(n_i-\sum_{r=1}^{s_i-1}(R_{ir}+1)\right)[g(X_{is_i})]^\alpha
}{
\sum_{r=1}^{j-1}(R_{ir}+1)[g(X_{ir})]^\alpha +
\left(n_i-\sum_{r=1}^{j-1}(R_{ir}+1)\right)[g(X_{ij})]^\alpha
}
\right].
$$

and

$$
S^{\mathcal{X}_i}(\alpha,\beta_i)
=2\beta_i\left[
\sum_{r=1}^{s_i-1}(R_{ir}+1)[g(X_{ir})]^\alpha +
\left(n_i-\sum_{r=1}^{s_i-1}(R_{ir}+1)\right)[g(X_{is_i})]^\alpha
\right].
$$

Then:

- \(P^{\mathcal{X}_i}(\alpha) \sim \chi^2_{2(s_i-1)}\),  
- \(S^{\mathcal{X}_i}(\alpha,\beta_i) \sim \chi^2_{2s_i}\),  
- They are statistically **independent**.

---

## **Theorem 2**

For each \(i = 1,2,\dots,m\), the pivotal quantity  
\(P^{\mathcal{X}_i}(\alpha)\) is **increasing in \(\alpha\)**.

---

# 📌 Combined Pivotal Quantity

Define:

$$
P(\alpha)=\sum_{i=1}^m P^{\mathcal{X}_i}(\alpha).
$$

Since the terms are independent:

$$
P(\alpha) \sim \chi^2_{2\sum_{i=1}^{m}(s_i-1)}.
$$

Because \(P(\alpha)\) is monotone increasing in \(\alpha\), the equation

$$
P(\alpha)=\alpha^*, \qquad 
\alpha^* \sim \chi^2_{2\sum_{i=1}^{m}(s_i-1)}
$$

has a **unique solution**:

$$
\alpha = h(\alpha^*,T).
$$

For each \(i\):

$$
\beta_i=\frac{S^{\mathcal{X}_i}}{G^{\mathcal{X}_i}(\alpha)},
$$

where:

$$
G^{\mathcal{X}_i}(\alpha)
=2\left[
\sum_{r=1}^{s_i-1}(R_{ir}+1)[g(X_{ir})]^\alpha +
\left(n_i-\sum_{r=1}^{s_i-1}(R_{ir}+1)\right)[g(X_{is_i})]^\alpha
\right].
$$

Thus the generalized pivotal quantity for \(\beta_i\) is:

$$
\frac{S^{\mathcal{X}_i}}{G^{\mathcal{X}_i}(h(\alpha^*,T))}.
$$

---

# 📌 Algorithm for Generalized Estimation

### **Algorithm**

1. **Initialize** \(t = 1\).

2. **Simulate \(\alpha\):**

   - Draw \(\alpha^* \sim \chi^2_{2\sum_{i=1}^m(s_i-1)}\).  
   - Solve \(P(\alpha)=\alpha^*\) to obtain \(\alpha^{(t)}\).

3. **Simulate \(\beta_i\)**:

   - Generate \(S^{\mathcal{X}_i} \sim \chi^2_{2s_i}\).  
   - Compute:

     $$
     \beta_i^{(t)} =
     \frac{S^{\mathcal{X}_i}}{G^{\mathcal{X}_i}(\alpha^{(t)})}.
     $$

4. **Combine \(\beta_i\)** using inverse-variance weighting:

$$
\beta^{(t)}
= \frac{
\sum_{i=1}^m \beta_i^{(t)} / \mathrm{Var}(\beta_i)
}{
\sum_{i=1}^m 1 / \mathrm{Var}(\beta_i)
}.
$$

5. **Compute reliability measures**:

- \(S(x_0;\alpha^{(t)},\beta^{(t)})\)  
- \(H(x_0;\alpha^{(t)},\beta^{(t)})\)  
- \(\text{MTTF}(\alpha^{(t)},\beta^{(t)})\)

6. **Repeat Steps 2–4** for \(N\) iterations to obtain:

$$
\zeta^{(1)}, \dots, \zeta^{(N)},
\qquad 
\zeta=\{\alpha,\beta,\mathbf{B}, S(x_0), H(x_0), \text{MTTF}\}.
$$

7. **Generalized Confidence Interval (GCI)**:

   - Sort the values:

     $$
     \zeta^{[1]} \le \cdots \le \zeta^{[N]}.
     $$

   - For confidence level \(1-\xi\):

     $$
     (\zeta^{[j]},\ \zeta^{[j + N - (\lfloor N\xi\rfloor+1)]}).
     $$

8. **Generalized Point Estimator**:

$$
\widehat{\zeta} = \frac{1}{N}\sum_{j=1}^{N} \zeta^{(j)}.
$$

---

# 📚 References

- Maswadah (2022). *Improved maximum likelihood estimation of the shape-scale family...*  
  **Journal of Applied Statistics**, 49(11), 2825–2844.

- Weerahandi (2004). *Generalized Inference in Repeated Measures.* Wiley.

- Singh et al. (2023). *Inference for reliability in a multicomponent stress–strength model...*  
  **Quality Technology & Quantitative Management**, 21(2), 147–176.

- Singh et al. (2024).  
  *Inference for unit inverse Weibull distribution under block progressive Type-II censoring.*  
  **Journal of Statistical Theory and Practice**, 18(3), 42.
