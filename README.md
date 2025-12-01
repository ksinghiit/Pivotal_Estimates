\title{Pivotal based inference}
 \author{Kundan Singh$^1$, Yogesh Mani Tripathi$^2$  and Liang Wang$^{3}$} 

\maketitle \noindent 
 
{Following \citet{Maswadah2022}}, the probability density function (PDF) of  a shape-scale family (SSF) distribution is given by
\begin{align}\label{eq1}
    f(x;\alpha,\beta)=\alpha\beta (g(x))^{\alpha-1} g^{'}(x) \exp{\left\{-\beta (g(x))^{\alpha}\right\}},\ x>0. \alpha,\beta>0,
\end{align}
{where $g^{'}(x)=\frac{dg(x)}{dx}$, and function $g(x)$ is differentiable and increases in $x$ with} $g(0^+)=0$ and $g(x) \to +\infty$ as $x \to +\infty$. The cumulative distribution function (CDF) is 
\begin{align}\label{eq2}
    F(x;\alpha,\beta)=1-\exp{\left\{-\beta (g(x))^{\alpha}\right\}},\ \alpha, \beta>0, x>0,
\end{align}
where $\alpha>0$ and $\beta>0$ are shape and scale  {parameters, respectively. For simplicity, this model is denoted as SSF$(\alpha,\beta)$}. 
Also, the survival function (RF) $S(\cdot)$, hazard rate function (HRF) $H(\cdot)$ and the mean time to failure (MDTF) of this family are expressed as {
\begin{align}\label{eq3}
\aligned
S(x;\alpha,\beta)&= \exp{\left\{-\beta (g(x))^{\alpha}\right\}},\\
H(x;\alpha,\beta)&=\alpha\beta (g(x))^{\alpha-1} g^{'}(x), \\
   MDTF(\alpha,\beta)&=\int_X xf(x;\alpha,\beta)dx.
\endaligned
\end{align}}
 {For clarity, is is noted that some well-known} lifetime distributions like exponential, Weibull, modified Weibull, Pareto, generalized Pareto, Lomax, Burr-type-XII, and unit inverse Weibull belong to this family according to  different expressions of $(g(x))^\alpha$.
 \\
 For constructing generalized point estimators (GPEs) and generalized confidence intervals (GCIs) for model parameters and reliability indices under block progressive censoring, the following pivotal quantities are proposed first

\begin{theorem}\label{thm4.1}
For $i=1,2,\ldots,m$, quantities 
\begin{align*}
    P^{\mathcal{X}_i}(\alpha)=2 \sum_{j=1}^{s_i-1}\ln\left[\frac{\sum_{r=1}^{s_i-1}(R_{ir}+1)(g(X_{ir}))^\alpha+\left(n_i-\sum_{r=1}^{s_i-1}(R_{ir}+1)\right)(g(X_{is_i}))^\alpha}{\sum_{r=1}^{j-1}(R_{ir}+1)(g(X_{ir}))^\alpha+\left(n_i-\sum_{r=1}^{j-1}(R_{ir}+1)\right)(g(X_{ij}))^\alpha}\right],
\end{align*}
and
\begin{align*}
S^{\mathcal{X}_i}(\alpha,\beta_i)=2\beta_i\left[\sum_{r=1}^{s_i-1}(R_{ir}+1)(g(X_{ir}))^\alpha+\left(n_i-\sum_{r=1}^{s_i-1}(R_{ir}+1)\right)(g(X_{is_i}))^\alpha\right],
\end{align*}
have chi-square distributions with $2(s_i-1)$ and $2s_i$ degrees of freedom respectively, and they are statistically independent.
\end{theorem}

\begin{theorem}\label{thm4.2}
    For each $i=1,2,\dots,m$, the pivotal quantity $P^{\mathcal{X}_i}(\alpha)$  {increases} in $\alpha$.
\end{theorem}

Now, consider the pivotal quantity
\begin{align}
    P(\alpha)=&\sum_{i=1}^{m}P^{\mathcal{X}_i}(\alpha)\nonumber\\
    =& ~2 \sum_{i=1}^{m}\sum_{j=1}^{s_i-1}\ln\left[\frac{\sum_{r=1}^{s_i-1}(R_{ir}+1)(g(X_{ir}))^\alpha+\left(n_i-\sum_{r=1}^{s_i-1}(R_{ir}+1)\right)(g(X_{is_i}))^\alpha}{\sum_{r=1}^{j-1}(R_{ir}+1)(g(X_{ir}))^\alpha+\left(n_i-\sum_{r=1}^{j-1}(R_{ir}+1)\right)(g(X_{ij}))^\alpha}\right].    
\end{align}
Since $P^{\mathcal{X}_i}(\alpha)$ are independent from one another for $i=1,2,\dots,m$. So $P(\alpha)$ has a chi-square distribution with degree of freedom $2\sum_{i=1}^{m}(s_i-1)$. Also,  {using Theorem \ref{thm4.2}, it is observed that quantity $P(\alpha)$ increases in $\alpha$}.

Note that $P(\alpha)=\alpha^*$ has a unique solution for $\alpha^*\sim \chi^2_{2\sum_{i=1}^{m}(s_i-1)}$, denoted as $h(\alpha^*,T)$. Moreover, $\beta_i$ can be written as
\begin{align}
   \beta_i=\frac{S^{\mathcal{X}_i}}{G^{\mathcal{X}_i}(\alpha)}, i=1,2,\dots,m,
\end{align}
where $G^{\mathcal{X}_i}(\alpha)=2\Big[\sum_{r=1}^{s_i-1}(R_{ir}+1)(g(X_{ir}))^\alpha+\left(n_i-\sum_{r=1}^{s_i-1}(R_{ir}+1)\right)(g(X_{is_i}))^\alpha\Big]$. Then the generalized pivotal quantity for $\beta_i,\ i=1,2\dots,m$, can be obtained as $\frac{S^{\mathcal{X}_i}}{G^{\mathcal{X}_i}(g(\alpha^{*},T))}$ (see, \citet{Weerahandi2004}). Based on these generalized pivotal quantities, an algorithm is proposed next, which is useful for obtaining estimates for reliability and DDF characteristics.

\textbf{Algorithm:}
\begin{itemize}[noitemsep]
    \item[Step 1.] Set $t=1$.
    \item[Step 2.] Generate $\alpha^*$ from $\chi^2_{2\sum_{i=1}^{m}(s_i-1)}$ distribution, and obtained $\alpha^{(t)}$ by solving $P(\alpha)=\alpha^*$.
    \item[Step 3.] For each $i=1,2.\dots,m$, generate $S^{\mathcal{X}_i}$ form $\chi^2_{2s_i}$ distribution and get an observation $\beta_i^{(t)}$ for $\beta_i$ using $\beta_i^{(t)}=\frac{S^{\mathcal{X}_i}}{G^{\mathcal{X}_i}(\alpha^t)}$.
    \item[Step 4.] Compute $\beta^{(t)}=\frac{\sum_{i=1}^{m}\frac{1}{Var(\beta_i)}\beta_i^t}{\sum_{i=1}^{m}\frac{1}{Var(\beta_i)}}$, where $Var(\beta_i)=\frac{1}{t}\sum_{r=1}^{t}\left(\beta_i^{(r)}-\frac{1}{t}\sum_{r=1}^{t}\beta_i^{(r)}\right)^2$. Obtained the reliability indices $S(x_0;\alpha^{(t)},\beta^{(t)}),H(x_0;\alpha^{(t)},\beta^{(t)})$ and MDTF$(\alpha^{(t)},\beta^{(t)})$.
    
  \item[Step 5.] Repeat Steps $2$ and $4,$ $\mathcal{N}$ times, one can obtain $\mathcal{N}$ values of $\zeta=\{\alpha,\mathbf{B}, \beta, S(x_0), H(x_0), MDTF\}$  as $\zeta^{(1)},\zeta^{(2)},\dots,\zeta^{(\mathcal{N})}$.
 
\item[Step 6.] To construct the generalized confidence interval (GCI) of $\zeta$, first arrange all estimates of $\zeta$ in an ascending order as $\zeta^{[1]},\zeta^{[2]},\dots,\zeta^{[N]}$. For arbitrary $0<\xi<1$, a $100(1-\xi)\%$ confidence interval of $\zeta$ can be obtained as
\begin{eqnarray}
\Big(\zeta^{[j]},\zeta^{[j+\mathcal{N}-[\mathcal{N}_{\xi}+1]]}\Big), j=1,2,\dots,[\mathcal{N}_{\xi}].\nonumber
\end{eqnarray}
where $[r]$ denotes the greatest integer less than or equal to $r$. Therefore, the $100(1-\xi)\%$ GCI of $\zeta$ can be constructed as the $j^*$th one satisfying
\begin{eqnarray}
\zeta^{[j+\mathcal{N}-[\mathcal{N}_{\xi}+1]]} -\zeta^{[j^*]}=\text{min}\big(\zeta^{[j+\mathcal{N}-[\mathcal{N}_{\xi}+1]]},\zeta^{[j]} \big).\nonumber
\end{eqnarray}
In addition, based on $\zeta^{(1)},\zeta^{(1)},\dots,\zeta^{(\mathcal{N})}$, a generalized point estimator for $\zeta$ is given by
\begin{eqnarray}
\widehat{\zeta} = \frac{1}{\mathcal{N}}\sum_{j=1}^{\mathcal{N}}\zeta^{(j)}.\nonumber
\end{eqnarray}
\end{itemize}


\begin{thebibliography}{Yogesh}

\bibitem[Maswadah(2022)]{Maswadah2022} Maswadah, M. (2022). Improved maximum likelihood estimation of the shape-scale family based on the generalized progressive hybrid censoring scheme. {\it Journal of Applied Statistics}, 49(11), 2825-2844.

 
\bibitem[Weerahandi(2004)]{Weerahandi2004} Weerahandi, S. (2004). {\it Generalized inference in repeated measures: Exact methods in MANOVA and mixed models} (Vol. 500). John Wiley \& Sons.

\bibitem[Singh et al.(2023)]{Singh2023} Singh, K., Mahto, A. K., Tripathi, Y., \& Wang, L. (2023). Inference for reliability in a multicomponent stress-strength model for a unit inverse Weibull distribution under  {Type-II} censoring. {\it Quality Technology \& Quantitative Management}, 21(2), 147-176.

\bibitem[Singh et al.(2025)]{Singh2025}Singh, K., Tripathi, Y. M., Lodhi, C., \& Wang, L. (2024). Inference for unit inverse Weibull distribution under block progressive type-ii censoring. {\it Journal of Statistical Theory and Practice}, 18(3), 42.
\end{thebibliography}


\end{document}
