Sample dirichlet weights:

Draw $n-1$ iid samples from a uniform distribution on $(0,1)$:

$$u_i \sim U(0,1) \;\;\; {\rm for } \;\;\; i=1,\dots,n-1.$$

Sort $\boldsymbol{u}$ so that $u_{(0)}=0 < u_{(1)} < \cdots < u_{(n-1)} < u_{(n)}=1$. Now the $n$ differences $x_i = u_i - u_{i-1}$ follow a Dirichlet$(0,\dots,0)$ distribution having $n-1$ parameters. Now let $v_i = nx_i$.

Now the Bayesian bootstrap proceeds by maximizing the weighted likelihood:

$$L(s) = \prod_{i=1}^n \exp[ v_i \; K_h(\|s_i-s\|) \; \sigma^{-2} \{ y_i - x_i\beta(s) \}^2 ]$$

The observations weights are the product of kernel weights $K_h(\cdot)$ and the Dirichlet weights $\boldsymbol{v}$. The weighted likelihood as written is for the linear regression setting, and can be maximized via weighted least squares. Thus, one realization of the BB estimate is:

$$\hat{\beta}^*(s) = (X^T W V^* X)^{-1} X^T W V^* Y$$

where $W$ is diagonal with $W_{ii} = K_h(\|s_i-s\|)$ and $V^*$ is diagonal with $V^*_{ii} = v_i$. As many samples can be drawn from the BB distribution as necessary.

