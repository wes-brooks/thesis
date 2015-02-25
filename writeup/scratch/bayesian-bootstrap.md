---
title: "Dirichlet bootstrap"
author: "Wesley Brooks"
date: "February 25, 2015"
output: html_document
---

#Introduction
Local adaptive grouped regularization (LAGR) is a method for estimating varying coefficient regression (VCR) models [@Brooks-Zhu-Liu-2014]. The method of LAGR simultaneously estimates local regression coefficients and does local variable selection by shrinking some local coefficients to exactly zero. For VCR models estimated by local polynomial approximations to the coefficient functions, pointwise confidence intervals can be estimated from the asymptotic normality of the local coefficient estimates [@Cai-Fan-Li-2000].

Local variable selection in the method of LAGR is done by an adaptive group lasso. Research into the distribution of lasso-type coefficient estimates has been an active subject recently [...] Because the distribution of lasso-type coefficient estimates has a point mass at zero, approximations based on the Gaussian distribution are not appropriate. Resampling-based methods such as the bootstrap are therefore common for lasso-type estimators. In the case of linear regression, the residual bootstrap has been shown to have some nice properties [@Chatterjee-Lahiri-2011].

In this paper we introduce a Dirichlet bootstrap (DB) approach to estimating the distribution of local coefficient estimates. The DB is a specific case of the Bayesian bootstrap of @Rubin-1981 where the Dirichlet prior has all parameters set equal to $1$.

#Methods

##Local adaptive grouped regularization


##Dirichlet bootstrap
The Dirichlet bootstrap is a version of the Bayesian bootstrap... The DB differs from the standard bootstrap in that it generates continuous weights for each observation, rather than resampling each an integer number of times. Thus, while the traditional bootstrap draws weights for the observations from a Multinomial$(1/n, \dots, 1/n)$ distribution, the DB draws the weights from a Dirichlet$(1, \dots, 1)$ distribution.

###Sample dirichlet weights:

Draw $n-1$ iid samples from a uniform distribution on $(0,1)$:

$$u_i \sim U(0,1) \;\;\; {\rm for } \;\;\; i=1,\dots,n-1.$$

Sort $\boldsymbol{u}$ so that $u_{(0)}=0 < u_{(1)} < \cdots < u_{(n-1)} < u_{(n)}=1$. Now the $n$ differences $x_i = u_i - u_{i-1}$ follow a Dirichlet$(1,\dots,1)$ distribution having $n-1$ parameters. Now let $v_i = nx_i$.

###Maximum weighted likelihood
The Bayesian bootstrap proceeds by maximizing the weighted likelihood:

$$L(s) = \prod_{i=1}^n \exp[ v_i \; K_h(\|s_i-s\|) \; \sigma^{-2} \{ y_i - x_i\beta(s) \}^2 ]$$

The observations weights are the product of kernel weights $K_h(\cdot)$ and the Dirichlet weights $\boldsymbol{v}$. The weighted likelihood as written is for the linear regression setting, and can be maximized via weighted least squares. Thus, one realization of the BB estimate is:

$$\hat{\beta}^*(s) = (X^T W V^* X)^{-1} X^T W V^* Y$$

where $W$ is diagonal with $W_{ii} = K_h(\|s_i-s\|)$ and $V^*$ is diagonal with $V^*_{ii} = v_i$. As many samples can be drawn from the BB distribution as necessary.

####Proposition
That the BB samples, $\hat{\beta}^*(s)$, are drawn from the same distribution as the MLE $\hat{\beta}(s)$.


#Simulation
