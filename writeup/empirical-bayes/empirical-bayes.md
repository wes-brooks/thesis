---
title: "An Empirical Bayes procedure for inference in local polynomial models"
author: "Wesley Brooks"
date: "February 25, 2015"
output: html_document
---

#Introduction
Local adaptive grouped regularization (LAGR) is a method for estimating varying coefficient regression (VCR) models [@Brooks-Zhu-Liu-2014]. The method of LAGR utilizes local polynomial regression, a form of kernel smoothing [@Fan-Gijbels-1996]. Kernel smoothing requires estimation of a bandwidth parameter $\hat{h}$, and then estimation and inference proceed conditional on $\hat{h}$.

However, there is uncertainty in the estimation of $\hat{h}$, which is not reflected when estimates are conditional on $\hat{h}$. Furthermore, the bandwidth is a nuisance parameter, so we would prefer to report results marginal to $h$ rather than conditional on $\hat{h}$.

On the other hand, we want to avoid making any assumptions about $h$ because such assumptions would be difficult to justify, difficult to check, and could affect estimation and inference for the parameters of interest. 

#Methods
