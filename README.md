<!-- MathJax -->
<script src="https://polyfill.io/v3/polyfill.min.js?features=es6"></script>
<script id="MathJax-script" async src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script>



# StiffKinetic.jl

![Lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)<!--
![Lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-stable-green.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-retired-orange.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-archived-red.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-dormant-blue.svg) -->
[![Build Status](https://travis-ci.com/tremelow/Kinetic.jl.svg?branch=master)](https://travis-ci.com/tremelow/Kinetic.jl)
[![codecov.io](http://codecov.io/github/tremelow/Kinetic.jl/coverage.svg?branch=master)](http://codecov.io/github/tremelow/Kinetic.jl?branch=master)

A Julia library to provide a structure and high-performance numerical
tools for the simulation of kinetic problems with stiff collisions. This
is my first library ever, so it won't be world changing!

This will eventually be implemented in context of simulations using the 
[DrWatson library](https://juliadynamics.github.io/DrWatson.jl/dev/).

## Mathematical context

As a preliminary to treating the Boltzmann equation, we consider
hyperbolic problems of the form

<p align="center">$$
\begin{cases}
    \partial_t u + \partial_x v = 0 , \\ \displaystyle
    \partial_t v + \frac{1}{\varepsilon^{2\alpha}} \partial_x p(u)
    = -\frac{1}{\varepsilon^{1+\alpha}} \left( v - f(u) \right) .
\end{cases}
$$</p>

Note that for the linearized Boltzmann equation, the particle density 
![f(t,x,v)](https://latex.codecogs.com/svg.latex?%5Cinline%20f%28t%2Cx%2Cv%29)
satisfies
<p align="center">$$
    \varepsilon \partial_t f + v \cdot \nabla_x f 
    = \frac{1}{\varepsilon} Lf
$$</p>

which in the case ![Assump. x, v](https://latex.codecogs.com/svg.latex?%5Cinline%20x%20%5Cin%20%5COmega%20%5Csubseteq%20%5Cmathbb%7BR%7D%2C%5C%20v%20%5Cin%20%5C%7B-1%2C%201%5C%7D)
can be written in the form of the telegraph equation. Setting the mass
![Def. rho](https://latex.codecogs.com/svg.latex?%5Cinline%20%5Crho%20%3D%20f%281%29%20&plus;%20f%28-1%29)
and the current ![Def.
j](https://latex.codecogs.com/svg.latex?%5Cinline%5Cvarepsilon%20j%20%3D%20f%281%29%20-%20f%28-1%29),
<p align="center">$$
\begin{cases} \displaystyle
    \partial_t \rho + \partial_x j = 0 , \\ \displaystyle
    \partial_t j + \frac{1}{\varepsilon^2} \partial_x \rho
    = -\frac{1}{\varepsilon^2} j ,
\end{cases}
$$</p>

which fits with the first hyperbolic problem.
