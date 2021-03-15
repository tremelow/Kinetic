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

As a preliminary to treating the Boltzmann equation, we consider
hyperbolic problems of the form

<p align="center"><img src="https://latex.codecogs.com/svg.latex?%5Cbegin%7Bcases%7D%20%5Cdisplaystyle%20%5Cpartial_t%20u%20&plus;%20%5Cpartial_x%20v%20%3D%200%20%2C%20%5C%5C%20%5Cdisplaystyle%20%5Cpartial_t%20v%20&plus;%20%5Cfrac%7B1%7D%7B%5Cvarepsilon%5E%7B2%5Calpha%7D%7D%20%5Cpartial_x%20p%28u%29%20%3D%20-%5Cfrac%7B1%7D%7B%5Cvarepsilon%5E%7B1&plus;%5Calpha%7D%7D%20%5Cleft%28%20v%20-%20f%28u%29%20%5Cright%29%20.%20%5Cend%7Bcases%7D" alt="Stiff hyperbolic relaxation problem"/>
<!--
$$
\left\{ \begin{array}
    \partial_t u + \partial_x v = 0 \\
    \partial_t v + \frac{1}{\varepsilon^{2\alpha}} \partial_x p(u) ,
    = -\frac{1}{\varepsilon^{1+\alpha}} \left( v - f(u) \right) .
\end{array} \right.
$$
-->
</p>

Note that for the linearized Boltzmann equation, the particle density 
![f(t,x,v)](https://latex.codecogs.com/svg.latex?%5Cinline%20f%28t%2Cx%2Cv%29)
satisfies
<p align="center"><img src="https://latex.codecogs.com/svg.latex?%5Cvarepsilon%20%5Cpartial_t%20f%20&plus;%20v%20%5Ccdot%20%5Cnabla_x%20f%20%3D%20%5Cfrac%7B1%7D%7B%5Cvarepsilon%7D%20Lf" alt="Boltzmann equation"/>
<!--
$$
    \varepsilon \partial_t f + v \cdot \nabla_x f 
    = \frac{1}{\varepsilon} Lf
$$
-->
</p>

which in the case ![Assump. x, v](https://latex.codecogs.com/svg.latex?%5Cinline%20x%20%5Cin%20%5COmega%20%5Csubseteq%20%5Cmathbb%7BR%7D%2C%5C%20v%20%5Cin%20%5C%7B-1%2C%201%5C%7D)
can be written in the form of the telegraph equation. Setting the mass
![Def. rho](https://latex.codecogs.com/svg.latex?%5Cinline%5Crho%20%3D%20f%281%29%20&plus;%20f%28-1%29)
and the current ![Def.
j](https://latex.codecogs.com/svg.latex?%5Cinline%5Cvarepsilon%20j%20%3D%20f%281%29%20-%20f%28-1%29),
<p align="center"><img src="https://latex.codecogs.com/svg.latex?%5Cinline%20%5Cbegin%7Bcases%7D%20%5Cpartial_t%20%5Crho%20&plus;%20%5Cpartial_x%20j%20%3D%200%20%2C%20%5C%5C%20%5Cdisplaystyle%20%5Cpartial_t%20j%20&plus;%20%5Cfrac%7B1%7D%7B%5Cvarepsilon%5E2%7D%20%5Cpartial_x%20%5Crho%20%3D%20-%5Cfrac%7B1%7D%7B%5Cvarepsilon%5E2%7D%20j%20%2C%20%5Cend%7Bcases%7D" alt="Telegraph equation"/>
<!--
$$
\begin{cases}
    \partial_t \rho + \partial_x j = 0 , \\
    \partial_t j + \frac{1}{\varepsilon^2} \partial_x \rho
    = -\frac{1}{\varepsilon^2} j ,
\end{cases}
$$
-->
</p>
which fits with the first hyperbolic problem.


This will eventually be implemented in context of simulations using the 
[DrWatson library](https://juliadynamics.github.io/DrWatson.jl/dev/).