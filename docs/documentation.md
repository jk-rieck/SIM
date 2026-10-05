# Introduction
--------------
 
This is the documentation for the McGill Sea Ice Model.

## Masks 
--------

The model uses a mask filled with integers to identify land from ocean. Ocean is marked with 1, and land is marked with 0. 

There are two masks that are used in the model: `maskC` and `maskB`. `maskC` is a $(n_{x} + 2) \times (n_{y} + 2)$ array that indicates the presence of land at the cell center. `maskB` is called the velocity mask and has one more element in each dimension.

![masks](source/mask.png "Masks")
  
      In the figure, dots stand for the ocean/land mask. The `x` stand for the `maskB`, also called the velocity mask. An element in `maskB` is set to 1 if all the surrounding points in `maskC` are equal to 1. Hence, all points on the borders of the grid are equal to 0. 


Eventually, it would be cleaner to use a logical mask instead of an integer mask. This would simplify some expressions and the code could then use the instrinsic functions in Fortran 90 that take logical mask as arguments. 

## Finite Differences
---------------------
The standard way to compute finite differences (numerical derivatives) is to use Taylors series expansion formula at different locations:

$
\begin{align}
f(x + \delta) & \approx f(x) + \delta f'(x) + \delta^{2} f''(x)/2 + O(\delta^{3}) \\
f(x - \delta) & \approx f(x) - \delta f'(x) + \delta^{2} f''(x)/2 + O(\delta^{3})
\end{align}
$

and solve for the desired expression. For instance, if we substract the second equation from the fist, we get the central difference formula for the first derivative:

$
\begin{align}
f'(x) = \frac{f(x+\delta) - f(x-\delta)}{2\delta} + O(\delta^{2})
\end{align}
$

Similarly, if we add both equations, we get the central difference formula for the second derivative:

$
\begin{align}
f''(x) = \frac{f(x + \delta) -2f(x) + f(x - \delta)}{2\delta^{2}} + O(\delta^{2})
\end{align}
$


This is all very fine, but in practice, we have situations where we know $f(x)$, $f(x + \delta)$ and $f(x - \delta/2)$ and we want the derivative $f'(x)`$.   
Using substitutions can get hairy. An elegant way out if this is to use Newton's divided difference interpolation formula: 

$
\begin{align}
f(x) \approx f_{0} + (x - x_{0}) f[x_{0}, x_{1}] + (x - x_{0})(x - x_{1}) f[x_{0}, x_{1}, x_{2}] + \ldots
\end{align}
$

where $f[x_{0}, x_{1}] = \frac{f_{1} - f_{0}}{x_{1} - x_{0}}$ and 

$
\begin{align}
f[x_{0}, \ldots, x_{k}] = \frac{f[x_{1}, \ldots, x_{k}] - f[x_{0}, \ldots, x_{k-1}]}{x_{k} - x_{0}}.
\end{align}
$

By differentiating $f(x)$, we can find the derivative at any place, knowing the value of the function at $x_{0}$, $x_{1}$ and $x_{2}$:

$
\begin{align}
f'(x) \approx f[x_{0}, x_{1}] + (2x - x_{1} - x_{0}) f[x_{0}, x_{1}, x_{2}] + \ldots
\end{align}
$

So, to being, let's assume $x_{i+1} = x_{i+h}$, we can simplify things quite a bit:

$
\begin{align}
f[x_{0}, x_{1}] & = \frac{f_{1} - f_{0}}{h} \\
f[x_{0}, x_{1}, x_{2}] & = \frac{f_{0} - 2f_{1} + f_{2}}{2h^{2}}
\end{align}
$

So that taking $x=x_{0}$ yields:

$
\begin{align}
f'_{0} & = \frac{f_{1} - f_{0}}{h} - h \frac{f_{0} - 2f_{1} + f_{2}}{2h^{2}} \\
       & = \frac{-3f_{0} + 4f_{1} -f_{2}}{2h} 
\end{align}
$

And similarly for $x=x_{1}$ and $x=x_{2}$:

$
\begin{align}
f'_{1} & = \frac{-f_{0} + f_{2}}{2h} \\
f'_{2} & = \frac{f_{0} - 4f_{1} + 3f_{2}}{2h}
\end{align}
$

The three last formulas correspong to the forward, central and backward difference formulas of the second order.  

Now this is very useful to compute the derivative at points that are not multiples of $h$. For instance, imagine we have ice velocities at points $i$ and $i+1$, but we know that at point $i-1/2$, the velocity is $0$ because of a land mass. The distance between points $i$ and $i+1$ is $h$, and between $i-1/2$ and $i$ the distance is $h/2$, so:

$
\begin{align}
f[x_{i-1/2}, x_{i}] & = \frac{2f_{i}}{h} \\
f[x_{i-1/2}, x_{i}, x_{i+1}] & = \frac{2}{3h^{2}}(f_{i+1} - 3f_{i})
\end{align}
$

and 

$
\begin{align}
f'_{i}  & = \frac{2f_{i}}{h} + \frac{h}{2}\frac{2}{3h^{2}}(f_{i+1} - 3f_{i}) \\
        & = \frac{6f_{i} + f_{i+1} - 3f_{i}}{3h} \\
        & = \frac{3f_{i} + f_{i+1}}{3h}
\end{align}
$

## Practical Cases
------------------

Assume u and v are on a staggered C grid. Let's say we want to compute $\partial_{x} u$ at the grid center. This can be computed simply using central differences: $(u_{i, j} - u_{i+1, j}) / \delta$. Cast in the divided difference notation, we would have: $f_{0} = u_{i, j}$, $f_{1} = $, $f_{2} = u{i+1,j}$ and $x_{0} = 0$, $x_{1} = \delta / 2$, $x_{2} = \delta$.
 
Now assume we want $\partial_{x} v$ at the node center. First, we compute the vertical average to obtain values at the nodes:  
$v'_{i, j} = (v_{i, j} + v_{i, j+1}) / 2$. Then we apply the divided difference formula with: $f_{0} = v'_{i-1, j}$, $f_{1} = v'_{i, j}$, $f_{2} = v'_{i, j+1}$ and $x_{0} = 0$, $x_{1} = \delta$, $x_{2} = \delta$. 

Now imagine there is a land mass on the node $i-1, j$, we would then have:  
$f_{0} = 0$, $f_{1} = v'_{i, j}$, $f_{2} = v'_{i+1, j}$ and $x_{0} = 0$, $x_{1} = \delta / 2$, $x_{2} = 3\delta / 2$. 

And similarly if the mask is on the right side. Now the idea is to write functions that deal gracefully with all these cases wihout too much special casing, which is numerical recipe for bugs. 

# Advection Schemes
-----------------


This section describes the upwind (upstream) advection scheme that is currently used in the model. Next, we describe other advection schemes, their advantages and disadvantages.  

Some of the material in this section is inspired from the [MITGCM](https://mitgcm.readthedocs.io/en/latest/).


## Definitions
--------------
-   **Courant number**  
The Courant number is defined as $C = u \Delta t / \Delta x$, the velocity times the time step divided by the grid resolution. It is used to characterize advection conditions. If $C = 1$, it means that a quantity is advected over an entire grid cell in a single time step, which can corrupt some advection schemes. It is hence important to check the typical values of the Courant number to make sure the advection scheme is appropriate.

- **Advection schemes**   
are computational methods used to solve the differential equation  
    $ \partial_{t} A + u \partial_{x} A = 0 $

Look at [Numerical Models of Oceans and Oceanic Processes](https://www.sciencedirect.com/bookseries/international-geophysics/vol/66/suppl/C) for a good reference. 

------------

## Upwind Advection Scheme
--------------------------

This is the simplest advection scheme. It solves the advection equation using backward space difference and forward time differences. This is what is currently used in the model to advect ice thickness and ice covered area. The main problem with this scheme is that it is diffusive. That is, if we write the Taylor series of the differential equation, we see that the solving the finite difference approximation amounts to solving an advection-diffusion equation. In other words, the scheme adds a diffusive term that smooths the advected quantity. 

# Free Drift Solution to the Sea Ice Momentum Equation
------------------------------------------------------
_David Huard, Bruno Tremblay  
June 2008_  



## Context
-------

The free drift solution to the momentum equation yields the velocity of the ice when neglecting the rheology of the ice, that is, the internal forces occurring in the ice. This solution is used in the model as a first approximation to the full solution, a step that is necessary to ensure that the non linear solver has an adequate guess to being its first iteration.

## Equations
---------

The free drift equation is written as:

$
\begin{align}
0 = -\rho_{i} h f \vec{k} \times \vec{u}_{i} + \tau_{a}  - \tau_{w} - \rho_{i} h g \nabla H_{d}
\end{align}
$

where $\rho_{i}$ is the ice density, $h$ the ice thickness, $\vec{k}$ the vector normal to the surface (pointing up), $\vec{u}_{i}$ the ice velocity vector, $\tau_{a}$ the wind shear stress, $\tau_{w}$ the water shear stress and $H_{d}$ the sea surface height. Note that the sign of the water stress is the inverse of that of the wind (air) stress.  

Let's break this down!

### Coriolis Force 

The first term on the left is the Coriolis force and its components in the $x$ and $y$ directions are given by:

$
\begin{align}
-\rho_{i} h f \vec{k} \times \vec{u}_{i} = 
   \left[\begin{matrix}
   \rho_{i} h f u_{y} \\
   -\rho_{i} h f u_{x} \\
   0
   \end{matrix} \right],
\end{align}
$

where $f$ is the Coriolis parameter and is given by $f = 2 \omega \sin \phi$, $\omega$ being Earth's angular velocity and $\phi$ the latitude. Since the model is currently restricted to the Arctic, the Coriolis parameter is taken constant. From now on, we will ignore the velocity in third dimension who is always 0 since we are considering motion on a plane.

### Surface stress

The second and third terms are surface stresses, and are defined in the model as quadratic laws. That is, the stress depends on the square of the velocity.  Both the air and water stresses can be written in the general form:

$
\begin{align}
\tau = \rho C_{d} |\vec{du}| (\vec{du} \cos \theta + \vec{k} \times \vec{du} \sin \theta),
\end{align}
$

where $\rho$ is the dragging medium density, $C_{d}$ the drag coefficient, $\vec{du}$ the relative speed between the forcing medium and the ice and `\theta` the turning angle. A simpler for surface stress can be obtained if we assume that the stress is related linearly to the forcing medium velocity : 

$
\begin{align}
\tau = \rho C_{d} (\vec{du} \cos \theta + \vec{k} \times \vec{du} \sin \theta), 
\end{align}
$

Since wind speeds are much higher than the ice speeds, we assume for the air stress that $\vec{u}_{a} - \vec{u}_{i} \approx \vec{u}_{a}$.

Specifically, for the quadratic law, the components of the air stress are:

$
\begin{align}
\tau_a =\left[ \begin{matrix}
   \rho_{a} C_{da} \sqrt{v_{a}^{2} + u_{a}^{2}} \left (\cos \theta_{a} u_{a} - \sin \theta_{a} v_{a} \right) \\
   \rho_{a} C_{da} \sqrt{v_{a}^{2} + u_{a}^{2}} \left (\cos \theta_{a} v_{a} + \sin \theta_{a} u_{a} \right)
   \end{matrix}  \right]
\end{align}
$

and the components of the water stress are:

$ 
\begin{align}
-\tau_w = \left[ \begin{matrix}
   -C_{dw} \rho_{w} \sqrt{\left(v_{i} - v_{w} \right)^{2} + \left(u_{i} - u_{w}\right)^{2}} \left(\left(u_{i} - u_{w} \right) \cos \theta_{w} - \left(v_{i} - v_{w} \right) \sin \theta_{w} \right) \\
   -C_{dw} \rho_{w} \sqrt{\left(v_{i} - v_{w} \right)^{2} + \left(u_{i} - u_{w}\right)^{2}} \left(\left(v_{i} - v_{w} \right) \cos \theta_{w} + \left(u_{i} - u_{w} \right) \sin \theta_{w} \right)
   \end{matrix} \right]
\end{align}
$


### Sea Surface Tilt 

The last term, the sea surface tilt, is expressed using the ocean geostrophic
currents. Geostrophic water currents are the currents that occur when we
neglect all forces except the Coriolis force and the gravity:

$ 
\begin{align}
f \vec{u}_{w} = g \vec{k} \times \nabla H.
\end{align}
$

We use this relation to express the force due to sea surface tilt by the geostrophic currents, that is,

$ 
\begin{align}
\vec{k} \times f \vec{u}_{w} & = \vec{k} \times (g \vec{k} \times \nabla H) \\
f \vec{k} \times \vec{u}_{w} & = -g \nabla H
\end{align}
$

so that the last term of the first equation can be written as

$ 
\begin{align}
- \rho_{i} h g \nabla H_{d} = \rho_{i} h f \vec{k} \times \vec{u}_{w} 
   = \left [ \begin{matrix}
   -\rho_{i} h f v_{w} \\
   \rho_{i} h f u_{w}
   \end{matrix} \right].
\end{align}
$

Bringing everything together, we have:

$ 
\begin{align}
0 =\left[ \begin{matrix}
   C_{dw} \rho_{w} \sqrt{v^{2} + u^{2}} \left(v \sin \theta_{w} - u \cos \theta_{w} \right) + \rho_{i} f h v + C_{da} \rho_{a} \sqrt{v_{a}^{2} + u_{a}^{2}} \left(u_{a} \cos  \theta_{a} - v_{a} \sin \theta_{a} \right)\\
   -C_{dw} \rho_{w} \sqrt{v^{2} + u^{2}} \left(u \sin \theta_{w} + v \cos \theta_{w} \right) - \rho_{i} f h u + C_{da} \rho_{a} \sqrt{v_{a}^{2} + u_{a}^{2}} \left(u_{a} \sin \theta_{a} + v_{a} \cos \theta_{a} \right)
   \end{matrix} \right]
\end{align}
$

where $u = u_{i} - u_{w}$ and $v = v_{i} - v_{w}$. If we call the right hand side of the momentum equation $F(\vec{u})$, our next task is to find $\vec{u}$ such that $F(\vec{u}) = 0$, that is, find the roots of the equations.


## Solving for $\vec{u}$
--------------------

To find a solution for $\vec{u}$ the ice velocity relative to the ocean currents, we use the Newton-Raphson (N-R) method. This method is a generic method used to find the zero of a function, that is, the value $x$ for which $f(x) = 0$. To do so, N-R uses a starting guess $x_{n}$, and computes the next guess by iteration using

$ 
\begin{align}
x_{n+1} = x_{n} - \frac{f(x_{n})}{f \prime(x_{n})}.
\end{align}
$

For functions of more than one variable, the idea is similar, but the derivative is replaced by the Jacobian matrix :

$ 
\begin{align}
\vec{u}_{n+1} = \vec{u}_{n} - J^{-1}(\vec{u}_{n}) F(\vec{u}_{n}).
\end{align}
$

For our particular case, we have a vector function $\vec F(u_{x}, u_{y})$ that has two components $F_{x}$ and $F_{y}$. The Jacobian can we written as follows:

$ 
\begin{align}
J =
   \left[ \begin{matrix}
   J_{11}  & J_{12} \\
   J_{21} & J_{22} \end{matrix} \right]
   =
   \left[ \begin{matrix}
   \frac{\partial F_{x}}{\partial u} &  \frac{\partial F_{x}}{\partial v} \\
   \frac{\partial F_{y}}{\partial u} &  \frac{\partial F_{y}}{\partial v}
   \end{matrix}\right].
\end{align}
$

For a quadratic surface stress, the terms of the Jacobian are:

$ 
\begin{align}
\partial_{u} F_{x} & = \rho_{w} C_{dw} \frac{\left[u v \sin \theta_{w} - \left(v^{2} + 2 u^{2} \right) \cos \theta_{w} \right]}{\sqrt{v^{2} + u^{2}}} \\
\partial_{v} F_{x} & = \rho_{w} C_{dw} \frac{\left(2v^{2} + u^{2} \right) \sin \theta_{w} - u v \cos \theta_{w}}{\sqrt{v^{2} + u^{2}}} + \rho_{i} h f \\
\partial_{u} F_{y} & = -\rho_{w} C_{dw} \frac{\left(v^{2} + 2u^{2} \right) \sin \theta_{w} + u v \cos \theta_{w}}{\sqrt{v^{2} + u^{2}}} - \rho_{i} h f\\
\partial_{v} F_{y} & = -\rho_{w} C_{dw} \frac{\left[u v \sin \theta_{w} + \left(2v^{2} + u^{2} \right) \cos \theta_{w} \right]}{\sqrt{v^{2} + u^{2}}}
\end{align}
$
 
while they are much simpler if a linear surface stress is assumed:

$ 
\begin{align}
\partial_{u} F_{x} & = -\rho_{w} C_{dw} \cos \theta_{w} \\
\partial_{v} F_{x} & = \rho_{w} C_{dw} \sin \theta_{w} + \rho_{i} h f\\
\partial_{u} F_{y} & = -\rho_{w} C_{dw} \sin \theta_{w} - \rho_{i} h f\\
\partial_{v} F_{y} & = -\rho_{w} C_{dw} \cos \theta_{w}
\end{align}
$

Knowing that for a $2 \times 2$ matrix $J$, the inverse is computed using:

$ 
\begin{align}
J^{-1} =
    \frac{1}{\det{J}}\left[ \begin{matrix}
    J_{22} & -J_{12} \\
    -J_{21} & J_{11}
    \end{matrix} \right].
\end{align}
$

where $\det{J} = J_{11} J_{22} - J_{21} J_{12}$. Assuming $M = J^{-1}$ we have everything we need to know to compute $\vec{u}$:

$ 
\begin{align}
\det{J} & = \partial_{u} F_{x} \partial_{v} F_{y} - \partial_{u} F_{y} \partial_{v} F_{x} \\
M_{11} & = \frac{\partial_{v} F_{y}}{\det{J}} \\
M_{12} & = \frac{-\partial_{y} F_{x}}{\det{J}} \\
M_{21} & = \frac{-\partial_{u} F_{y}}{\det{J}} \\
M_{22} & = \frac{\partial_{u} F_{x}}{\det{J}}
\end{align}
$

so that

$ 
\begin{align}
u_{n+1} & = u_{n} - (M_{11} F_{x} (\vec{u_{n}}) + M_{12} F_{y} (\vec{u_{n}}))  \\
v_{n+1} & = v_{n} - (M_{21} F_{x} (\vec{u_{n}}) + M_{22} F_{y} (\vec{u_{n}}))
\end{align}
$

Remember that this is the solution for the relative ice velocity,not the ice velocities themselves. 

## Original Numerical Implementation
---------------------------------

In the original code `UV_solve_NR_B`, here are the steps performed to find a solution to the free drift equation:

  1. Set $\vec{u} = 0$

  2. Compute the relative linear free drift solution on the B-grid. The linear free drift solution refers to the solution to the momentum equation when the rheology term is absent and where the surface stresses are proportional to the medium's speed. That is, instead of the quadratic law for the stress,  we have $\tau = \rho C (\vec{u} \cos \theta + \vec{k} \times \vec{u} \sin \theta$. Since the function is linear, the differentation is straightforward and the N-R step leads us directly to the analytical solution for $\vec{u}$. This solution will be the first guess for the nonlinear solution. 

  3. Loop over the solution to the nonlinear free drift equation until a convergence criteria is met or a maximum number of loops is attained. This solution is described above, and consists in computing, for each $u$ and $v$ in the matrix, the coefficients of the inverse Jacobian matrix, the solution to $\vec{F}(\vec{u})$ at the last iteration and corresponding $\delta \vec{u}$. 

  4. Interpolate the solution on the C-grid. 

  5. Apply the boundary conditions: the velocity and its derivative along the two directions are 0. 


## Proposed Numerical Implementation
---------------------------------

A modular approach is to code a simple primitive procedure for each stress component: Coriolis, linear surface stress, quadratic surface stress, sea surface tilt. A free_drift_momentum_eq subroutine can then be written using calls to the primitive functions. Since there are no interactions between adjacent grid cells, the primitive procedures can be declared **elemental**. The main advantage is that those procedures can be reused in other parts of the code without modification. A `linear_freedrift_solve` function would solve the momentum equation for a linear drag law. The result would then be used as a starting guess in an iterative `quadratic_freedrift_solve`. The latter subroutine would itself call a `quadratic_freedrift_step` which would compute a $\Delta \vec{u}$ for one step using the momentum
equation and the computed Jacobian at the current $\vec{u}$. 

- **elemental**
Elemental statements are applied to procedure operating on scalars to tell the compiler that it has no side effect and can be safely applied in a loop to all the elements of an array. In other words, it vectorizes a scalar procedure. 

## Notes
-----

In this text, $M_{11}$, $M_{12}$, $M_{21}$ and $M_{22}$ refer to the opposite of variables `A11`, `A12`, `A21` and `A22` in the code. I chose this different notation
since there is mix up potential: in the paper, A is used both for the ice covered area and a $2 \times 2$ matrix multipling the relative ice velocity; in the code, Aij is the negative of the inverse of the Jacobian.


# Rheology
----------

## Definitions
--------------
    
- **Constitutive law**  
A constitutive law is an equation relating the strain of some material to the imposed stress.

- **Deviatoric stress tensor $\sigma \prime_{ij}$**  
The difference between the stress tensor and the isotropic pressure $p$. Note that this mechanical pressure cannot always be identified with the thermodynamic pressure. 

- **Isotropic pressure $p$**  
The average normal pressure.
$$
p = -\frac{1}{2}(\sigma_{kk}).
$$

- **Principal stresses**   
The principal stresses are the eigenvalues $\sigma_{1}$, $\sigma_{2}$ of the stress tensor with the convention that $\sigma_{1} > \sigma_{2}$. This can be understood as the diagonal elements of the stress tensor when the coordinate system is rotated such that the off-diagonal elements are 0:

$$
\sigma^{\prime_{ij}} = \left[ \begin{matrix}
   \sigma_{1} & 0 \\
   0 & \sigma_{2}
   \end{matrix}\right]
$$
<br><br><br>

![principal_stress](source/principal_stress.png "Principal stress")

          Principal stresses.


- **Rheology**  
Rheology is the study of flow and matter under the influence of a stress. 


- **Shear strain rate**  
The rate at which two sides of an elemental square close together:

$$
\dot{\gamma}_{ij} = \left(\frac{\partial u_{i}}{\partial x_{j}} + \frac{\partial u_{j}}{\partial x_{i}} \right)
$$


- **Shear stress $\tau_{ij}$**  
The pure shearing stress. 


- **Stress**  
A stress is a force per unit area on a surface. The stress field is the distribution of internal stresses that balance a set of external forces. 

- **Stress tensor $\sigma_{ij}$**    
A matrix describing the stresses on an infinitesimal surface (or cube in  3D). The first index $i$ refers to the direction of the normal to the place  the stress acts on. The second index $j$ refers to the direction of the stress. For example, $\sigma_{xx}$ refers to the stress on the $y$-plane in the $x$ direction, while $\sigma_{xy}$ refers to the stress on the $y$-plane in  the $y$ direction. A normal stress refers to a stress perpendicular to the  surface ($\sigma_{xx}, \sigma_{yy}$), while a shear stress refers to a stress parallel to the surface ($\sigma_{xy}$). By convention, tensile  stresses are positive while conpressive ones are negative. In a 2-D coordinate system, stresses are described using a $2 \times 2$ matrix $\sigma_{ij}$ 
$$
\begin{align}
\sigma_{ij} = \left[ \begin{matrix}
   \sigma_{11} & \sigma_{12} \\
   sigma_{21} & \sigma_{22}
   \end{matrix}\right]
\end{align}
$$<br>
To satisfy the equilibium of moments (no rotation), the stress tensor needs to be symmetric, ie, $\sigma_{ij} = \sigma_{ji}$. 
<br><br>

![stress](source/stress.png "Stress")

          State of stress. 

- **Stress invariants**   
The stress invariantes are defined as 
$$
\begin{align}
\dot{\epsilon}_{I} & \equiv \dot{\epsilon}_{1} + \dot{\epsilon}_{2} = \text{divergence} \\
\dot{\epsilon}_{II} & \equiv \dot{\epsilon}_{1} - \dot{\epsilon}_{2} = \text{maximum shear rate}
\end{align}
$$<br>
If we assume that $\dot{\epsilon}_{I}$ and $\dot{\epsilon}_{II}$ are the real and complex components of a complex variable, and $\theta$ the angle between the complex vector and the real axis, then $\theta = 0$, $\pi / 4$, $\pi / 2$, $3\pi / 4$ and $\pi$ correspond to pure divergence, uniaxial extension, pureshear, uniaxial contraction and pure convergence respectively. If you find the last statement puzzling, please stop and think about it before going on.   

- **Strain $\epsilon$**  
The deformation of a material $\frac{dl}{l}$. In 3D you can imagine the strain as the displacement of the corner of a cube under stress. 
<br>

![types_of_motion](source/types_of_motion.png "Types of motion")

<br>

- **Strain rate $\dot{\epsilon}$**  
The rate of deformation over time $\frac{1}{l}\frac{dl}{dt} = \frac{v}{l}$. In 2-D, this is described by the strain rate tensor:
$$
\dot{\epsilon}_{ij} = \frac{1}{2}\left(\frac{\partial u_{i}}{\partial x_{j}} + \frac{\partial u_{j}}{\partial x_{i}}\right)
$$
The diagonal elements of the strain rate tensor are the normal strain rates. The off-diagonal elements are one half the shear strain rate components ($\dot{\gamma}_{ij}$). 
Note that some authors denote the strain rate simply by $\epsilon$. 

- **Velocity gradient tensor**  
The velocity gradient tensor is written as  $\partial u_{i} / \partial x_{j}$. The symmetrical part of the tensor is the strain rate tensor, and the antisymmetrical part is the rotation tensor, which is related to the vorticity ($\omega$).   
$$
\frac{\partial u_{i}}{\partial x_{j}}   = \dot{\epsilon}_{ij} + \Omega_{ij}
$$
where $\Omega = \frac{1}{2} \nabla \times \vec{u} = \frac{1}{2} \omega$. 


- **Yield stress**  
The minimum stress that must be applied to initiate significant flow and a significant drop in viscosity. 

- **Newtonian fluid**  
A fluid whose viscosity is independent of the shear conditions. The resistance between two layers of fluid is proportional to the difference in speed of these fluids. The constitutive law for a compressible Newtonian fluid is 
$$
\sigma_{ij} = -2 \eta \dot{\epsilon}_{ij} + \left(\frac{2\eta}{3} - \zeta \right) \dot{\epsilon}_{kk} \delta_{ij}
$$

- **Elastic material**  
A material is elastic if the deformation follows the applied stress. A perfect elastic material has a stress proportional to the strain. 

- **Viscous material**  
A material is viscous if the deformation increases linearly at a constant stress. In other words, you need to keep pushing on it to make it move. 

- **Yield curve**  
Curve describing the stress at which a material starts to yield.    

- **Bulk viscosity $\zeta$**  
Viscosity associated with volume expansion, also called volume viscosity or Lamé's constant.

- **Shear viscosity $\eta$**  
Viscosity when the applied stress is a shear stress. 


## Elliptical Yield Curve
----------------------

$
\begin{align}
\zeta & = P / 2 \Delta \\
\eta & = \zeta / e^{2} \\
\Delta & = [(\dot{\epsilon}_{11}^{2} + \dot{\epsilon}_{22}^{2})(1 + 1 / e^{2}) + 4e^{-2}\dot{\epsilon}_{12}^{2} + 2 \dot{\epsilon}_{11} \dot{\epsilon}_{22}(1 - 1 / e^{2})]^{1/2}
\end{align}
$

## Viscous Plastic Constitutive Law
--------------------------------

$
\begin{align}
2 \eta \dot{\epsilon}_{ij} + [\zeta - \eta] \dot{\epsilon}_{kk} \delta_{ij} - P\delta_{ij} / 2
\end{align}
$

## Internal Stress Force
---------------------
The force components due to internal stresses in the ice are calculated from $F_{i} = \partial \sigma_{ij} / \partial x_{j}$ (summation is implied on the repeated indices.)

# Solving the Momentum Equation
------------------------------
*David Huard
August 5, 2008*

You should read the section on solving the free drift equation before tackling this section on solving the complete momentum equation (ME).  


The ME is a non-linear equation. Furthermore, we do not have an analytical formula for the Jacobian, hence, we cannot use the same technique as in the free drift equation (Newton's algorithm). The current strategy to solve for $u$ is to iterate over the following steps:

 1. Linearize the equation. This means that all non-linear terms in $u$ are set to the value of $u$ at the last iteration. 
 2. Once the non-linear terms are fixed, we have a linear equation in $u$. We can use a linear solver to find $u$ pretty efficiently. This is now done using GMRES. 
 3. Use the solution to the linear equation to update $u$, and return to step 1. The way this is currently done is rather simple, but we are working on implementing a Jacobian Free Newton Method to improve this. 

The number of iteration (linearize, solve the linear equation, update) is called an outer loop in our jargon. To get a fully converged solution to the non-linear equation can take hundreds of outer loops. 


## Constant Terms
-----------------

The ME contains terms that are independent of the ice velocity $u$. These are:
 
 * The wind stress, which only depends on the wind velocity. 
 * The sea ice tilt, which only depends on the water velocity. 
 * The pressure term, $\nabla \cdot P$. 

## Linear Terms
---------------

The linear terms in the momentum equations are 

 * The Coriolis force. 
 * The linear component of the water stress.

## Non-Linear Terms
------------------

The non-linear terms in the ME are
 
 * The non-linear component of the water stress, $|\vec{u} - \vec{u_w}|$. 
 * The strain rate terms in the constitutive law.  
$
\nabla \cdot (2 \eta \dot{\epsilon}_{ij} + [\zeta - \eta] \dot {\epsilon}_{kk} \delta_{ij})
$

## Relation to the Code
----------------------

Although I am not 100% sure to understand everything in the code, here is rough breakdown of what happens before and during an outer loop: 

**Before**: 
 * Load wind velocity, compute wind stress and the sea ice tilt. 
 * Compute the stress related to the pressure. 
 * Put those terms (air stress, sea ice tilt, pressure) in the variables `R1p` and `R2p`.   

**During**: 
 * Compute the non-linear terms in the water stress (`r1ppr2pp`).
 * Compute the viscous coefficients (`viscouscoefficients`). 
 * Solve the linear equation (`prep_fgmres`). 
 * Update the ice velocity (`r1ppr2pp`).


To facilitate the implementation of other solvers and clean up the code, I propose to separate those chunks into well identified individual procedural subroutines:
 
 * Subroutine for the wind stress (done).
 * Subroutine for the sea ice tilt (done). Since this does not vary at all, ever, store that in a variable. 
 * Subroutine for the ice pressure (done, untested). 
 * Subroutine to gather the independent terms (todo).
 * Subroutine to compute the non linear term in the water stress. 
 * Subroutine to compute the linear term in the water stress. 
 * Subroutine to compute the viscous coefficients (`viscouscoefficient`, convert to procedural).
 * Subroutine to compute the left hand side (`matvec`, convert to procedural)
 * Subroutine to gather the right hand side (`r1ppr2pp`, clean up). 
 * Subroutine to compute the linear and non-linear residuals (convert to procedural). 