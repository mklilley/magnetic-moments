---
jupyter:
  jupytext:
    formats: ipynb,src//md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.2'
      jupytext_version: 1.3.3
  kernelspec:
    display_name: Julia 1.3.1
    language: julia
    name: julia-1.3
---

# Stern and Gerlach simulation


It was discovered experimentally,[in 1922 by Stern and Gerlach](https://www.feynmanlectures.caltech.edu/II_35.html#Ch35-S2) (SG), that the spin of an electron along any direction is quantised, taking values of either $\hbar/2$ or $-\hbar/2$. The experiment involved firing silver atoms into an inhomogeneous magnetic field and observing a splitting of the beam in the direction of the field inhomogeneity.

Since then, people have described an idealised system of many SG systems in parallel that are capable of splitting apart and recombining the particle beams together (often called a Stern-Gerlach interferometer) to create interference effects that are commonly discussed in the two slit experiment (see e.g. the [Feynman lectures](https://www.feynmanlectures.caltech.edu/III_05.html#Ch5-S4)). 

To my surprise, it is only [very recently](https://arxiv.org/abs/1801.02708) that Stern-Gerlach interferometer experiments are being considered practical enough to implement and yet they are used frequently as a pedagogical example. It is therefore worth first developing an intuition for what one might expect from such SG interferometer experiments using only classical physics. That is the purpose of this notebook.


## Libraries

```julia
# need this to use the cross product later
using LinearAlgebra 

# using Plots http://docs.juliaplots.org/latest/install/
# https://gist.github.com/gizmaa/7214002 - examples uing PyPlot in Julia
using PyPlot
```

## Physical constants

```julia
# The g-factor https://en.wikipedia.org/wiki/G-factor_(physics)
g=2

# Bohn magneton
mu_B = 9.27400968e-24;

# electron charge
e = 1.602176634e-19

# electron mass (not needed but I have it here in case I need it later)
m_e = 9.1093837015e-31

# magnetic moment of the electron (ignoring g-factor corrections for now, i.e. g=2)
mu_e =  mu_B

# mass of silver atom
m_Ag = 107.9*1.66053906660e-27

# gyromagnetic ratio https://en.wikipedia.org/wiki/Gyromagnetic_ratio
gyro =  g*e/(2.0*m_Ag)

# magnetic moment of silver (same as a singleelectron because Ag has 1 electron in outer shell)
mu_Ag =  mu_e;
```

## Experimental set-up
A 2009 paper by França on "[The Stern–Gerlach Phenomenon According to Classical Electrodynamics](https://link.springer.com/article/10.1007/s10701-009-9338-1)" describes a lot of the set-up and parameters, e.g.

![Stern Gerlach setup](SG-setup.png)


An approximation for the magnetic field inside the electromagnet (length 3.5cm) is described by França's Eq 23:

$$
B = (-\beta x, 0, B_0 + \beta z)
$$

with $ \beta = 1800$ explicitly stated at the bottom of page 1186 and $B_0$ inferred from França's Fig 3 and checked in Fig. 5 of [Rabi's original paper](https://link.springer.com/article/10.1007/BF01339837).

```julia
# max B field inside the electromagnet
B0 = 1.4

# x and z B field gradient inside the electromagnet
beta = 1800

# electromagnet dimesion in y
# note, I cannot find dimensions for x and z
machine_dim_y = 0.035;
```

To start with, we'll simulate a particle moving through the electromagnet only and deal with what happens at the magnet entrance later on.

Let's explore the B field given França's parameters.

```julia
function B(r)
    B = zeros(3)
    B[1] = -beta*r[1]
    B[3] = B0 + beta*r[3]
    return B
end;
```

We can visualise at B field lines and also the strength of B in the x,z plane:

```julia
r = [[j,0,i] for i=-0.01:0.0001:0.01, j=-0.01:0.0001:0.01];
x = [x for (x,y,z) in r]
z = [z for (x,y,z) in r];
```

```julia
B_calc = B.(r)
B_calc_x = [Bx for (Bx,By,Bz) in B_calc] 
B_calc_z = [Bz for (Bx,By,Bz) in B_calc]
B_norm = norm.(B.(r));
```

```julia
st = streamplot(x*1000, z*1000 ,B_calc_x, B_calc_z, linewidth=3*Bnorm/maximum(Bnorm), color=B_calc_z)
cb=colorbar(st.lines)
cb.set_label(L"B_z")
xlabel("x (mm)")
ylabel("z (mm)");
```

```julia
con = contour(x.*1000, z.*1000, B_norm, levels=20)
xlabel("x (mm)")
ylabel("z (mm)")
clabel(con, inline=1, fontsize=10);
```

This field is very strong! To get a sense of what equipment generates these fields, below are figures taken  from the book `An Introduction to Quantum Physics` by French. This book is referenced in the França paper as a "modern" SG experiment.

![Modern Stern Gerlach setup 1](SG-french-01.png)

![Modern Stern Gerlach setup 2](SG-french-02.png)


## Equations of motion for a magnetic moment in a B field

If we consider an infinitesimal current loop, then the magnetic moment $\mu$ is:

$$
\mu = IA
$$

where $I$ is the current of the loop and $A$ is its area. The force on the loop is:

$$
\nabla(\mu\cdot B)
$$

and the torque is:

$$
\tau = \mu \times B
$$

cf. [wikipedia](https://en.wikipedia.org/wiki/Magnetic_moment#Effects_of_an_external_magnetic_field).

Recalling that:
1. a current loop has [angular momentum](https://en.wikipedia.org/wiki/Magnetic_moment#Relation_to_angular_momentum) $L=\mu/\gamma$ where $\gamma$ is the gyromagnetic ratio
2. the magnetic moment is not a function of position

we can write some equations of motion:


$$
\frac{dv}{dt} = \frac{1}{m} \mu \nabla B
$$


$$
\frac{d\mu}{dt} = \gamma \mu\times B
$$


where $\nabla B$ is the [gradient of the vector]( https://en.wikipedia.org/wiki/Gradient#Gradient_of_a_vector) B (i.e. a tensor) with components $\nabla B_{ij} = \partial B_i / \partial x _i$ which we can calculate using the function below.

```julia
function gradB(r)
    gradB = zeros(3,3)
    gradB[3,3] = beta
    gradB[1,1] = -beta
    return gradB
end;
```

From the equations of motion above we can expect:
1. the magnetic moment to move in the direction of gradients in the field
2. the orientation of the magnetic moment will precess about the direction of the magnetic field (see [Larmor precession](https://en.wikipedia.org/wiki/Larmor_precession))


## Set-up initial values for the simulation

We'll start the particles off with velocity only in the y direction at $v_0 = 600 \ m/s$ (see end of page 1181 of França article).

We'll initialise a random orientation of the magnetic moment using a spherical coordinate system:

$$ \mu = [\sin(\theta)\cos(\phi), \ \sin(\theta)\sin(\phi), \ \cos(\theta)]\mu_{Ag} $$

where $\theta$ and $\phi$ are the [spherical angle coordiantes](https://en.wikipedia.org/wiki/Spherical_coordinate_system) corresponding to the ISO convention. In order to ensure a uniform distribution we must choose

$$
\theta = \arccos(2\times rand()-1)
$$


$$
\phi = 2\pi\times rand()
$$

as dicussed [here](https://math.stackexchange.com/a/1586583).

Because we expect some precessional motion, we should choose a time-step that resolves it. The [Larmor frequency](https://en.wikipedia.org/wiki/Larmor_precession#Larmor_frequency) $\omega = \gamma B$ is a good guide for this. We therefore need a timestep $\Delta t < 2\pi/\omega$:



```julia
2*pi/(gyro * B0)
```

```julia
v0 = 600.0
v = [0,v0,0]
r = [0,0,0]

theta = acos(2*rand()-1)
phi = rand()*2*pi
mu = [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)]*mu_Ag

dt = 1.0e-7
t_max = machine_dim_y/v0
times = collect(0:dt:t_max)


r_save = zeros(length(times),3)
mu_save = zeros(length(times),3);
```

## Solving the equations of motion


To numerically solve equations of motion:


$$
\frac{dv}{dt} = \frac{1}{m} \mu \nabla B
$$


$$
\frac{d\mu}{dt} = \gamma \mu\times B
$$


We employ an algorithm that's a hybrid of [leapfrog](https://en.wikipedia.org/wiki/Leapfrog_integration) and the [Boris method](https://www.particleincell.com/2011/vxb-rotation/).

```julia
function frog_boris(r,v,mu)
    
    # leapfrog kick
    v += (0.5*dt*(mu'*gradB(r))/m_Ag)'
    
    # Boris
    Br = B(r)
    s = 2.0/(1+(norm(Br)*gyro*dt*0.5)^2)
    mu_prime = mu + 0.5*dt*gyro*cross(mu,Br)
    mu += 0.5*dt*gyro*s*cross(mu_prime,Br)
    
    # leapfrog drift
    r += dt*v
    
    # leapfrog kick
    v += (0.5*dt*(mu'*gradB(r))/m_Ag)'
    
    return r, v, mu
end;
```

```julia
for (i,t) in enumerate(times)
    
    r_save[i,:] = r
    mu_save[i,:] = mu
    
    r, v, mu = frog_boris(r, v, mu)
    
end
```

## Visualise the single particle results

```julia
# In case we simulate lots of points and we don't want to plot them all because of memory,
# we can plot every 'n' points
n = 1
times_plot = times[1:n:end]./1e-6
r_save_plot = r_save[1:n:end,:]./1e-3
mu_save_plot = mu_save[1:n:end,:]./mu_B;
mu_norm_plot = [norm(mu_save_plot[i,:]) for i=1:size(mu_save_plot,1)]

figure(figsize=(8,10))

subplot(411)
plot(times_plot,(r_save_plot[:,3]), label="z");
ylabel("z (mm)")

subplot(412)
plot(times_plot,(r_save_plot[:,2]), label="y");
ylabel("y (mm)")

subplot(413)
plot(times_plot,(r_save_plot[:,1]), label="x");
ylabel("x (mm)")

subplot(414)
plot(times_plot,mu_save_plot[:,1], label=L"$\mu_x$")
plot(times_plot,mu_save_plot[:,2], label=L"$\mu_y$")
plot(times_plot,mu_save_plot[:,3], label=L"$\mu_z$")
plot(times_plot,mu_norm_plot, label=L"$\mu$")
xlabel(L"time ($\mu$s)")
ylabel(L"magnetic moment ($\mu_B$)")
legend(loc="right");
```

Because we initialised the magnetic moment at $x=0$, where $B_x=0$, the precession started around the z direction and so we see an oscillation of the $\mu_x$ and $\mu_y$. Because the $\mu_x$ oscillates, the force from the gradient of B in x also oscillates which creates the ripple in the x position.


## Simulating many particles


There are many things we can change about the initial conditions when simulating many particles, we can create a distribution of positions, of speeds, of magnetic moments etc. Let's begin with creating a distribution of magnetic moments. This simply means running the code above many times, picking a different random $\mu$ every time.

Because most of the interest in the SG experiment is in the final z distribution of the particles after they have gone through the magnet, we'll store the final z position and visualise how they are distributed. 

```julia
# number of particles
np = 1000

# place to store the final z-position
zs = zeros(np)

for p in 1:np
    
    v = [0,v0,0]
    r = [0,0,0]
    theta = acos(2*rand()-1)
    phi = rand()*2*pi
    mu = [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)]*mu_Ag


    for (i,t) in enumerate(times)

        r,v,mu = frog_boris(r,v,mu)

    end
    
    zs[p] = r[3]

end
```

```julia
hist(zs/1e-3, normed=true);
xlabel("z (mm)");
```

We can see that there is a pretty even distribution of z's as one might expect - nothing controversial.


## Including edge effects


The França paper makes the case that the SG results can be explained classically as being caused by the inhomogeneity in the magnetic field along y. Specifically, the magnetic field goes from being zero outside the magnet to a maximum inside the magnet over a very short length scale (cf Fig 3 of França):

![B field gradient in y](SG-B-gradient-y.png)

França shows that if energy of the particle and electromagnet are considered together, then the work done by (or on) the magnet, as the particle enters, leads to a change of$\theta$ (the angle of $\mu$ with respect to z axis) according to:

$$
\dot \theta = - \frac{v_y}{B_z(y)}\frac{\partial B_z(y)}{\partial y} \cot(y)
$$


One can begin to think about trying to verify this through simulation by first creating some approximations for the edge field, e.g.

```julia
z = LinRange(-0.006,0.004,100)

subplot(211)
B_edge = B0.*(0.5.*tanh.(800.0.*z.+2.0).+0.5)
plot(z,B_edge)
ylabel(L"B_z (T)");

subplot(212)
B_edge_grad = B0.*0.5.*sech.(800.0.*z.+2.0).^2
plot(z,B_edge_grad)
xlabel("y (mm)");
ylabel(L"\partial B_z /\partial y \ (T/m)");
```

We can re-create the field and gradient functions to include the new inhomogeneity. 

```julia
function B(r)
    B = zeros(3)
    B[1] = -beta*r[1]
    B[3] = B0.*(0.5.*tanh.(800.0.*r[3].+2.0).+0.5);
    return B
end;
```

```julia
function gradB(r)
    gradB = zeros(3,3)
    gradB[3,3] = beta
    gradB[1,1] = -beta
    gradB[3,2] = B0.*0.5.*sech.(800.0.*r[3].+2.0).^2
    return gradB
end;
```

and then re-run the simulation code for many particles:

```julia
# number of particles
np = 1000

# place to store the final z-position
zs = zeros(np)

for p in 1:np
    
    v = [0,v0,0]
    r = [0,0,0]
    theta = acos(2*rand()-1)
    phi = rand()*2*pi
    mu = [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)]*mu_Ag


    for (i,t) in enumerate(times)

        r,v,mu = frog_boris(r,v,mu)

    end
    
    zs[p] = r[3]

end
```

```julia
hist(zs/1e-3, normed=true);
xlabel("z (mm)");
```

There is again nothing striking in the above figure, the distribution is again pretty uniform. We might wonder whether we have sufficient resolution of the inhomogeneity, if not then maybe that's why we don't see any difference. Let's see:

```julia
# time to traverse the Bz ramp
ramp_time = 4e-3 / v0
```

```julia
# how many time steps do we have while the particle is in the Bx ramp
ramp_time/dt
```

We seem to have enough resolution, so that's not the problem.


Perhaps the problem lies in our current equations of motion:
$$
\frac{dv}{dt} = \frac{1}{m} \mu \nabla B
$$


$$
\frac{d\mu}{dt} = \gamma \mu\times B
$$


There does not appear to be an obvious way that the inhomogeneity in the field  can affect $\mu$ in the way described by França. Our equations are however not complete - there are several options for modifying the equations of motion to include more physics. We could:
1. Include self consistent equations for the magnetic field - this would be in the spirit of the França paper.
2. Go beyond the infinitesimal magnetic dipole approximation, i.e. beyond $\tau = \mu \times B$. Rather than considering only the dipole contributions to the torque perhaps we could also look at the [quadrupole terms](https://physics.stackexchange.com/a/208922).
3. Add in relativistic effects arising from the accelerating motion of the particle, so called [Thomas precession](https://en.wikipedia.org/wiki/Thomas_precession).

We'll consider relativity first.


## Thomas precession

In addition to the precession of $\mu$ about the magnetic field that we've described in our equations of motion, there is another term called [Thomas precession](https://en.wikipedia.org/wiki/Thomas_precession). This precession comes from the fact that the particle, whose magnetic movement we are trying to describe, is accelerating. When the particle accelerates, its instantaneous rest frame turns out to rotate with angular velocity $\omega$ given approximately (for non-relativistic speeds) by:
$$
\omega \approx \frac{1}{2c^2}a\times v
$$
where $a$ is the acceleration and $v$ the velocity.

In rotating frames, the equations of motion need to be modified because the corodinate system itself is time-dependent. [Working through the details](https://www.southampton.ac.uk/~stefano/courses/PHYS2006/chapter4.pdf), one finds:
$$
\frac{d}{dt}_{rot} = \frac{d}{dt}_{rest} + \omega \times
$$
The equation of motion for $\mu$ therefore needs to be modified to give:
$$
\frac{d \mu}{dt} = \gamma\mu \times B +\frac{1}{2c^2}(a\times v)\times \mu
$$
Substituting the particle acceleration $\frac{dv}{dt} = \frac{1}{m} \mu \nabla B$, we get:
$$
\frac{d \mu}{dt} = \gamma \mu \times \left[ B -\frac{1}{2mc^2}(\mu\nabla B)\times v) \right]
$$


---
**To be continued**
