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
    display_name: Python 3
    language: python
    name: python3
---

# Stern and Gerlach simulation


It was discovered experimentally,[in 1922 by Stern and Gerlach](https://www.feynmanlectures.caltech.edu/II_35.html#Ch35-S2) (SG), that the spin of an electron along any direction is quantised, taking values of either $\hbar/2$ or $-\hbar/2$. 

Since then, people have described an idealised system of many SG systems in parallel that are capable of splitting apart and recombining the particle beams together (often called a Stern-Gerlach interferometer) to create interference effects that are commonly discussed in the two slit experiment (see e.g. the [Feynman lectures](https://www.feynmanlectures.caltech.edu/III_05.html#Ch5-S4)). 

To my surprise, it is only [very recently](https://arxiv.org/abs/1801.02708) that Stern-Gerlach interferometer experiments are being considered practical enough to implement and yet they are used frequently as a pedagogical example. It is therefore worth first developing an intuition for what one might expect from such SG interferometer experiments using only classical physics. That is the purpose of this notebook.




```python
import numpy as np
import sys
import matplotlib.pyplot as plt
%matplotlib inline
np.set_printoptions(threshold=sys.maxsize)
import scipy.constants as const
```

```python
const.find("dipole")
```

```python
machine_dim_x = 1.0
machine_dim_y = 0.1
machine_dim_z = 0.1

Bmin = 1.0e-4
Bmax = 10.0e-4
deltaB = Bmax - Bmin
gradB0 = deltaB / machine_dim_z

e = const.e
m = const.m_e
gyro =  e/(2.0*m)
mu_e = const.physical_constants["atomic unit of mag. dipole mom."][0]
```

```python
def B(r):
    B = np.zeros(3)
    B[2] = (Bmin + deltaB*r[2]/machine_dim_z)
    return B
```

```python
def gradB(r):
    gradB = np.zeros((3,3))
    gradB[2,2] = gradB0
    return gradB
```

```python
dt = 1.0e-10
t_max = 1.0e-5
times = np.arange(0,t_max,dt)

mu = np.array([0.0,1/np.sqrt(2),1/np.sqrt(2)])*mu_e
r = np.array([0,0.05,0.05])
v = np.array([600.0,0,0])
r_save = np.zeros((times.size,3))
mu_save = np.zeros((times.size,3))
```

```python
# evolution for mu inspired by https://www.particleincell.com/2011/vxb-rotation/
          
for i,t in enumerate(times):
    v += 0.5*dt*np.matmul(mu,gradB(r))/m            
    Br = B(r)
    s = 2.0/(1+(np.linalg.norm(Br)*gyro*dt*0.5)**2.0)
    mu_prime = mu + 0.5*dt*gyro*np.cross(mu,Br)
    mu += 0.5*dt*gyro*s*np.cross(mu_prime,Br)            
    r += dt*v
    v += 0.5*dt*np.matmul(mu,gradB(r))/m
    
    r_save[i,:] = r
    mu_save[i,:] = mu
```

```python
# plt.plot(times,r_save[:,0], label="x")
# plt.plot(times,r_save[:,1], label="y")
plt.plot(times,r_save[:,2], label="z")
plt.legend();
```

```python
plt.plot(times,mu_save[:,0], label="$\mu_x$")
plt.plot(times,mu_save[:,1], label="$\mu_y$")
plt.plot(times,mu_save[:,2], label="$\mu_z$")
plt.legend(loc="right");
```

```python
r[0]
```

```python

```

```python

```

## Adding inhomogeneous field

```python
s = 0.5*np.tanh(-10.0*(x-machine_dim_x))+0.5
```

```python
plt.plot(x,s)
```

```python

```
