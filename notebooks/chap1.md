# 1.4 Computing Actions


```python
from sympy import symbols, Rational, integrate, Function
from sympy import Matrix

x = Function('x')
y = Function('y')
z = Function('z')
t = symbols('t')
m = symbols('m')

q = Matrix([x(t), y(t), z(t)])

def velocity(local):
    return Matrix(local[4:7])


def l_free_particle(mass):
    def inner(local):
        v = velocity(local)
        return Rational(1,2) * mass * (v.transpose() * v)[0]
    
    return inner


def gamma(q, t):
    # Returns the local tuple
    return Matrix([t, q, q.diff()])


def lagrangian_action(L, q, t1, t2):
    return integrate(L(Gamma(q,t)), (t, t1, t2))


```


```python
test_path = Matrix([4 * t + 7, 3 * t + 5, 2 * t + 1])
test_path
```




$\displaystyle \left[\begin{matrix}4 t + 7\\3 t + 5\\2 t + 1\end{matrix}\right]$




```python
l_free_particle(3.0)(gamma(q,t))

```




$\displaystyle 1.5 \left(\frac{d}{d t} x{\left(t \right)}\right)^{2} + 1.5 \left(\frac{d}{d t} y{\left(t \right)}\right)^{2} + 1.5 \left(\frac{d}{d t} z{\left(t \right)}\right)^{2}$




```python
lagrangian_action(l_free_particle(3.0), test_path, 0.0, 10.0)
```




$\displaystyle 435.0$




```python

```
