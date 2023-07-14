The support for the standard definition of Legendre polynomials is the range $[-1,1]$. So whent the standard Legendre polynomials are used, the underlying random variable is assumed to have the interval $[-1,1]$, not $[a,b]$.

If we would like to represent a uniform random variable on $[a,b]$ and not $[-1,1]$, then we need to find a new set of Legendre polynomials that is "scaled" to be orthogonal on $[a,b]$.

If we have a random variable $x$ which is on $[a,b]$, we have to map it to a random variable $z$ which is on $[-1,1]$. This is done by the equation

$$ x = \frac{b-a}{2}z + \frac{a+b}{2} $$

The expansion of a function on the interval $[a,b]$ in Legendre polynomials is given by:

$$g(x) = \sum_{n=0}^{∞} c_n \cdot P_n(z),	x∈[a,b]$$

and 

$$ c_n = \frac{2n+1}{2} \int_{-1}^{1} g(\frac{b-a}{2}z + \frac{a+b}{2})P_n(z)\mathrm{d}z$$

If we consider $g(x) = x, x \sim U(-0.9,0.9)$ then:

$$x =\frac{b-a}{2}z + \frac{a+b}{2} = \frac{1.8}{2}z $$

$$ c_0 = \frac{1}{2} \int_{-1}^{1} \frac{1.8}{2}zP_0(z)\mathrm{d}z  $$

$$ c_1 = \frac{3}{2} \int_{-1}^{1} \frac{1.8}{2}zP_1(z)\mathrm{d}z = 1.8/2 $$

$$ ... $$

And we get:
$$ x = \frac{1.8}{2}z, z \sim U(-1,1) $$

So we can write x as:
$$ x = x_0\Psi_0 + x_1\Psi_1 + x_2\Psi_2 + ... $$
Where $x_0 = 0, x_1 = \frac{1.8}{2}, x_2 = 0, x_3 = 0 ...$

