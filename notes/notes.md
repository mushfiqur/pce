Set of legendre polynomials:

$$ f_0(x) = 1 $$
$$ f_1(x) = x $$
$$ f_2(x) = \frac{1}{3}(3x^2 - 1)$$
$$ ... $$

These can be generated up to any order with this relationship:

$$
(n+1)f_{n+1} (x) = (2n + 1)xf_n(x) - nf_{n-1}(x)
$$

We need multivariate a legendre polynomial basis set.

To calculate this set, generate $n$ sets of polynomials, where $n$ is the number of noise sources.

$$ A =  \{1, x_1, \frac{1}{3}(3{x_1}^2 - 1)\}$$
$$ B = \{1, x_2, \frac{1}{3}(3{x_2}^2 - 1)\}$$
$$...$$

Then calculate the tensor product of the sets:

$$P = A ⊗ B$$

The resulting set is:
$$ P_0 = 1 \cdot 1 $$
$$ P_1 = 1 \cdot x_2 $$
$$ P_2 = 1 \cdot \frac{1}{3}(3{x_2}^2 - 1) $$
$$ P_3 = x_1 \cdot 1 $$
$$ P_4 = x_1 \cdot x_2 $$
$$ P_5 = x_1 \cdot \frac{1}{3}(3{x_2}^2 - 1) $$
$$ P_6 = \frac{1}{3}(3{x_1}^2 - 1) \cdot 1 $$
$$ P_7 = \frac{1}{3}(3{x_1}^2 - 1) \cdot x_2 $$
$$ P_8 = \frac{1}{3}(3{x_1}^2 - 1) \cdot \frac{1}{3}(3{x_2}^2 - 1) $$

Then, keep only the polynomials up to a specified order $p$
$$ P_0 = 1 $$
$$ P_1 =  x_2 $$
$$ P_2 = \frac{1}{3}(3{x_2}^2 - 1) $$
$$ P_3 = x_1 $$
$$ P_4 = x_1 x_2 $$
$$ P_5 = \frac{1}{3}(3{x_1}^2 - 1) $$

Lastly, we need to calculate the expected values $E[P_iP_jP_k]$ where $i,j,k ∈ \{0,1,2,3,4,5\}$. And also, $E[{P_0}^2],E[{P_1}^2],E[{P_2}^2],E[{P_3}^2],E[{P_4}^2],E[{P_5}^2]$

Total number of terms for $n$ noise sources up to order $p$ is:

$$
\frac{(n+p)!}{n! \cdot p!}
$$

For $p=2$:
$$
\frac{(n+1)(n+2)}{2}
$$

With a 100 noise sources, the polynomial basis set has $5151$ polynomials and generating the look up table requires $45,583,355,552$ polynomial multiplications. 

$$ P_1 = {x_3}$$
$$ P_8 = {x_1}{x_2}$$
$$ $$
$$ E[P_1 P_8 P_k] = {x_3} \cdot \ {x_1}{x_2} \cdot P_k$$

$$ E[\Psi_0\cdot\Psi_1\cdot\Psi_3] $$

$$
(\frac{1}{\alpha^{n+1} - 1})^2 \cdot \frac{1}{3} 
$$