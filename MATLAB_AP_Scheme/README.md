Author: [Hélène Hivert](http://perso.ec-lyon.fr/helene.hivert/)
Here is a short summary of the function:

birth, death, tau, m (respectively) : birth rate, death rate, tau function for gene transfer, gaussian kernel (may be useless in this last version of my codes - I am not sure).

munNeumann_x, punNeumann_y : computes u_{i-1} (mun=minus one) and u_{i+1} (pun=plus one) with "Neumann" boundary conditions (not really Neumann though here : the value at the boundary is chosen according the slope at the end of the interval).

u_init : the initial data for u

SchemaAP_Essai2 : the AP scheme. I think it works correctly for $\varepsilon \in [10^{-6}, 1]$. I am not sure for $10^{-7}$, but it is close to the square root of the numerical precision, so it's hard to say.
SchemaExplicite_f : an explicit scheme for f. Not robust when $\varepsilon\to 0$.
SchemaExplicite_u : an explicit scheme for u. The code is such that it is as robust as possible when $\varepsilon\to 0$. But it does not enjoy the AP property.
SchemaLimite_Implicite_u : the scheme for the limit equation, with the implicit treatment of $\rho$.

All the Test_*... are scripts calling the corresponding function for numerical tests.
