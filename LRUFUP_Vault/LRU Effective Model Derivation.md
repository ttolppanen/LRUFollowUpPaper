We would like to simulate the 5-qubit code with our LRU attached, but this would total to 15 qutrits. To reduce the numerical complexity, we could try to make an effective model of our LRU, tracing out the degrees of freedom related to the unwanted sites.
	Tracing out the second and the third site directly would give us just a model where we do not
see the dissipation on the last site, and our effective LRU would not work. What we can do is to change the basis for the operators, so that
$$
\hat{H} = \sum_i w_i \hat{n}_i +  J(\hat{a}_i\hat{a}^{\dagger}_{i + 1} + h.c.) = \sum_i\omega_i\hat{b}^{\dagger}_i\hat{b}_i.
$$
We will go to a new basis, in which there is no interaction between the "sites". Now in the new basis tracing out the extra degrees of freedom is more justified.
### Splitting the Qutrit Hamiltonian In Two
Note that the above Hamiltonian describes coupled Harmonic oscillators. I think for the Bose-Hubbard model we would not be able to do this due to the anharmonicity. Anyways, we can justify working with two level systems by considering the stack dynamics. For large enough anharmonicity, we know that the dynamics are such that they conserve the number of stacks. This means that a stack of excitation will not split into single excitations, or that single excitations will not form a stack. From this it follows that for large U the Bose-Hubbard Hamiltonian is block diagonal for different number of stacks. It means that there is no interaction between stacks and single excitations, so we can consider two separate systems instead. We will first derive the effective model for single excitations, working with the two level Bose-Hubbard model. Then we will derive the effective model for the stack dynamics, starting from the effective stack model.
### Finding the Basis
Now we will derive the new lowering operators $\hat{b}_i$. To start, we can write the Hamiltonian as
$$
\hat{H} = \sum_i w_i \hat{n}_i +  J(\hat{a}_i\hat{a}^{\dagger}_{i + 1} + h.c.) 
=
\begin{pmatrix}
\hat{a}_1^\dagger & \hat{a}_2^\dagger & \hat{a}_3^\dagger
\end{pmatrix}

\begin{pmatrix}
w_1 & J   & 0   \\
J   & w_2 & J   \\
0   & J   & w_3
\end{pmatrix}

\begin{pmatrix}
\hat{a}_1 \\
\hat{a}_2 \\
\hat{a}_3
\end{pmatrix}
= \vec{a}^{\dagger} A \vec{a},
$$
and we would like to diagonalize the matrix $A$. The matrix $A$ is little bit complicated, and the solution for it is involved, so for now we will make some simplifications and work with that. We express the Hamiltonian in the units of $J$, so that we can set $J = 1$. Furthermore, we assume that the qubits have fixed frequency difference so that $w_1 = w + \Delta$, $w_2 = w$ and $w_3 = w - \Delta$. Then the matrix $A$ is
$$
A = \begin{pmatrix}
w + \Delta & 1   & 0   \\
1   & w & 1   \\
0   & 1   & w - \Delta
\end{pmatrix}.
$$
We can use mathematica to solve the eigenvectors and eigenvalues. From the eigenvectors we can construct the unitary matrix $U$ (**missing normalization**) that will diagonalize the matrix $A$
$$
U = 
\begin{pmatrix}
\vec{e}_1 & \vec{e}_2 & \vec{e}_3  
\end{pmatrix}
= 
\begin{pmatrix}
1 + \Delta \tilde{\Delta}_+ & -1       & 1 + \Delta \tilde{\Delta}_- \\
\tilde{\Delta}_+            & \Delta   & \tilde{\Delta}_-            \\
1                           & 1        & 1
\end{pmatrix},
\quad 
\tilde{\Delta}_{\pm} = \Delta \pm \sqrt{2 + \Delta^2}.
$$
The lengths of the eigenvectors needed for normalization are 
$$
\begin{align}
|\vec{e}_1| &= \sqrt{2}\sqrt{(2 + \Delta^2)(1 + \Delta \tilde{\Delta}_+}), \\
|\vec{e}_2| &= \sqrt{2 + \Delta^2}, \\ 
|\vec{e}_3| &= \sqrt{2}\sqrt{(2 + \Delta^2)(1 + \Delta \tilde{\Delta}_-}). 
\end{align}
$$
We can now diagonalize the matrix $A$
$$
UAU^\dagger = D =
\begin{pmatrix}
w + \sqrt{2 + \Delta^2} & 0   & 0   \\
0   & w & 0   \\
0   & 0   & w - \sqrt{2 + \Delta^2}
\end{pmatrix}.
$$
This implies that
$$
\vec{a}^{\dagger} A \vec{a} = \vec{a}^{\dagger} U^\dagger U A U^\dagger U \vec{a}
= \vec{b}^{\dagger} D \vec{b},
$$
where
$$
\vec{b} = U \vec{a} 
= 
\begin{pmatrix}
\frac{1 + \Delta \tilde{\Delta}_+}{|\vec{e}_1|} & \frac{-1}{|\vec{e}_2|}       & \frac{1 + \Delta \tilde{\Delta}_-}{|\vec{e}_3|} \\
\frac{\tilde{\Delta}_+}{|\vec{e}_1|}            & \frac{\Delta}{|\vec{e}_2|}   & \frac{\tilde{\Delta}_-}{|\vec{e}_3|}           \\
\frac{1}{|\vec{e}_1|}                            & \frac{1}{|\vec{e}_2|}        & \frac{1}{|\vec{e}_3|}
\end{pmatrix}
\begin{pmatrix}
\hat{a}_1 \\
\hat{a}_2 \\
\hat{a}_3
\end{pmatrix}
=
\begin{pmatrix}
\frac{1 + \Delta \tilde{\Delta}_+}{|\vec{e}_1|} \hat{a}_1 + \frac{-1}{|\vec{e}_2|} \hat{a}_2 + \frac{1 + \Delta \tilde{\Delta}_-}{|\vec{e}_3|} \hat{a}_3 \\
\frac{\tilde{\Delta}_+}{|\vec{e}_1|} \hat{a}_1 + \frac{\Delta}{|\vec{e}_2|} \hat{a}_2 + \frac{\tilde{\Delta}_-}{|\vec{e}_3|} \hat{a}_3 \\
\frac{1}{|\vec{e}_1|} \hat{a}_1 + \frac{1}{|\vec{e}_2|} \hat{a}_2 + \frac{1}{|\vec{e}_3|} \hat{a}_3
\end{pmatrix}.
$$
We can also show that
$$
\vec{a} = U^\dagger \vec{b} 
=
\begin{pmatrix}
\frac{(1 + \Delta \tilde{\Delta}_+)\hat{b}_1 + \tilde{\Delta}_+ \hat{b}_2 + \hat{b}_3}{|\vec{e}_1|} \\
\frac{-\hat{b}_1 + \Delta \hat{b}_2 + \hat{b}_3}{|\vec{e}_2|} \\
\frac{(1 + \Delta \tilde{\Delta}_-)\hat{b}_1 + \tilde{\Delta}_- \hat{b}_2 + \hat{b}_3}{|\vec{e}_3|}
\end{pmatrix}.
$$
We will need $\hat{a}_3$ from this expression later.

Now we have what we need to diagonalize the Hamiltonian, it simply reads as
$$
\hat{H} = \sum_i w_i \hat{n}_i +  J(\hat{a}_i\hat{a}^{\dagger}_{i + 1} + h.c.) 
= \vec{a}^{\dagger} A \vec{a}
= \vec{b}^{\dagger} D \vec{b}
= (w + \sqrt{2 + \Delta^2}) \hat{b}^\dagger_1\hat{b}_1 + w \hat{b}^\dagger_2\hat{b}_2 + (w - \sqrt{2 + \Delta^2}) \hat{b}^\dagger_3\hat{b}_3
$$
We can test our calculation, by noting that when the system is uncoupled to begin with, the new operators should reduce to the old. This means that the unitary matrix $U$ should become diagonal. We take the limit of $J \to 0$ , since everything is expressed in the units of $J$, if follows that $\Delta \to \infty$. We can take the limit with mathematica, it gives
$$
\begin{align}
\vec{e}_1 / |\vec{e}_1| &\xrightarrow{\Delta \to \infty} (1, 0, 0)^T, \\
\vec{e}_2 / |\vec{e}_2| &\xrightarrow{\Delta \to \infty} (0, 1, 0)^T, \\
\vec{e}_3 / |\vec{e}_3| &\xrightarrow{\Delta \to \infty} (0, 0, 1)^T. \\
\end{align}
$$
I checked some of these, but deriving them here is kinda lengthy and besides the point. Note here that as $J \to 0$, $\hat{b}_1 \to \hat{a}_1$, so that $\hat{b}_1$ is describes the first site.
### New Operators Describe Harmonic Oscillators
To show that the new operators $\hat{b}_i$ indeed describe harmonic oscillators, they need to fulfill the [commutation relation](https://en.wikipedia.org/wiki/Ladder_operator#General_formulation)
$$
[\hat{b}^\dagger_i, \hat{b}_i] = -1.
$$
To show this, notice that the unitary matrix $U$ has orthonormal column and **row** vectors, following directly from it being unitary. With this in mind we can write
$$
\vec{b} = U\vec{a}
=\begin{pmatrix}
\vec{r}_1 \cdot \vec{a} \\
\vec{r}_2 \cdot \vec{a} \\
\vec{r}_3 \cdot \vec{a}
\end{pmatrix}.
$$
For the commutation relation, we need to calculate $\hat{b}^\dagger_i \hat{b}_i$ and  $\hat{b}_i \hat{b}^\dagger_i$.
$$
\hat{b}^\dagger_i \hat{b}_i = (\vec{r}_i \cdot \vec{a}^\dagger)(\vec{r}_i \cdot \vec{a})
=
(\sum_{j} r_{ij}\hat{a}^\dagger_j)(\sum_k r_{ik}\hat{a}_k)
=
\sum_{jk} r_{ij} r_{ik} \hat{a}^\dagger_j \hat{a}_k.
$$
Similarly we get
$$
\hat{b}_i \hat{b}^\dagger_i
=
\sum_{jk} r_{ij} r_{ik} \hat{a}_j \hat{a}^\dagger_k.
$$
Before we write $[\hat{b}^\dagger_i, \hat{b}_i]$, note that since
$$
[\hat{a}_i^\dagger, \hat{a}_j] = 0 \quad \text{if } i \neq j
$$
then in the commutator $[\hat{b}^\dagger_i, \hat{b}_i]$ all terms where $i \neq j$ will vanish
$$
\begin{align}
[\hat{b}^\dagger_i, \hat{b}_i]
&=
\hat{b}^\dagger_i \hat{b}_i - \hat{b}_i \hat{b}^\dagger_i
=
\sum_{jk} r_{ij} r_{ik} \hat{a}^\dagger_j \hat{a}_k - \sum_{jk} r_{ij} r_{ik} \hat{a}_j \hat{a}^\dagger_k
=
\sum_{jk} r_{ij} r_{ik} (\hat{a}^\dagger_j \hat{a}_k - \hat{a}_k \hat{a}^\dagger_j) \\
&=
\sum_{jk} r_{ij} r_{ik} [\hat{a}^\dagger_j, \hat{a}_k]
=
\sum_{j} r_{ij} r_{ij} [\hat{a}^\dagger_j, \hat{a}_j]
=
\sum_{j} r_{ij} r_{ij} (-1)
=
\vec{r}_i \cdot \vec{r}_i (-1)
=
-1,
\end{align}
$$
since the row vectors $\vec{r}_i$ are orthonormal due to the unitary of $U$, $[\hat{a}^\dagger_i, \hat{a}_i] = -1$, and also
$$
\sum_{jk} r_{ij} r_{ik} \hat{a}_j \hat{a}^\dagger_k
=
\sum_{kj} r_{ik} r_{ij} \hat{a}_k \hat{a}^\dagger_j
=
\sum_{jk} r_{ij} r_{ik} \hat{a}_k \hat{a}^\dagger_j.
$$
So, as long as $U$ is unitary, the operators $\hat{b}_i$ are valid lowering operators, and describe a harmonic oscillator.
### New Basis
work in progress...
With the new operators comes a new basis. It would be insightful to know the connection between the original basis and the new basis. To derive the new basis, we can calculate the eigenstates of the new operators. With mathematica we get
$$
\begin{align}
|0\rangle_b &= | 0\rangle + ...,\\
|1\rangle_b &= .
\end{align}
$$
$$
\langle \Psi| a^\dagger a|\Psi\rangle = \langle \Psi| U^\dagger U a^\dagger a U^\dagger U|\Psi\rangle = \langle \Psi_b| b^\dagger b|\Psi_b\rangle
$$
##### Excitations in the New Operators
### Effective Master Equation
To derive the effective model, we start from the LRU master equation where instead of measurement we use dissipation
$$
\frac{d\hat\rho}{dt} = -i[\hat{H}, \hat{\rho}] + \Gamma \left( \hat{a}_3\rho\hat{a}_3^\dagger - \frac{1}{2}\{{\hat{a}_3^\dagger\hat{a}_3, \hat{\rho}}\} \right).
$$
Since we know the Hamiltonian in the new basis, we need to only focus on the dissipation part. We assume that in the new basis the initial state is not entangled, and since the Hamiltonian in the basis does not have interaction, the state will stay un-entangled. The initial state is thus
$$
\hat{\rho} = \hat{\rho}_1 \otimes \hat{\rho}_2 \otimes \hat{\rho}_3
$$
Now we should insert
$$
\hat{a}_3 = \frac{(1 + \Delta \tilde{\Delta}_-)\hat{b}_1 + \tilde{\Delta}_- \hat{b}_2 + \hat{b}_3}{\sqrt{2}\sqrt{(2 + \Delta^2)(1 + \Delta \tilde{\Delta}_-})} = \alpha\hat{b}_1 + \beta\hat{b}_2 + \gamma\hat{b}_3,
$$
to the dissipation term, and take the partial trace
$$
\frac{d\hat{\rho}}{dt} \rightarrow Tr_{2,3}(\frac{d\hat{\rho}}{dt}).
$$
Note that
$$
Tr_{2, 3} (\hat{\rho}\hat{M}_2) = \hat{\rho}_1 \otimes Tr_{2}(\hat{\rho}_2\hat{M}_2) \otimes Tr_{3}(\hat{\rho}_3) = \langle \hat{M}_2 \rangle\hat{\rho}_1. 
$$
First let's look at $\hat{a}_3^\dagger \hat{a}_3$
$$
\begin{align}
\hat{a}_3^\dagger \hat{a}_3 =& (\alpha\hat{b}^\dagger_1 + \beta\hat{b}^\dagger_2 + \gamma\hat{b}^\dagger_3) (\alpha\hat{b}_1 + \beta\hat{b}_2 + \gamma\hat{b}_3) \\

=& \alpha^2 \hat{b}^\dagger_1 \hat{b}_1 + \alpha\beta \hat{b}^\dagger_1 \hat{b}_2 + \alpha\gamma \hat{b}^\dagger_1 \hat{b}_3 \\

+& \alpha\beta \hat{b}^\dagger_2 \hat{b}_1 + \beta^2 \hat{b}^\dagger_2 \hat{b}_2 + \beta\gamma \hat{b}^\dagger_2 \hat{b}_3 \\

+& \alpha\gamma \hat{b}^\dagger_3 \hat{b}_1 + \beta\gamma \hat{b}^\dagger_3 \hat{b}_2 + \gamma^2 \hat{b}^\dagger_3 \hat{b}_3,
\end{align}
$$
so that
$$
\begin{align}
Tr_{2,3}(\hat{a}_3^\dagger\hat{a}_3\hat{\rho})

=& (\alpha^2 \hat{b}^\dagger_1 \hat{b}_1 \\

+& (\alpha\beta \langle \hat{b}_2 \rangle + \alpha\gamma \langle \hat{b}_3 \rangle )\hat{b}^\dagger_1 \\

+& (\alpha\beta\langle \hat{b}^\dagger_2 \rangle  + \alpha\gamma\langle \hat{b}^\dagger_3 \rangle)\hat{b}_1 \\

+& \beta^2 \langle\hat{b}^\dagger_2 \hat{b}_2\rangle + \beta\gamma \langle\hat{b}^\dagger_2\rangle \langle\hat{b}_3\rangle + \beta\gamma \langle\hat{b}^\dagger_3\rangle \langle\hat{b}_2\rangle + \gamma^2 \langle\hat{b}^\dagger_3 \hat{b}_3\rangle )\hat{\rho},
\end{align}
$$
Now for the first term we get
$$
\begin{align}
\hat{a}_3 \hat{\rho}\hat{a}^\dagger_3
=& \alpha^2 \hat{b}_1 \hat{\rho} \hat{b}^\dagger_1 + \alpha\beta \hat{b}_1 \hat{\rho} \hat{b}^\dagger_2 + \alpha\gamma \hat{b}_1 \hat{\rho} \hat{b}^\dagger_3 \\

+& \alpha\beta \hat{b}_2 \hat{\rho} \hat{b}^\dagger_1 + \beta^2 \hat{b}_2 \hat{\rho} \hat{b}^\dagger_2 + \beta\gamma \hat{b}_2 \hat{\rho} \hat{b}^\dagger_3 \\

+& \alpha\gamma \hat{b}_3 \hat{\rho} \hat{b}^\dagger_1 + \beta\gamma \hat{b}_3 \hat{\rho} \hat{b}^\dagger_2 + \gamma^2 \hat{b}_3 \hat{\rho} \hat{b}^\dagger_3,
\end{align}
$$
so that
$$
\begin{align}
Tr_{2, 3}(\hat{a}_3 \hat{\rho}\hat{a}^\dagger_3)
=& \alpha^2 \hat{b}_1 \hat{\rho} \hat{b}^\dagger_1 \\ 

+& (\alpha\beta \langle \hat{b}^\dagger_2 \rangle + \alpha\gamma  \langle \hat{b}^\dagger_3 \rangle)\hat{b}_1 \hat{\rho} \\

+& (\alpha\beta \langle \hat{b}_2 \rangle  + \alpha\gamma \langle \hat{b}_3 \rangle)\hat{\rho} \hat{b}^\dagger_1\\

+& (\beta^2 \langle \hat{b}^\dagger_2 \hat{b}_2 \rangle + \beta\gamma \langle \hat{b}_2 \rangle \langle \hat{b}^\dagger_3 \rangle + \beta\gamma \langle \hat{b}_3 \rangle \langle \hat{b}^\dagger_2 \rangle + \gamma^2 \langle \hat{b}^\dagger_3 \hat{b}_3 \rangle) \hat{\rho},
\end{align}
$$
Now using these we can calculate
$$
\begin{align}
Tr_{2,3} \left( \hat{a}_3\rho\hat{a}_3^\dagger - \frac{1}{2}\{{\hat{a}_3^\dagger\hat{a}_3, \hat{\rho}}\} \right) 
=& 
\alpha^2 \hat{b}_1 \hat{\rho} \hat{b}^\dagger_1 - \frac{1}{2} \alpha^2 \{\hat{b}^\dagger_1 \hat{b}_1 , \hat{\rho} \} \\

+& (\alpha\beta \langle \hat{b}^\dagger_2 \rangle + \alpha\gamma  \langle \hat{b}^\dagger_3 \rangle)\hat{b}_1 \hat{\rho} 

- \frac{1}{2}(\alpha\beta\langle \hat{b}^\dagger_2 \rangle  + \alpha\gamma\langle \hat{b}^\dagger_3 \rangle)(\hat{b}_1\hat{\rho} + \hat{\rho}\hat{b}_1) \\

+& (\alpha\beta \langle \hat{b}_2 \rangle  + \alpha\gamma \langle \hat{b}_3 \rangle)\hat{\rho} \hat{b}^\dagger_1 

- \frac{1}{2} (\alpha\beta \langle \hat{b}_2 \rangle + \alpha\gamma \langle \hat{b}_3 \rangle ) (\hat{b}^\dagger_1\hat{\rho} + \hat{\rho}\hat{b}^\dagger_1) \\

=& \alpha^2 \hat{b}_1 \hat{\rho} \hat{b}^\dagger_1 - \frac{1}{2} \alpha^2 \{\hat{b}^\dagger_1 \hat{b}_1 , \hat{\rho} \} \\

+& (\alpha\beta \langle \hat{b}^\dagger_2 \rangle + \alpha\gamma  \langle \hat{b}^\dagger_3 \rangle)[\hat{b}_1, \hat{\rho}] \\

+& (\alpha\beta \langle \hat{b}_2 \rangle  + \alpha\gamma \langle \hat{b}_3 \rangle)[\hat{\rho}, \hat{b}^\dagger_1].
\end{align}
$$
We can assume that the terms with the expectation values vanish. This in general depends on the initial state, but we can assume that it is so here. Now comparing to the original dissipation term, we see that dissipation on the last physical site amounts to dissipation with $\hat{b}_1$, where the dissipation rate has changed

$$
\Gamma \rightarrow \alpha^2 \Gamma, \quad \alpha^2 = \frac{1+\Delta \tilde\Delta_-}{2(2 + \Delta^2)} = \frac{1+\Delta^2 - \Delta \sqrt{2 + \Delta^2}}{2(2 + \Delta^2)}
$$
Here is a plot for $\alpha^2$ as a function of $\Delta$.
![[alpha_squared.png]]
In our paper we have $\Delta \approx 20$, so the dissipation vanishes on the first site very quickly. This is for the qubit, so it means that our first site is quite decoupled from the dissipation. With mathematica we can do a series expansion for $\alpha^2$ assuming that $\Delta / J >> 0$
$$
\alpha^2 = \frac{1}{4\Delta^4} + O[\frac{1}{\Delta}]^6 = \frac{J^4}{4\Delta^4} + O[\frac{J}{\Delta}]^6.
$$
To highlight, the point here is that **in a system of three transmons, the last site dissipation transforms to dissipation on a single system, where the dissipation rate is reduced by a factor of $\alpha^2$**.