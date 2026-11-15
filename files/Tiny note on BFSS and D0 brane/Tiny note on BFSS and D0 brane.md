There are plenty of arXiv preprints regarding the solution of BFSS matrix model since last week, for example 2611.09837, 2611.14512, 2611.17325. This topic is quite interesting and increasingly important, therefore I'm planning to learn with this tiny note.

# What is BFSS

The Banks–Fischler–Shenker–Susskind matrix quantum mechanics is the unique interacting quantum mechanics with maximal supersymmetry and $Spin(9)_R$ obtaining from dimensional reduction of $10d\ \mathcal{N}=1$ theory, with d.o.f described by a $Spin(9)$ vector $X^I$ and a $Spin(9)$ spinor $\theta^{\alpha}$, both in $U(N)$ adjoint representation. This describes $N$ interacting D0 brane in $Type IIA$ string theory, dual to a KK momentum around the compactified circle.

Finding the threshold bound state is a problem of constructing the gauge-invariant $Spin(9)$ singlet ground state of the $SU(N)$ BFSS (after decoupling the $U(1)$ sector from $U(N)$). Since gapless, it's quite difficult.

The Hamiltonian is 
$$
H = \mathrm{Tr}\left[\frac{1}{2}P_iP_i-\frac{1}{4}[X^i,X^j]^2 - \frac{1}{2}\theta^T\gamma_i[X^i,\theta]\right]
$$
Since $V(X) = 0$ implies the flat direction $[X^i,X^j]=0$ , this means we actually can separate D0 brane to arbitrary far places without contributing to its potential energy, which exists scattering state has energy close to 0, making it gapless.
# 2611.09837

This paper solve the ground (bound) state wavefunction of $SU(2)$ BFSS matrix theory in a very high numerical precision. They use the Rayleigh-Ritz method to an auxiliary operator $H_{aux}$, specially they use representation theory for the d.o.f reduction specially in $SU(2)$.

According to our definition of the theory, it's not hard to obtain the Hilbert space of the theory. For bosonic part, we have a infinite dimensional space $L^2(\mathrm{R}^{27})$ so we denote it in the position picture $(X)$. For fermionic d.o.f, the Hilbert space is in finite dimension, we note it in Fock picture by defining
$$
f_j = \frac{\theta_{2j-1}+i\theta_{2j}}{\sqrt{2}},\quad j=1,\dots,24
$$
the total Hilbert space is
$$
L^2(\mathbb{R}^{27})\otimes \mathbb{C}^{2^{24}}
$$
which is too massive to compute.

They first apply the singular value decomposition $X = OD(x)V^{T}$ to the matrix $X_{ai},\quad a=1,2,3,\ i=1,\cdots,9$, this means there are only $x_1,x_2,x_3$ as singular values. Since the ground state as a gauge and $Spin(9)$ singlet, we should have
$$
| \Psi(OXV^T)\rangle = U(O,V)|\Psi(X)\rangle
$$
Moreover, they examine the stable subgroup maintaining the three singular values, the subgroup consists of the discrete rotation of first three directions and the rotation of last six directions, together
$$
G_D = K_4\times SO(6)\qquad K_4 = \mathbb{Z}_2\times \mathbb{Z}_2
$$
so fermion wave function belongs to the invariant subspace of this subgroup
$$
\mathcal{F}_* = \mathcal{F}^{K_4\times SO(6)}
$$
Using the Schur polynomial, the author obtain that
$$
\dim \mathcal{F}_* = 2066
$$
so the wavefunction can be denoted as
$$
|\Psi(D(x))\rangle = \sum_{\mu=1}^{2066}\phi_{\mu}(x_1,x_2,x_3)|e_{\mu}\rangle
$$
Why this reduction is strictly restricted to $SU(2)$? For general $SU(N)$ BFSS, we have the number of bosonic variables as $9(N^2-1)$, and the overall dimension of gauge group together with R symmetry group is $(N^2-1)+36 = N^2 + 35$, only when $N=2$ the first is lower then the latter, otherwise it's not possible to cancel most bosonic variables.

Normally for the Rayleigh-Ritz method, we construct a testing wavefunction $|\Psi\rangle$ and minimize
$$
E_{\Psi} = \frac{\langle\Psi|H|\Psi\rangle}{\langle\Psi|\Psi\rangle}
$$
However, this cannot work for BFSS since it has a continuous spectrum starting from $E=0$. The paper constructs an auxiliary Hamiltonian 
$$
H_{aux} = \frac{1}{16}\sum_{\alpha=1}^{16}Q_{\alpha}R^2Q_{\alpha}\qquad R^2=\sum_{a,i}X^2_{ai}
$$
This auxiliary Hamiltonian has the same zero energy sector with the original Hamiltonian, however, it is proved that this auxiliary Hamiltonian has a gap, which enables this method to work.
$$
\Delta_{aux}\ge \frac{9-2\sqrt{14}}{4}\approx 0.379
$$
How is this wonderful relation proved? 

First, a small example showing why introducing $R^2$ would possibly make the theory gapped. We all know $H_0=p^2$ is gapless, for $H = px^2p$, we can prove that
$$
0\le \int_0^{\infty}\left|x\psi^{\prime}(x)+\frac{1}{2}\psi(x)\right|^2dx = \int_0^{\infty}x^2\left|\psi^{\prime}(x)\right|^2dx-\frac{1}{4}\int_0^{\infty}\left|\psi(x)\right|^2dx
$$
which means adding such $R^2$ has possibilities to make the theory gapped. For $|\eta\rangle$ in $\mathbf{16}$ of $Spin(9)$, the paper proved
$$
\langle\eta|H|\eta\rangle \ge c\langle\eta|R^{-2}|\eta\rangle
$$
where $c=\frac{9-2\sqrt{14}}{4}$. Then defines $F=\langle\psi|P_{\ker H}|\psi\rangle$, the paper also proved
$$
1-F\le\frac{1}{c}\langle H_{aux}\rangle_{\psi}
$$
this means $\langle H_{aux}\rangle_{\psi}\ge c$.

After obtaining the wave function, the author checked the $\langle H\rangle\le 10^{-8}$ while $F\ge 0.9999978$ , since the $SU(2)$ BFSS has a unique zero energy ground state, this means this wave function has a very large overlap with it. Also, the author checked the spatial distribution of the wave function, by defining $r$ as the largest singular value of $X$, they checked 
$$
p(r)\sim 1890 \sqrt{2} r^{-10}
$$
which can be derived from the three graviton amplitude of M theory(2312.12592) .
# 2611.14512

For the other two papers, we focus on their difference between the first paper, especially numerical solving methods. Since they both used the same d.o.f reduction of the theory, the solution also becomes using the testing function
$$
|\Omega\rangle = \sum_{t=1}^{2066}f_t(s_1,s_2,s_3)|t\rangle
$$
It's worth noting that the singular value is slightly different with the last paper's choice, the author use a gauge fixing here regarding the gauge group $SU(2)$
$$
X^i = s_iT_i\quad i=1,2,3\qquad X^4=\cdots=X^9=0
$$
where the scalar potential is
$$
V(s_1,s_2,s_3) = s_1^2s_2^2+s_1^2s_3^2+s_2^2s_3^2
$$

However, the functional method use a different function to minimize in this paper
$$
\frac{\mathcal{F}(\Omega)}{\mathcal{N}(\Omega)}\quad\mathcal{F}(\Omega) = \sum_{\alpha=1}^{16}\langle\Omega|Q_{\alpha}w(s)Q_{\alpha}|\Omega\rangle\quad \mathcal{N}(\Omega)=\langle\Omega|\Omega\rangle
$$
where
$$
w(s) = \left(1+\frac{s_1^2}{c_w}\right)^2, \qquad c_w=14.29
$$
the testing wave function
$$
\Omega(s)= \frac{1}{(1+s^2/c)^3} \exp\left( -\frac{V(s)}{\sqrt{2s^2}} \right) \sum_t\chi_t(s)|t\rangle\quad s^2=s_1^2+s_2^2+s_3^2,\qquad c=\frac{20}{2^{4/3}}
$$
on actual calculation, the paper separate the $s_1$ direction to core and tail part and do polynomial interpolation to each part, with $8\times 10^{5}$ nodes and $1.7\times 10^9$ variables.

Not quite interested in the scattering analysis afterwards, skip.
