A.Kitaev and C.Laumann


# 2 Topological phenomena in 1D: boundary modes in the Majorana chain

## 2.1 Definition
这章主要讨论对象是两个有$\mathbb{Z}_2$对称性并有topological degeneracy的例子，分别是transverse field Ising model（TFIM）以及spin-polarized superconductor（SPSC）。这两个模型的定义如下

- TFIM的定义是一个spin-$\frac{1}{2}$模型，其Hamiltonian为
$$
H_S = -J\sum_{j=1}^{N-1}\sigma^x_j\sigma^x_{j+1}-h_z\sum_{j=1}^N\sigma^z_j
$$
其中我们可以定义global symmetry算符
$$
U = \prod_{j=1}^N \sigma^z_j
$$
利用$\sigma^z\sigma^x\sigma^z = -\sigma^x$，我们显然可以使得
$$
[H,U]=0\quad U^2=I
$$

- SPSC模型的定义是一个1D fermionic system
$$
H_F = \sum_{j=1}^{N-1}\left(-w(a_j^{\dagger}a_{j+1}+a_{j+1}^{\dagger}a_j)+\Delta a_ja_{j+1}+\Delta^*a_{j+1}^{\dagger}a_j^{\agger}\right)-\mu\sum_{j=1}^N\left(a_j^{\dagger}a_j-\frac{1}{2}\right)
$$
我们一般会假设$w = \Delta = \Delta^*$，则
$$
H_F = w\sum_{j=1}^{N-1}(a_j-a_j^{\dagger})(a_{j+1}+a_{j+1}^{\dagger})-\mu\sum_{j=1}^N(a_j^{\dagger}a_j-\frac{1}{2})
$$
其中我们可以定义
$$
P_F = (-1)^{\sum_ja_j^{\dagger}a_j}
$$
利用$Pa_jP^{-1}=-a_j\quad Pa_j^{\dagger}P^{-1}=-a_j^{\dagger}$，得到
$$
[H,P]=0\quad P^2=I
$$

 对于TFIM，在$h_z=0$的时候，基态一般写为
$$
|\leftarrow\cdots \leftarrow\rangle\quad |\rightarrow\cdots\rightarrow\rangle
$$
其中$\sigma^x |\leftarrow\rangle = |\leftarrow\rangle \quad \sigma^x |\rightarrow\rangle = -|\rightarrow\rangle$ ，由于$\sigma^z$会反转$|\leftarrow\rangle$和$|\rightarrow\rangle$，因此如果打开$h_z$，会存在一个domain逐步反转这些态，这个反转需要进行N次，其振幅应该有$t\sim e^{-N/\xi}$形式。因此此时有
$$
H_{eff} = \begin{pmatrix}
  0&-t \\
  -t&0
\end{pmatrix}
$$
这使得$h_z\ne 0$让基态的简并被破缺，出现一个$2t$的间隙，只有当$N\to\infty$的时候才会回到简并。

## 2.2 两个理论的关系
 我们可以把TFIM和SPSC连接起来，这个所使用的是Jordan-Wigner变换。首先需要让自旋上和下的态分别对应到费米子数0和1
$$
|\uparrow\rangle\longleftrightarrow |0\rangle, \qquad |\downarrow\rangle\longleftrightarrow |1\rangle
$$
此时有对应$\sigma_j^z=(-1)^{a_j^\dagger a_j}$，或者可以写作
$$
\sigma^z_j = 1-2a_j^{\dagger}a_j
$$
回顾自旋升降算符$\sigma^\pm_j$，此时
$$
\sigma_j^+=\frac{\sigma_j^x+i\sigma_j^y}{2}, \qquad \sigma_j^-=\frac{\sigma_j^x-i\sigma_j^y}{2}
$$
此时可以自然的对应$\sigma^{+}_j=a_j$，以及$\sigma^{-}_j=a^{\dagger}_j$。可是这里有个巨大的问题是$\sigma^{\pm}_j$没有反对易关系，Jordan-Wigner的解决方法是在上面乘上一条由所有左侧自旋构成的string
$$
a_j = \left(\prod_{k=1}^{j-1}\sigma^z_k\right)\sigma_j^{+}\quad a_j^{\dagger} = \left(\prod_{k=1}^{j-1}\sigma^z_k\right)\sigma_j^{-}
$$
这个String常常记为$S_j$，即$a_j = S_j\sigma_j^{+}\quad a_j^{\dagger} = S_j\sigma_j^{-}$ 。现在我们可以把$\sigma^x_j = \sigma^+_j+\sigma^-_j = S_j(a_j+a_j^{\dagger})$代入TFIM的Hamiltonian，有
$$
\sigma^x_j\sigma^x_{j+1} = S_j(a_j+a_j^{\dagger}) S_j\sigma^z_j(a_j+a_j^{\dagger}) = (a_j+a_j^{\dagger}) \sigma^z_j(a_{j+1}+a_{j+1}^{\dagger})
$$
代入$\sigma^z_j = (1-2a^\dagger_ja_j)$，得到
$$
-(a_j-a^{\dagger}_j)(a_{j+1}+a^{\dagger}_{j+1})
$$
因此参数对应关系是
$$
w = J,\quad \mu = -2h_z
$$
两个对称性此时也显然是对应的


## 2.3 Majorana

我们现在定义$2N$个Majorana变量
$$
c_{2j-1} = a_j+a_j^{\dagger}\quad c_{2j}=-i(a_j-a_j^{\dagger})
$$
此时有代数关系$\{c_k,c_l\}=2\delta_{kl}$。此时Hamilitonian为
$$
H = \frac{i}{2}\left(  v\sum_{j=1}^Nc_{2j-1}c_{2j}+w\sum_{j=1}^{N-1}c_{2j}c_{2j+1}\right)
$$
对称性为
$$
P = \prod_{k=1}^{N}(-ic_{2k-1}c_{2k})
$$
在$v=0$的时候，只有后面一项有贡献，两侧的$c_1,c_{2N}$此时完全不会在Hamiltonian里面出现（这些和Hamiltonian对易的称为Majorana零模，此时会形成简并态）。若我们一般定义左右侧边界的Majorana算符为$b_l,b_r$，定义
$$
a=\frac{b_l+ib_r}{2},\qquad a^\dagger=\frac{b_l-ib_r}{2}
$$
在$v\ne0$并$|v|<w$，此时的effective Hamiltonian为
$$
H_{eff} = \frac{i}{2}\epsilon b_lb_r = \epsilon \left(a^{\dagger}a-\frac{1}{2}\right)
$$

对一般的Quadratic Majorana Hamiltonian，A反对称
$$
H(A) = \frac{i}{4}\sum_{j,k=1}^{2N}A_{jk}c_jc_k
$$
其中A一般可以写为如下形式
$$
A =  Q
\begin{pmatrix} 0&\epsilon_1&&\\ -\epsilon_1&0&&\\ &&0&\epsilon_2\\ &&-\epsilon_2&0\\ &&&&\ddots \end{pmatrix}Q^T\quad Q\in O(2N),\epsilon_k\ge 0
$$
对之前的模型来说，$\epsilon_{4j-3}$都是$v$，$\epsilon_{4j-1}$都是$w$。我们如果此时希望找到Majorana零模的话，即找到$[H,b]=0$的话，代入$b=\sum_lu_lc_l$，此时得到
$$
[H,b] = -i\sum_{j,k}u_jA_{jk}c_k
$$
因此此时的$uA=0$可以解得
$$
u = (1,0,\frac{v}{w},0,\left(\frac{v}{w}\right)^2,\cdots)
$$
即
$$
b_l\propto c_1+\frac vw c_3 +\left(\frac vw\right)^2c_5+\cdots
$$
不过，现实中可能存在相互作用，此时Hamiltonian不再是严格的二次型，此时边界态是否还能继续存在？作者给出了一个确切的构造，此时我们有费米子宇称算符
$$
P = \prod_{k=1}^{N}(-ic_{2k-1}c_{2k})
$$
进一步，构造
$$
X = c_{2j+1}\prod_{k=1}^{j}(-ic_{2k-1}c_{2k})
$$
我们可以证明
$$
[H,P]=[H,X]=\{X,P\}=0
$$
这使得$X|\psi_+\rangle$是一个宇称不同的基态。这确保了简并性。



# 3 The two-dimensional toric code

这个其实我觉得我很熟了，不熟的建议去玩小游戏https://liukelyu.github.io/toric-code-playground/，但为了完整性还是看看。
$$
H = -J_e \sum_s A_s - J_m \sum_p B_p
$$
$$
A_s = \prod_{j\in v}\sigma^x_j\quad B_p = \prod_{j\in\partial p}\sigma^z_{j}
$$
利用$\{\sigma^i,\sigma^j\}=2\delta^{ij}I$，可以证明任意两个$A_s,B_p$交换。我们可以定义plaquette flux，对每个边有$s_j=\pm1$
$$
w_p(s) = \prod_{j\in\partial p}s_j
$$
如果$w_p=-1$，则说$p$上有一个vortex。

现在考虑这个模型的基态，首先要求
$$
A_s|\Psi\rangle = |\Psi\rangle
$$
我们选取$\sigma^x$的本征态$|\pm\rangle$ ，同时我们定义$|+\rangle$是所谓的空，而$|-\rangle$是所谓的“有弦占据”。此时即要求每个顶点连接了偶数个弦。现在回顾$\sigma^z$会把$|\pm\rangle$反转，这意味着在进一步要求
$$
B_p |\Psi\rangle = |\Psi\rangle
$$
的话，我们必须要求在把每个闭合方格的基态都反转后仍然是基态的一部分，这要求所有弦回路都是闭合的，并且对所有闭合构型求和。
$$
|GS\rangle = \sum|all\ closed\ string\ loops\rangle
$$
如果我们用$\sigma^z$本征态表达的话，则是
$$
|GS\rangle = \mathcal{N}\prod_s(1+A_s)|all\ \uparrow\rangle
$$
我们可以从此理解为什么genus为$g$的Riemann面上有$2^{2g}$个独立基态，因为对应了cycle的数量。

我们现在可以理解一下激发态了，以弦图像（即$\sigma^x$本征态图像）演示，此时一个最简单的激发态是连接两个点之间的弦（及其在拓扑形变下的所有构型）
$$
W^{(e)}_l = \prod_{j\in l}\sigma^z_l
$$
这与所有$B_p$对易，也与除了端点之外的$A_s$对易，这个态的能量显然是$2J_e$。我们同时可以定义对偶的“m弦”
$$
W^{(m)}_{l^*} = \prod_{j\in l^*}\sigma^x_j
$$
这相当于在对偶格点里面的弦，给出一个$2J_m$的激发态。

# 4 Abelian anyons and quasiparticle statistics

## 4.1 fusion rules
众所周知，我们三次元里面只有费米子和玻色子两种统计，但是二次元他不一样，这里面有anyon。交换两个粒子本质上就是进行一个worldline braiding。我们可以先考虑一下刚才的toric code里面两种粒子，就是e charge和vortex m，我们只需要把e在m外面绕一圈就行。此时的态是
$$
|\xi\rangle\to\left(\prod_{j\in l}\sigma^z_j\right)|\xi\rangle = \left(\prod_{p\ in\ l}B_p\right)|\xi\rangle
$$
显然此时会变成$-|\xi\rangle$。这说明e和m交换两次会出一个负号。

e和e以及m和m显然是交换不变的。我们现在可以考虑一种e和m的复合粒子，画图可知这种粒子和自己交换两次也可以出一个负号。我们把这种粒子定义为$\epsilon$。

我们有一个fusion rules：
$$
e\times e=1\quad e\times m=\epsilon\\
m\times m=1\quad e\times \epsilon =m\\
\epsilon\times\epsilon = 1\quad m\times \epsilon = e
$$
## 4.2 Toric code in a field

我们现在考虑一种不是精确可解的模型，这将toric code放在一个磁场里面
$$
H = -J_e \sum_s A_s - J_m \sum_p B_p - \sum_j(h_x \sigma^x_j+h_z \sigma^x_j)
$$
我们显然也可以在这个模型中引入symmetry operators（不过，显然需要进行一些形式上的改变，因为你可以发现之前的$A_s$和$B_p$此时都和H不再对易了），但是如果要定义之前的loop operator的话需要对基态有一个完整的理解，但这个非常难算。我们在这里选择研究这个理论的gauge invariance来理解这个模型。

怎么做呢？我们现在给每个顶点也放一个自旋$\mu_v$，这相当于引入了新的一堆自由度，此时的算符应该为
$$
\sigma^z_{uv} = \tilde{\mu}^z_u\tilde{\sigma}^z_{uv}\tilde{\mu}^z_v\\
\sigma^x_{uv} = \tilde{\sigma}^x_{uv}\\
\mu^z_u = \tilde{\mu}^z_u\\
\mu^x_u=\tilde{\mu}^x_u\prod_{j\in star(u)}\tilde{\sigma}^x_j = \tilde{\mu}^x_u\tilde{A}_u
$$
为了消除我们新引入的这些自由度，我们引入一个规范
$$
G_u|\Psi\rangle := \tilde{\mu}^x_u\tilde{A}_u|\Psi\rangle = |\Psi\rangle
$$
这个是一个$\mathbb{Z}_2$ gauge constraint，这个作用在一个顶点上之后，显然把顶点上以及四条边上的z特征值全部反转。我们可以发现这样的反转不会改变$\sigma^z_{uv}$的值。同时$G^2_v = 1$，因此这是一个$\mathbb{Z}_2$对称性。

我们此时可以重写Hamiltonian为
$$
H={}-J_e\sum_v\widetilde\mu_v^x -J_m\sum_p\widetilde B_p\sum_{\langle u,v\rangle} \left( h_x\widetilde\sigma_{uv}^x +h_z\widetilde\mu_u^z \widetilde\sigma_{uv}^z \widetilde\mu_v^z \right)
$$

# 5 The honeycomb lattice model

## 5.1 Solving the honeycomb lattice model
我们现在考虑另一个可以完全求解的二维模型，即Honeycomb lattice model，这个模型定义为
$$
H = -J_x\sum_{x\ links}\sigma^x_j\sigma^x_k-J_y\sum_{y\ links}\sigma^y_j\sigma^y_k-J_z\sum_{z\ links}\sigma^z_j\sigma^z_k
$$
这个模型的自由度是honeycomb的vertice上的spin，但是这个模型的特殊之处是其相互作用是与方向有关的（link orientation dependent）。我们定义plaqutte operator
$$
W_p = \sigma^x_1\sigma^y_2\sigma^z_3\sigma^x_4\sigma^y_5\sigma^z_6
$$
显然可以验证$[H,W_p]=0$。我们现在引入若干Majorana算符$c,b^x,b^y,b^z$，并且定义
$$
\tilde{\sigma}^x = ib^xc\quad \tilde{\sigma}^y = ib^yc\quad
\tilde{\sigma}^z = ib^zc
$$
这个操作其实是增加了自由度（本来一个vertice上只有一个复费米子自由度，现在变成了俩），所以我们现在依旧要取一个规范来消除冗余。首先观察得到
$$
(\tilde{\sigma}^x)^2 = (ib^xc)^2=1
$$
与
$$
(\sigma^x)^2 = 1
$$
一致，但是如果我们检查代数关系的话
$$
\tilde{\sigma}^x\tilde{\sigma}^y = b^xb^y\ne i\tilde{\sigma}^z
$$
因此我们现在引入规范
$$
D = b^xb^yb^z c
$$
因为
$$
\tilde{\sigma}^x\tilde{\sigma}^y = iD\tilde{\sigma}^z
$$
所以我们要求所有物理态需要满足
$$
D|\Psi\rangle = |\Psi\rangle
$$
显然$D^2 =1$，因此这也是一个$\mathbb{Z}_2$规范。我们现在对所有节点都引入4个Majorana以及这个规范
$$
D_j|\Psi\rangle = |\Psi\rangle
$$
投影算符为
$$
\prod_j = \frac{1+D_j}{2}
$$
现在我们把刚才引入的新变量代回原来的理论，比如
$$
H^{(x)} = -J_x \sigma^x_j\sigma^x_k
$$
代入Majorana表示为
$$
H^{(x)} = -J_x(ib^x_jc_j)(ib^x_kc_k) = -J_xb^x_jb^x_kc_jc_k
$$
定义$\hat{u}_{jk} = ib^x_jb^x_k$，则
$$
H^{(x)} = iJ_x\hat{u}_{jk}c_jc_k
$$
或者更一般的写为$\hat{u}_{jk} = ib^{\alpha(j,k)}_jb^{\alpha(j,k)}_k$，其中$\alpha(j,k)$随$j,k$之间的属性确定为$x,y,z$。此时整个Hamiltonian可以写为
$$
H = \frac{i}{4}\sum_{j,k}\hat{A}_{j,k}c_jc_k\quad\hat{A}_{jk}=2J_{\alpha(j,k)}\hat{u}_{jk}
$$
我们可以发现
$$
[\hat{u}_{jk},H]=0\quad[\hat{u}_{jk},\hat{u}_{lm}]=0\quad \hat{u}^2_{jk}=1
$$
我们可以固定每条边上的$u_{jk}$，此时Hamiltonian变为一个简单的二次型，而$A$变为反对称矩阵，这就是第二章所提到的形式。之前的$D_j$此时可以理解为将每个顶点上的四个Majorana全部翻转
$$
b_j^\alpha\longrightarrow-b_j^\alpha, \qquad c_j\longrightarrow-c_j
$$
此时$\tilde{\sigma}$显然不变，但是$u_{jk}$会变号，这说明单独的$u_{jk}$并不是规范不变的物理可观测量，而真正有物理意义的是
$$
W_p = \prod_{\partial p}u_{jk}
$$
我们于是可以将$W_p$理解为六边形中的$\mathbb{Z}_2$通量。现在求解整个模型的方法已经明确，首先选定每个六边形的通量$w_p=\pm 1$，然后选择一组边变量$u_{jk}$满足
$$
w_p = \prod_{\partial p}u_{jk}
$$
在这组固定的$u_{jk}$下求解对应通量扇区的能量，得到相应的态之后使用投影算符回到真正的物理态空间。对于基态，Lieb 1994证明了所有$w_p=1$，即此时没有vortex。

## 5.2 Anyons in the honeycomb lattice model

我们首先可以把honeycomb里面的格点分为两类，暂且称为A类和B类，这种分类使得每一个边都连接一个A格点和一个B格点。如果我们只看A或者B格点，他们都分别组成一个边长是现在$\sqrt{3}$倍的三角晶格。

每个原胞现在有一个A和一个B格点，因此有两个Majorana
$$
c_{A,\mathbf{R}}\quad c_{B,\mathbf{R}}
$$
我们现在选择一种规范，对$u_{jk}$，定义从A指向B的时候值为$1$。现在整个Hamiltonian为
$$
H=i\sum_{\langle j,k\rangle_\alpha}^{j\in A,k\in B} J_\alpha c_jc_k
$$
我们记原胞平移矢量为$\mathbf{n}_1$和$\mathbf{n}_2$，此时对一个位于$\mathbf{R}$的A顶点，所连接的三个B顶点显然是$\mathbf{R}$，$\mathbf{R}+\mathbf{n}_1$，$\mathbf{R}+\mathbf{n}_2$，因此Hamiltonian为
$$
H=i\sum_{\mathbf R}\Big[ J_z c_{A,\mathbf R}c_{B,\mathbf R}+J_x c_{A,\mathbf R}c_{B,\mathbf R+\mathbf n_1}+J_y c_{A,\mathbf R}c_{B,\mathbf R+\mathbf n_2} \Big]
$$
我们进行Fourier变换
$$
c_{A,\mathbf R} =\frac1{\sqrt M}\sum_{\mathbf q} e^{i\mathbf q\cdot\mathbf R}c_{A,\mathbf q}\qquad c_{B,\mathbf R} =\frac1{\sqrt M}\sum_{\mathbf q} e^{i\mathbf q\cdot\mathbf R}c_{B,\mathbf q}
$$
以$c_{A,\mathbf{q}},c_{B,\mathbf{q}}$为基，我们得到一个矩阵
$$
iA(\mathbf q)= \begin{pmatrix} 0&if(\mathbf q)\\ -if^*(\mathbf q)&0 \end{pmatrix}
$$
其中
$$
f(\mathbf q)=2\left[ J_z+ J_xe^{i\mathbf q\cdot\mathbf n_1} +J_ye^{i\mathbf q\cdot\mathbf n_2} \right]
$$
此时能谱就是解$f(\mathbf q)=0$的方程，这即要求三个长度分别为$J_x,J_y,J_z$的向量形成一个三角形。显然，当三角不等式无法满足的时候，不存在零能费米子激发，此时就可以打开能隙。

我们现在研究这个模型里面的anyons。首先研究
$$
J_x=J_y=0\quad J_z>0
$$
的条件，此时理论变为若干独立dimer（即双自旋对），其中每一个dimer有四个状态，两个简并基态。我们可以把此时的两个状态理解为一个有效spin-$\frac{1}{2}$
$$
|\Uparrow\rangle :=|\uparrow\uparrow\rangle\quad|\Downarrow\rangle :=|\downarrow\downarrow\rangle
$$
现在我们打开很小的$J_x,J_y$，对x方向考虑简并微扰论，我们最终在四阶微扰处得到一个相互作用
$$
H_{\mathrm{eff}}^{(4)} =\mathrm{const} -\frac{J_x^2J_y^2}{16J_z^3}\sum_p Q_p
$$
此时$Q_p$是$\tau^y_1\tau^x_2\tau^y_3\tau^x_4$，其中$\tau$均为有效dimer自旋。因此在对两个方向均考虑微扰论并且旋转一下有效dimer自旋的情况下，我们可以得到Toric code Hamiltonian！此时我们发现在$J_z>>J_x,J_y$得到的是一个Abelian anyon理论。

我们可以从另一个方法得到non Abelian anyon，此时需要进入gapless phase。但是后两节我实在没看懂，之后看其他资料吧。