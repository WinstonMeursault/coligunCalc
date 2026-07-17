# Numerical Simulation of Single/Multi-Stage Synchronous Induction Coilgun

## References

1. **Ni et al.** (2015) — "Calculation of Temperature Rise of Armature in Multi-Stage Synchronous Inductive Coilgun Based on Current Filament Method", *High Power Laser and Particle Beams*, 27(9), 095001.
2. **Xiang H.** (2015) — *Principles and Technology of Electromagnetic Induction Coilgun*, Weapon Industry Press, Beijing.
3. **Li X., Xu Z., Li L.** (2015) — *Principles of Electromagnetic Coil Launch*, Weapon Industry Press, Beijing.
4. **Wu S.** (2003) — "Inductance Calculation of Air-Cored Cylindrical Coils", M.Eng. Thesis, Zhengzhou University.

Throughout this document, superscripts refer to the reference number in this list. For example, $F_x = \frac{dM}{dx} I_1 I_2$ [2].

---

## 1. Physical Principles

### 1.1 Working Principle

A synchronous induction coilgun converts electrical energy stored in capacitors into kinetic energy of a projectile (armature) via electromagnetic induction. The system consists of:

- **Pulse power supply**: capacitor banks charged to high voltage
- **Trigger switch**: high-voltage, high-current closing switch
- **Driving coil**: multi-turn helical coil (fixed stator)
- **Armature**: cylindrical conductor (typically aluminum or copper)
- **Projectile/payload**: attached to the armature

When the trigger switch closes, the capacitor discharges into the driving coil, producing a pulsed current $I_{\mathrm{d}}(t)$ and a transient magnetic field. By Faraday's law, this time-varying field induces eddy currents in the armature. According to Lenz's law, the induced current $I_{\mathrm{p}}$ flows in the opposite direction to $I_{\mathrm{d}}$. The interaction between the driving coil's magnetic field and the armature's induced current produces a Lorentz force that accelerates the armature along the barrel axis [1,2].

For multi-stage systems, additional driving coils are triggered sequentially as the armature passes each stage, continuously accelerating the projectile [1].

> **Figure 1.** Principle of coilgun (magnetic flux density in space and induced current in armature). Source: [1] Fig. 1.

### 1.2 Force Calculation by the Inductance Method

Two methods exist for computing the electromagnetic force on the armature: the Ampere force method and the inductance (energy) method. The inductance method, based on the principle of virtual work, provides both accuracy and physical insight [2,3].

For a system of two coaxial coils with self-inductances $L_{\mathrm{d}}$, $L_{\mathrm{p}}$ and mutual inductance $M$, the total stored magnetic energy is [2,3]:

$$
W_{\mathrm{m}} = \frac{1}{2} L_{\mathrm{d}} I_{\mathrm{d}}^2 + \frac{1}{2} L_{\mathrm{p}} I_{\mathrm{p}}^2 + M I_{\mathrm{d}} I_{\mathrm{p}} \tag{1.1}
$$

Since only $M$ varies with armature position $x$, the axial force is the derivative of the magnetic energy:

$$
\boxed{F_x = \frac{\mathrm{d}W_{\mathrm{m}}}{\mathrm{d}x} = \frac{\mathrm{d}M}{\mathrm{d}x} I_{\mathrm{d}} I_{\mathrm{p}}} \tag{1.2}
$$

The force is proportional to the mutual inductance gradient $\mathrm{d}M/\mathrm{d}x$ and the product of the two currents. The sign convention: when $\mathrm{d}M/\mathrm{d}x$ is negative and $I_{\mathrm{d}}, I_{\mathrm{p}}$ have opposite directions, the armature experiences a forward (positive $x$) thrust [2].

> **Figure 2.** Mutual inductance $M$ and inductance gradient $\mathrm{d}M/\mathrm{d}x$ vs. axial distance between coils. Source: [2] Fig. 2-2. The gradient peaks at the "inductance length" $l_{\mathrm{m}} \leqslant r_{\mathrm{d}}$ from the center plane. When the armature center passes the driving coil center plane ($\mathrm{d}M/\mathrm{d}x = 0$), the force reverses — this must be avoided by removing or reversing one of the currents before crossover [3].

---

## 2. Current Filament Method

### 2.1 Motivation

In a solid cylindrical armature, induced eddy currents are not uniformly distributed due to the skin effect. The skin depth is [2]:

$$
\delta = \sqrt{\frac{2}{\omega \sigma \mu}} \tag{2.1}
$$

where $\omega$ is the driving coil discharge current frequency, $\sigma$ is the armature conductivity, and $\mu$ is the armature permeability.

To capture this non-uniform current distribution, the current filament method (CFM) discretizes the solid armature into $m \times n$ concentric current-carrying ring elements — $m$ divisions axially and $n$ divisions radially [1,2].

> **Figure 3.** Subdivision of solid armature into $m$ axial slices $\times$ $n$ radial layers. Source: [2] Fig. 2-4.

### 2.2 Assumptions

1. Each current filament ring has a sufficiently small cross-section that the induced current is uniformly distributed across it.
2. Each filament behaves as an independent single-turn coil with its own resistance $R_{ij}$ and self-inductance $L_{ij}$.
3. The armature filaments move rigidly together; relative positions between filaments are fixed.
4. Mutual inductances between armature filaments $M_{aa,ij,kl}$ are constant (filaments are stationary relative to each other).
5. Mutual inductances between driving coils and armature filaments $M_{\mathrm{d},ij}$ vary with armature position $x$.

> **Figure 4.** Equivalent circuit model of the current filament method. Source: [1] Fig. 2 and [2] Fig. 2-5.

### 2.3 Filament Geometry

For a cylindrical armature with inner radius $r_{\mathrm{ai}}$, outer radius $r_{\mathrm{ae}}$, length $l_{\mathrm{a}}$, and center position $x$:

**Radial discretization** ($n$ layers, $j = 1, 2, \ldots, n$):

$$
\delta r = \frac{r_{\mathrm{ae}} - r_{\mathrm{ai}}}{n}
$$

Inner radius of filament $j$:

$$
r_{\mathrm{a},j}^{\mathrm{(in)}} = r_{\mathrm{ai}} + (j - 1) \cdot \delta r
$$

Outer radius of filament $j$:

$$
r_{\mathrm{a},j}^{\mathrm{(out)}} = r_{\mathrm{ai}} + j \cdot \delta r
$$

Mean radius of filament $j$:

$$
\bar{r}_{\mathrm{a},j} = r_{\mathrm{ai}} + \left(j - \frac{1}{2}\right) \delta r
$$

**Axial discretization** ($m$ slices, $i = 1, 2, \ldots, m$):

$$
\delta l = \frac{l_{\mathrm{a}}}{m}
$$

Axial position (absolute coordinate) of filament $(i,j)$:

$$
x_{\mathrm{a},ij} = x - \frac{l_{\mathrm{a}}}{2} + \left(i - \frac{1}{2}\right) \delta l
$$

The total number of filament loops is $N_{\mathrm{f}} = m \times n$.

---

## 3. Circuit Equations

### 3.1 Single-Stage, Single Filament Equivalent Current Ring Model

For initial understanding, treat the entire armature as a single equivalent current ring. The circuit equations are [2]:

**Driving coil loop:**

$$
L_{\mathrm{d}} \frac{\mathrm{d}I_{\mathrm{d}}}{\mathrm{d}t} - M \frac{\mathrm{d}I_{\mathrm{p}}}{\mathrm{d}t} - v_{\mathrm{p}} I_{\mathrm{p}} \frac{\mathrm{d}M}{\mathrm{d}x} + R_{\mathrm{d}} I_{\mathrm{d}} = U_{\mathrm{C}} \tag{3.1}
$$

**Armature loop:**

$$
L_{\mathrm{p}} \frac{\mathrm{d}I_{\mathrm{p}}}{\mathrm{d}t} - M \frac{\mathrm{d}I_{\mathrm{d}}}{\mathrm{d}t} - v_{\mathrm{p}} I_{\mathrm{d}} \frac{\mathrm{d}M}{\mathrm{d}x} + R_{\mathrm{p}} I_{\mathrm{p}} = 0 \tag{3.2}
$$

**Capacitor voltage** (decay equation):

$$
U_{\mathrm{C}}(t) = U_{\mathrm{C}}(0) - \frac{1}{C} \int_0^t I_{\mathrm{d}} \, \mathrm{d}\tau \tag{3.3}
$$

In matrix form:

$$
\boxed{[\dot{I}] = \bigl([L] - [M]\bigr)^{-1} \Bigl([U] + v_{\mathrm{p}}\Bigl[\frac{\mathrm{d}M}{\mathrm{d}x}\Bigr][I] - [R][I]\Bigr)} \tag{3.4}
$$

where:

$$
[L] = \begin{bmatrix} L_{\mathrm{d}} & 0 \\ 0 & L_{\mathrm{p}} \end{bmatrix}, \quad
[M] = \begin{bmatrix} 0 & M \\ M & 0 \end{bmatrix}, \quad
[R] = \begin{bmatrix} R_{\mathrm{d}} & 0 \\ 0 & R_{\mathrm{p}} \end{bmatrix}, \quad
[U] = \begin{bmatrix} U_{\mathrm{C}} \\ 0 \end{bmatrix}, \quad
[I] = \begin{bmatrix} I_{\mathrm{d}} \\ I_{\mathrm{p}} \end{bmatrix}
$$

Note: $\frac{\mathrm{d}}{\mathrm{d}t}(M I_{\mathrm{p}}) = M \frac{\mathrm{d}I_{\mathrm{p}}}{\mathrm{d}t} + I_{\mathrm{p}} \frac{\mathrm{d}M}{\mathrm{d}x} v_{\mathrm{p}}$, where $v_{\mathrm{p}} = \mathrm{d}x/\mathrm{d}t$ [2].

### 3.2 Single-Stage, Current Filament Model

With $N_{\mathrm{f}} = m \times n$ armature filaments, the system has $1 + N_{\mathrm{f}}$ circuit equations [2]:

**Driving coil** (equation index 0):

$$
R_{\mathrm{d}} I_{\mathrm{d}} + L_{\mathrm{d}} \frac{\mathrm{d}I_{\mathrm{d}}}{\mathrm{d}t} - \sum_{i=1}^{m} \sum_{j=1}^{n} \frac{\mathrm{d}}{\mathrm{d}t} \bigl(M_{\mathrm{d},ij} I_{ij}\bigr) = U_{\mathrm{C}} \tag{3.5}
$$

**Armature filament** $(i,j)$:

$$
R_{ij} I_{ij} + L_{ij} \frac{\mathrm{d}I_{ij}}{\mathrm{d}t} - \frac{\mathrm{d}}{\mathrm{d}t} \bigl(M_{\mathrm{d},ij} I_{\mathrm{d}}\bigr) + \sum_{\substack{k=1 \\ (k,l) \neq (i,j)}}^{m} \sum_{l=1}^{n} M_{aa,ij,kl} \frac{\mathrm{d}I_{kl}}{\mathrm{d}t} = 0 \tag{3.6}
$$

Expanding the time-derivative of position-dependent mutual inductance [1,2]:

$$
\frac{\mathrm{d}}{\mathrm{d}t} \bigl(M_{\mathrm{d},ij} I_{\mathrm{d}}\bigr) = M_{\mathrm{d},ij} \frac{\mathrm{d}I_{\mathrm{d}}}{\mathrm{d}t} + I_{\mathrm{d}} \frac{\mathrm{d}M_{\mathrm{d},ij}}{\mathrm{d}x} v_{\mathrm{p}} \tag{3.7}
$$

Collecting all equations into matrix form:

$$
\boxed{[\dot{I}] = \bigl([L] - [M_{\mathrm{I}}]\bigr)^{-1} \Bigl([U] + v_{\mathrm{p}}\Bigl[\frac{\mathrm{d}M_{\mathrm{I}}}{\mathrm{d}x}\Bigr][I] - [R][I] - [M][I]\Bigr)} \tag{3.8}
$$

**Matrix definitions** (all dimensions $(N_{\mathrm{f}}+1) \times (N_{\mathrm{f}}+1)$):

- $[R] = \mathrm{diag}(R_{\mathrm{d}}, R_{11}, \ldots, R_{1n}, \ldots, R_{mn})$ — diagonal resistance matrix
- $[L] = \mathrm{diag}(L_{\mathrm{d}}, L_{11}, \ldots, L_{1n}, \ldots, L_{mn})$ — diagonal self-inductance matrix
- $[U] = [U_{\mathrm{C}}, 0, 0, \ldots, 0]^{\mathsf{T}}$ — voltage source vector (capacitor voltage at index 0)
- $[I] = [I_{\mathrm{d}}, I_{11}, \ldots, I_{1n}, \ldots, I_{mn}]^{\mathsf{T}}$ — current state vector

**Mutual inductance matrices** (all with zeros on the diagonal) [2]:

- $[M_{\mathrm{I}}]$ — driving coil-to-armature filament mutual inductances. The first row and column contain $M_{\mathrm{d},ij}$ terms; all other entries are zero (armature filaments are stationary relative to each other, and $M_{\mathrm{I}}$ is reserved specifically for the coil-armature coupling):

$$
[M_{\mathrm{I}}] = \begin{bmatrix}
0 & M_{\mathrm{d},11} & \cdots & M_{\mathrm{d},1n} & \cdots & M_{\mathrm{d},mn} \\
M_{\mathrm{d},11} & & & & & \\
\vdots & & & \text{0} & & \\
M_{\mathrm{d},1n} & & & & & \\
\vdots & & & & & \\
M_{\mathrm{d},mn} & & & & &
\end{bmatrix}
$$

- $\bigl[\frac{\mathrm{d}M_{\mathrm{I}}}{\mathrm{d}x}\bigr]$ — gradient of $[M_{\mathrm{I}}]$ with respect to position $x$, same sparsity pattern, contains $\frac{\mathrm{d}M_{\mathrm{d},ij}}{\mathrm{d}x}$ terms. This term represents the motional EMF due to the relative motion between driving coil and armature.
- $[M]$ — armature filament-to-filament mutual inductances $M_{aa,ij,kl}$ (for $i \neq k$ or $j \neq l$; zero on diagonal and in the first row/column). These are constant since filaments are fixed relative to each other.

Key physical insight: $v_{\mathrm{p}}[\mathrm{d}M_{\mathrm{I}}/\mathrm{d}x][I]$ in Eq. (3.8) represents the motional EMF (back-EMF) generated by the relative motion between driving coil and armature. The term $[M][\dot{I}]$ does not appear explicitly in Eq. (3.8) because $[M]$ is on the RHS — the filament-to-filament mutual coupling contributes to the effective resistance seen by each filament through the $[M][I]$ term (when combined with the inductive coupling through diagonal elements of $[L]$).

### 3.3 Multi-Stage Coilgun Circuit Equations

For an $s$-stage system, each stage $i$ has its own capacitor $C_i$, initial voltage $U_{i0}$, driving coil $L_{\mathrm{d}i}$, and resistance $R_{\mathrm{d}i}$. Not all stages are active simultaneously — only stages that have been triggered have non-zero current [1,2].

**Driving coil $i$ (active):**

$$
R_{\mathrm{d}i} I_{\mathrm{d}i} + L_{\mathrm{d}i} \frac{\mathrm{d}I_{\mathrm{d}i}}{\mathrm{d}t} - \sum_{j=1}^{m} \frac{\mathrm{d}}{\mathrm{d}t} \bigl(M_{\mathrm{d}i,aj} I_{aj}\bigr) + \sum_{\substack{k=1 \\ k \neq i}}^{n_{\mathrm{active}}} M_{\mathrm{d}i,\mathrm{d}k} \frac{\mathrm{d}I_{\mathrm{d}k}}{\mathrm{d}t} = U_{i} \tag{3.9}
$$

**Driving coil $i$ (not yet triggered):** $I_{\mathrm{d}i} = 0$.

**Armature filament $j$:**

$$
R_{aj} I_{aj} + L_{aj} \frac{\mathrm{d}I_{aj}}{\mathrm{d}t} + \sum_{k=1}^{n_{\mathrm{active}}} \frac{\mathrm{d}}{\mathrm{d}t} \bigl(M_{\mathrm{d}k,aj} I_{\mathrm{d}k}\bigr) + \sum_{\substack{l=1 \\ l \neq j}}^{m} M_{aa,jl} \frac{\mathrm{d}I_{al}}{\mathrm{d}t} = 0 \tag{3.10}
$$

**Capacitor voltage for stage $i$:**

$$
U_i = U_{i0} - \frac{1}{C_i} \int_{t_i}^{t} I_{\mathrm{d}i} \, \mathrm{d}\tau \quad (t \geq t_i) \tag{3.11}
$$

where $t_i$ is the trigger time of stage $i$.

**Full matrix form** (after the convention of [1], Eq. 6):

$$
\begin{bmatrix} QU \\ 0 \end{bmatrix} =
\begin{bmatrix} QR_{\mathrm{c}} Q^{\mathsf{T}} & v Q \dot{M}_{\mathrm{ca}} \\
v \dot{M}_{\mathrm{ac}} Q^{\mathsf{T}} & R_{\mathrm{a}} \end{bmatrix}
\begin{bmatrix} Q I_{\mathrm{c}} \\ I_{\mathrm{a}} \end{bmatrix} +
\begin{bmatrix} Q L_{\mathrm{c}} Q^{\mathsf{T}} & Q M_{\mathrm{ca}} \\
M_{\mathrm{ac}} Q^{\mathsf{T}} & L_{\mathrm{a}} \end{bmatrix}
\begin{bmatrix} Q \dot{I}_{\mathrm{c}} \\ \dot{I}_{\mathrm{a}} \end{bmatrix} \tag{3.12}
$$

where:
- $R_{\mathrm{c}}, L_{\mathrm{c}}$: $n \times n$ diagonal matrices of driving coil resistances and self-inductances (off-diagonals are inter-coil mutual inductances $M_{\mathrm{cc},ki}$ for $L_{\mathrm{c}}$)
- $\dot{M}_{\mathrm{ca}}$: $n \times m$ matrix of mutual inductance gradients between driving coils and armature filaments
- $M_{\mathrm{ca}}$: $n \times m$ matrix of mutual inductances between driving coils and armature filaments
- $\dot{M}_{\mathrm{ac}} = \dot{M}_{\mathrm{ca}}^{\mathsf{T}}$, $M_{\mathrm{ac}} = M_{\mathrm{ca}}^{\mathsf{T}}$
- $R_{\mathrm{a}}$: $m \times m$ diagonal resistance matrix of armature filaments
- $L_{\mathrm{a}}$: $m \times m$ inductance matrix of armature filaments (diagonal: self-inductances; off-diagonal: inter-filament mutual inductances)
- $Q$: truncation matrix of dimension $i \times n$ (at time when $i$ stages are active):

$$
Q = \begin{bmatrix}
1 & 0 & \cdots & 0 \\
0 & 1 & \cdots & 0 \\
\vdots & \vdots & \ddots & \vdots \\
0 & 0 & \cdots & 1
\end{bmatrix}_{i \times n}, \quad i = 1, 2, \ldots, n \tag{3.13}
$$

The $Q$ matrix keeps the system dimension manageable as additional stages are triggered — only active stages contribute non-zero equations [1].

### 3.4 Inclusion of Crowbar (Freewheeling) Diode

In a practical coilgun, each driving coil's RLC discharge circuit is underdamped ($R_{\mathrm{d}} < 2\sqrt{L_{\mathrm{d}}/C}$) to maximize peak current. Without countermeasures, the capacitor voltage oscillates — after the first half-cycle it reverses polarity, causing the driving coil current to reverse direction. When $I_{\mathrm{d}}$ and $I_{\mathrm{p}}$ have the same sign, the product $\frac{\mathrm{d}M}{\mathrm{d}x} I_{\mathrm{d}} I_{\mathrm{p}}$ becomes negative, producing a **braking force** that decelerates the armature [2, §2.5].

The standard solution is a **crowbar (freewheeling) diode** $D$ connected in parallel with the capacitor, as shown in [2] Fig. 2-12.

#### 3.4.1 Circuit Model

The diode $D$ has two operating states:

**State 1 — Diode OFF** ($u_{AB} > 0$, capacitor voltage positive):

The diode is reverse-biased. The circuit is a standard series RLC. Driving coil current equals capacitor current, diode current is zero:

$$
I_{\mathrm{C}} = I_{\mathrm{d}}, \qquad I_{D} = 0 \tag{3.14}
$$

$$
U_{\mathrm{C}}(t) = U_0 - \frac{1}{C} \int_0^t I_{\mathrm{d}}(\tau) \, \mathrm{d}\tau \tag{3.15}
$$

The voltage at the capacitor terminals is $u_{AB}(t) = U_{\mathrm{C}}(t)$.

**State 2 — Diode ON** ($u_{AB}$ reaches zero from positive side):

The diode becomes forward-biased and conducts, short-circuiting the capacitor terminals. Hence $u_{AB} \approx 0$ (ideal diode: zero forward voltage). The capacitor is effectively removed from the circuit. The current commutates:

$$
I_{\mathrm{C}} = 0, \qquad I_{D} = I_{\mathrm{d}} \tag{3.16}
$$

The driving coil circuit becomes a simple **RL decay**:

$$
L_{\mathrm{d}} \frac{\mathrm{d}I_{\mathrm{d}}}{\mathrm{d}t} - \sum_{j} \frac{\mathrm{d}}{\mathrm{d}t}\bigl(M_{\mathrm{d},j} I_{aj}\bigr) + R_{\mathrm{d}} I_{\mathrm{d}} = 0 \tag{3.17}
$$

The current decays exponentially from its value at the moment of commutation:

$$
I_{\mathrm{d}}(t) \approx I_{\mathrm{d}}(t_{\mathrm{switch}}) \cdot \exp\!\left(-\frac{R_{\mathrm{d}}}{L_{\mathrm{d}}}(t - t_{\mathrm{switch}})\right) \tag{3.18}
$$

(the exact evolution accounts for the mutual coupling with armature filaments via the full matrix equation below).

> **Figure.** Circuit model of induction coilgun with crowbar diode. Source: [2] Fig. 2-12b.

#### 3.4.2 Matrix Formulation with Crowbar Diode

For each stage $i$, introduce a state flag $s_i \in \{0, 1\}$ where $s_i = 0$ means diode OFF (RLC mode) and $s_i = 1$ means diode ON (RL mode). The transition occurs when the capacitor voltage would become negative:

$$
\text{if } s_i = 0 \text{ and } U_i + \Delta U_i < 0 \text{, then } s_i \leftarrow 1, \; U_i \leftarrow 0 \tag{3.19}
$$

In matrix form, the system equation (3.8) is modified as follows. For stages with $s_i = 1$, the corresponding entry in the voltage vector $[U]$ is set to zero (the capacitor is shorted), but the driving coil current continues through the diode path:

$$
\boxed{[\dot{I}] = \bigl([L] - [M_{\mathrm{I}}]\bigr)^{-1} \Bigl([U_{\mathrm{eff}}] + v_{\mathrm{p}} [\dot{M}_{\mathrm{I}}] [I] - [R][I] - [M][\dot{I}]\Bigr)} \tag{3.20}
$$

where $[U_{\mathrm{eff}}]_0 = U_{\mathrm{C}}$ for the first stage diode-OFF, and $[U_{\mathrm{eff}}]_0 = 0$ when diode-ON. The capacitor voltage update in State 1 uses Eq. (3.15); in State 2 it remains clamped at zero.

**Physical effect**: With the diode, the driving coil current is always positive (unidirectional). The capacitor discharges completely during the first current pulse and does not recharge. This eliminates the braking force and improves both efficiency and capacitor lifetime [2, §2.5].

---

## 4. Inductance Calculations

### 4.1 Self-Inductance of a Hollow Cylindrical Coil

The self-inductance is derived from the magnetic vector potential $\mathbf{A}$ for an air-cored solenoid, assuming axisymmetric geometry and uniform current density [4]. This formula applies to **both** the multi-turn driving coil and the single-turn armature filaments — the only difference is the turns density $n_{\mathrm{c}}$.

For a coil with inner radius $R_1$, outer radius $R_2$, axial length $D$, and turns density $n_{\mathrm{c}}$ (turns per unit cross-sectional area):

$$
\boxed{L = 2\pi \mu_0 n_{\mathrm{c}}^2 R_1^5 \, T(q, p)} \tag{4.1}
$$

**For the multi-turn driving coil**: $n_{\mathrm{c}} = N_{\mathrm{d}} / [(R_2 - R_1) \cdot D]$, where $N_{\mathrm{d}}$ is the total number of turns.

**For each single-turn armature filament**: $n_{\mathrm{c}} = 1 / [(\delta r) \cdot (\delta l)]$, where $\delta r$ and $\delta l$ are the radial thickness and axial length of the filament (derived in §2.3). The formula remains valid because the current density is assumed uniform over the rectangular cross-section of the filament — an assumption justified by the small filament dimensions in the current filament method.

where $\mu_0 = 4\pi \times 10^{-7}$ H/m, and the shape parameters are:

$$
p = \frac{R_2}{R_1}, \quad q = \frac{D}{R_1} \tag{4.2}
$$

The dimensionless function $T(q,p)$ is given by the semi-infinite integral:

$$
\boxed{T(q, p) = \int_0^{\infty} U^2(1, p, x) \bigl(q x + e^{-q x} - 1\bigr) \, \mathrm{d}x} \tag{4.3}
$$

where $U(1, p, x)$ involves the first-order Bessel function $J_1$:

$$
U(R_1, R_2, \lambda) = \frac{1}{\lambda^3} \int_{\lambda R_1}^{\lambda R_2} t \, J_1(t) \, \mathrm{d}t \tag{4.4}
$$

After normalization with $x = \lambda R_1$:

$$
U(1, p, x) = \frac{1}{x^3} \int_x^{p x} t \, J_1(t) \, \mathrm{d}t \tag{4.5}
$$

**Derivation path** (from [4]):

1. Solve the boundary value problem for the magnetic vector potential $A_\phi$ for a single-turn circular loop.
2. Superpose to obtain the vector potential for the entire coil: $A_\phi(\rho, z)$.
3. Compute total flux linkage: $\Psi = \frac{1}{I} \int_{V'} \mathbf{A} \cdot \mathbf{J}_{\mathrm{c}} \, \mathrm{d}V'$.
4. Inductance: $L = \Psi / I$.

**Practical computation of $T(q,p)$** [4]:

- **Numerical integration**: Use Gauss-Laguerre quadrature with $\geq 15$ nodes to evaluate the infinite integral in Eq. (4.3).
- **Lookup table + bilinear interpolation**: For the parameter range $0.05 \leq q \leq 4.0$ and $1.05 \leq p \leq 4.0$, precomputed $T(q,p)$ tables are available. For values between table entries, use bilinear interpolation:

$$
\begin{aligned}
T(q, p) &= T(q_0, p_0) + \frac{T(q_1, p_0) - T(q_0, p_0)}{q_1 - q_0} (q - q_0) \\
&\quad + \frac{T(q_0, p_1) - T(q_0, p_0)}{p_1 - p_0} (p - p_0)
\end{aligned} \tag{4.6}
$$

where $q_0 < q < q_1$, $p_0 < p < p_1$.

> **Table 1.** Excerpt of the $T(q,p)$ function table (from [4] Table 4.1). Full table covers $q \in [0.05, 4.0]$ in steps of $0.05$, $p \in [1.05, 4.0]$ in steps of $0.05$. Maximum absolute error: $5 \times 10^{-5}$.

| $q \;/\; p$ | 1.50 | 1.55 | 1.60 | 1.80 | 2.00 | 2.50 | 3.00 |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| 0.10 | 0.0012 | 0.0014 | 0.0018 | 0.0032 | 0.0050 | 0.0122 | 0.0253 |
| 0.50 | 0.0202 | 0.0248 | 0.0303 | 0.0590 | 0.0956 | 0.2424 | 0.4920 |
| 1.00 | 0.0723 | 0.0884 | 0.1075 | 0.2133 | 0.3543 | 0.8894 | 1.7135 |
| 2.00 | 0.2133 | 0.2614 | 0.3119 | 0.5839 | 0.9359 | 2.1480 | 3.8932 |
| 3.00 | 0.3675 | 0.4540 | 0.5560 | 1.0715 | 1.8029 | 4.8452 | 10.2341 |
| 4.00 | 0.5246 | 0.6416 | 0.7814 | 1.5373 | 2.6963 | 8.2170 | 20.6373 |

### 4.2 Mutual Inductance Between Two Coaxial Circular Loops

For two infinitesimally thin circular loops of radii $a$ and $b$, coaxially separated by axial distance $h$, the mutual inductance has a closed-form expression using complete elliptic integrals [4]:

First define the elliptic modulus $k$:

$$
\boxed{k = \sqrt{\frac{4ab}{(a+b)^2 + h^2}}} \tag{4.7}
$$

Then the mutual inductance is:

$$
\boxed{M_{\text{loop}}(a, b, h) = \mu_0 \sqrt{ab} \left[ \left(\frac{2}{k} - k\right) K(k) - \frac{2}{k} E(k) \right]} \tag{4.8}
$$

where $K(k)$ and $E(k)$ are the complete elliptic integrals of the first and second kind, respectively:

$$
K(k) = \int_0^{\pi/2} \frac{\mathrm{d}\alpha}{\sqrt{1 - k^2 \sin^2 \alpha}}, \qquad
E(k) = \int_0^{\pi/2} \sqrt{1 - k^2 \sin^2 \alpha} \, \mathrm{d}\alpha \tag{4.9}
$$

**Derivation**: Start from Neumann's formula $M = \frac{\mu_0}{4\pi} \oint_{l_2} \oint_{l_1} \frac{\mathrm{d}\mathbf{l}_1 \cdot \mathrm{d}\mathbf{l}_2}{r}$, substitute the parametric forms for two coaxial circular loops, reduce via trigonometric identities, and integrate to the elliptic form [4].

> **Figure 5.** Coaxial circular loop geometry for mutual inductance calculation. Source: [4] Fig. 4.4.

### 4.3 Mutual Inductance Between Two Hollow Cylindrical Coils (Coil-to-Coil)

A practical coil with finite cross-section is treated as a continuous distribution of infinitesimal circular filament loops. The coil-to-coil mutual inductance is obtained by integrating the filament-level formula (Eq. 4.8) over both coil cross-sections [4]:

$$
M = \frac{N_1 N_2}{S_1 S_2} \iint_{S_1} \iint_{S_2} M_{\text{loop}}(r_1, r_2, |z_1 - z_2|) \, \mathrm{d}S_1 \, \mathrm{d}S_2 \tag{4.10}
$$

where $N_1$, $N_2$ are the turns of each coil, and $S_1$, $S_2$ are their cross-sectional areas. Parameterizing each coil's geometry explicitly:

$$
\boxed{M = \frac{N_1 N_2}{S_1 S_2} \int_{r_{11}}^{r_{12}} \mathrm{d}r_1 \int_{z_{11}}^{z_{12}} \mathrm{d}z_1 \int_{r_{21}}^{r_{22}} \mathrm{d}r_2 \int_{z_{21}}^{z_{22}} M_{\text{loop}}(r_1, r_2, |z_1 - z_2|) \, \mathrm{d}z_2} \tag{4.11}
$$

**Normalized coordinates for numerical integration** [4]:

Apply the coordinate transformation to the standard hypercube $[-1, 1]^4$:

$$
\begin{aligned}
r_k &= \frac{r_{k2} + r_{k1}}{2} + \frac{r_{k2} - r_{k1}}{2} r_k', \quad r_k' \in [-1, 1], \; k = 1, 2 \\
z_k &= \frac{z_{k2} + z_{k1}}{2} + \frac{z_{k2} - z_{k1}}{2} z_k', \quad z_k' \in [-1, 1], \; k = 1, 2
\end{aligned}
$$

The Jacobian determinant is:

$$
|J| = \frac{(r_{12} - r_{11})(z_{12} - z_{11})(r_{22} - r_{21})(z_{22} - z_{21})}{16} = \frac{S_1 S_2}{16} \tag{4.12}
$$

Substituting into Eq. (4.10) yields the normalized form:

$$
\boxed{M = \frac{N_1 N_2}{16} \int_{-1}^{1} \int_{-1}^{1} \int_{-1}^{1} \int_{-1}^{1} M_{\text{loop}}'(r_1', z_1', r_2', z_2') \, \mathrm{d}r_1' \, \mathrm{d}z_1' \, \mathrm{d}r_2' \, \mathrm{d}z_2'} \tag{4.13}
$$

where $M_{\text{loop}}'$ is $M_{\text{loop}}$ evaluated at the mapped coordinates. This 4D integral is the computationally dominant step in the simulation; it is evaluated using `scipy.integrate.nquad` or Gauss-Legendre quadrature in each dimension [4].

**Rectangular quadrature approximation** (for engineering use) [4]:

For $n$ quadrature nodes per dimension:

$$
M = N_1 N_2 \sum_{i=1}^{n} \sum_{j=1}^{n} \omega_i \omega_j M_{ij} \tag{4.14}
$$

> **Table 2.** Common quadrature nodes and coefficients for rectangular domain integration [4, Table 4.2].

| $n$ | Nodes $(x_k, y_k)$ | Weight $\omega_k$ | Error |
|:---:|:---|:---:|:---:|
| 4 | $(\pm\sqrt{1/3},\; \pm\sqrt{1/3})$ | $1/4$ | $O(h^6)$ |
| 9 | $(0,0),\ (\pm1,\pm1),\ (\pm1,0),\ (0,\pm1)$ | $4/9,\ 1/36,\ 1/9,\ 1/9$ | $O(h^6)$ |
| 9 | $(0,0),\ \bigl(\pm\sqrt{3/5},\ \pm\sqrt{3/5}\bigr),\ \bigl(\pm\sqrt{3/5},\ 0\bigr),\ \bigl(0,\ \pm\sqrt{3/5}\bigr)$ | $16/81,\ 25/324,\ 10/81,\ 10/81$ | $O(h^8)$ |

### 4.4 Mutual Inductance Gradient $\mathrm{d}M/\mathrm{d}x$

The force calculation (Eq. 1.2) requires the spatial derivative of mutual inductance with respect to axial displacement. **Reference [2] §3.1.3 derives both M and dM/dz for two coaxial circular loops** (Eqs. 3-46 and 3-47).

The mutual inductance $M(R_{\mathrm{d}}, R_{\mathrm{a}}, d)$ between two coaxial circular loops is given by Eq. (4.8) (replacing $a, b, h$ with $R_{\mathrm{d}}, R_{\mathrm{a}}, d$). The mutual inductance gradient is obtained by differentiating with respect to the axial separation $d$:

$$
\boxed{\frac{\mathrm{d}M}{\mathrm{d}z}(R_{\mathrm{d}}, R_{\mathrm{a}}, d) = \frac{\mu_0 \, k \, d}{4\,(1 - k^2)\,\sqrt{R_{\mathrm{d}} R_{\mathrm{a}}}} \Bigl[ 2\,(1 - k^2)\,K(k) - (2 - k^2)\,E(k) \Bigr]} \tag{4.15}
$$

where $k$ is the elliptic modulus defined in Eq. (4.7). $K(k)$ and $E(k)$ are complete elliptic integrals of the first and second kind with **modulus** $k$ [2, Eq. 3-47].

**Convention note**: Many computational libraries (including `scipy.special.ellipk` and `ellipe`) accept the **parameter** $m = k^2$ rather than the modulus $k$. When using such libraries, substitute $k = \sqrt{m}$ into Eq. (4.15):

$$
\frac{\mathrm{d}M}{\mathrm{d}z} = \frac{\mu_0 \sqrt{m} \, d}{4\,(1 - m)\,\sqrt{R_{\mathrm{d}} R_{\mathrm{a}}}} \Bigl[ 2\,(1 - m)\,K(m) - (2 - m)\,E(m) \Bigr] \tag{4.15b}
$$

**Physical sign**: $M$ is even in $d$, so $\mathrm{d}M/\mathrm{d}d$ is odd: $\mathrm{d}M/\mathrm{d}d < 0$ for $d > 0$ and $\mathrm{d}M/\mathrm{d}d > 0$ for $d < 0$. The total force sums over all driving coil–armature filament pairs (Eq. 5.1).

**Coil-to-coil mutual inductance gradient** (4D integration, same scheme as Eq. 4.13):

Substitute $\mathrm{d}M_{\text{loop}}/\mathrm{d}z$ (Eq. 4.15) for $M_{\text{loop}}$ in the 4D integration [2, §3.1.3]:

$$
\boxed{\frac{\mathrm{d}M}{\mathrm{d}x} = \frac{N_1 N_2}{16} \int_{-1}^{1} \int_{-1}^{1} \int_{-1}^{1} \int_{-1}^{1} \frac{\mathrm{d}M_{\text{loop}}}{\mathrm{d}z}\bigl(r_1', z_1', r_2', z_2'\bigr) \, \mathrm{d}r_1' \, \mathrm{d}z_1' \, \mathrm{d}r_2' \, \mathrm{d}z_2'} \tag{4.16}
$$

**Elliptic integral computation** [2, §3.1.3]: The arithmetic-geometric mean (AGM) method is recommended. Starting from $m_0 = 1$, $n_0 = \sqrt{1-k^2}$, iterate:

$$
m_i = \frac{m_{i-1} + n_{i-1}}{2}, \qquad n_i = \sqrt{m_{i-1} n_{i-1}}
$$

until convergence ($|m_N - n_N| / m_N < \varepsilon$). Then $K(k) \approx \pi / (2 m_N)$ and $E(k) \approx K(k) \cdot m_N$. Modern libraries (`scipy.special`) provide optimized implementations of $K(k)$, $E(k)$.

#### 4.4.1 Computational Symmetries and Optimization

The following symmetries can reduce the number of 4D integral evaluations by up to 50%:

**Reciprocity of $M$**: From Eq. (4.7), $k$ depends symmetrically on $a, b$ and on $d^2$ (even function of $d$). Hence:

$$
\boxed{M_{\text{loop}}(a, b, d) = M_{\text{loop}}(b, a, d) = M_{\text{loop}}(a, b, -d)} \tag{4.17}
$$

**Antisymmetry of $\mathrm{d}M/\mathrm{d}z$**: From Eq. (4.15), the numerator contains the factor $d$, making it odd in the signed separation:

$$
\boxed{\frac{\mathrm{d}M_{\text{loop}}}{\mathrm{d}z}(a, b, d) = -\frac{\mathrm{d}M_{\text{loop}}}{\mathrm{d}z}(a, b, -d)} \tag{4.18}
$$

Furthermore, swapping $a$ and $b$ whilst negating $d$:

$$
\frac{\mathrm{d}M_{\text{loop}}}{\mathrm{d}z}(a, b, d) = -\frac{\mathrm{d}M_{\text{loop}}}{\mathrm{d}z}(b, a, -d) \tag{4.19}
$$

In the 4D integration for coil-to-coil quantities, the driving coil is fixed and only armature filaments move. When building $[M_{\mathrm{I}}]$ and $[\dot{M}_{\mathrm{I}}]$ matrices at each time step, the integration is over the driving coil cross-section and each armature filament cross-section. The symmetry $M(a,b,d) = M(b,a,d)$ is always applicable between any two filaments regardless of which coil they belong to. For $[M]$ (armature inter-filament), the integration result is symmetric and computed once. For $[M_{\mathrm{I}}]$, the full 4D integral must be evaluated separately for each driving coil–armature filament pair at each time step because the axial separation $d$ changes with armature position.

**Edge cases for elliptic integrals**:

- **$k \to 0$** (filaments far apart, $d^2 \gg 4ab$): $K(k) \to \pi/2$, $E(k) \to \pi/2$, $M_{\text{loop}} \to 0$, $\mathrm{d}M_{\text{loop}}/\mathrm{d}z \to 0$. The dipole approximation $M \approx \mu_0 \pi a^2 b^2 / [2(a^2+b^2+d^2)^{3/2}]$ provides a numerically stable alternative.
- **$k \to 1$** (filaments very close, $d \to 0$ and $a \approx b$): $K(k)$ diverges logarithmically. In practice, clamp $k \leq 1 - \varepsilon$ with $\varepsilon \approx 10^{-12}$ when filaments approach zero separation. For filaments that coincide exactly ($d=0$, $a=b$), the mutual inductance formula is singular and the self-inductance formula (§4.1) must be used instead.

#### 4.4.2 AGM Algorithm for Elliptic Integrals

For reference implementations that compute $K(k)$ and $E(k)$ from scratch, the arithmetic-geometric mean (AGM) method converges quadratically in $\sim 5-8$ iterations to machine precision [2, §3.1.3].

Initialize: $m_0 = 1$, $n_0 = \sqrt{1-k^2}$.

Iterate for $i = 1, 2, \dots$:

$$
m_i = \frac{m_{i-1} + n_{i-1}}{2}, \qquad n_i = \sqrt{m_{i-1} \cdot n_{i-1}} \tag{4.20}
$$

Stop when $|m_i - n_i| / m_i < \varepsilon$ (e.g., $\varepsilon = 10^{-14}$).

Then the complete elliptic integral of the first kind is:

$$
K(k) = \frac{\pi}{2 \cdot m_i} \tag{4.21}
$$

For the second kind, the iterative scheme additionally tracks auxiliary sequences. Starting from $R_{a,0} = 1$, $R_{d,0} = 1-k^2$, $t_0 = 1/\sqrt{1-k^2}$:

$$
\begin{aligned}
t_i &= \frac{n_i}{4 m_i}\left(t_{i-1} + 2 + \frac{1}{t_{i-1}}\right) \\
R_{a,i} &= \frac{R_{a,i-1} + R_{d,i-1}}{2} \\
R_{d,i} &= \frac{R_{a,i-1} + R_{d,i-1} \cdot t_{i-1}}{1 + t_{i-1}}
\end{aligned} \tag{4.22}
$$

Convergence criteria are combined: $|R_{a,N} - R_{d,N}|/R_{a,N} < \varepsilon$, $|m_N - n_N|/m_N < \varepsilon$, $|t_N - 1| < \varepsilon$. Then:

$$
E(k) = \frac{\pi R_{a,N}}{2 m_N} \cdot K(k) \quad\text{(normalized)} \tag{4.23}
$$

In practice, `scipy.special.ellipk(m)` and `scipy.special.ellipe(m)` (with parameter $m = k^2$) should be used — they provide highly optimized, vectorized implementations.

---

## 5. Force and Motion Equations

### 5.1 Electromagnetic Force

From the principle of virtual work [1,2], the total axial force on the armature is the sum over all driving coil–armature filament pairs:

$$
\boxed{F_z = \sum_{i=1}^{n_{\mathrm{active}}} \sum_{j=1}^{m \times n} I_{\mathrm{d}i} \, I_{aj} \, \frac{\mathrm{d}M_{\mathrm{d}i, aj}}{\mathrm{d}z}} \tag{5.1}
$$

where $n_{\mathrm{active}}$ is the number of currently conducting driving coils, $m \times n$ is the total number of armature filaments, and $\mathrm{d}M_{\mathrm{d}i, aj}/\mathrm{d}z$ is the mutual inductance gradient between driving coil $i$ and armature filament $j$.

### 5.2 Equations of Motion

From Newton's second law [1]:

$$
\boxed{m_{\mathrm{a}} \frac{\mathrm{d}v}{\mathrm{d}t} = F_z = \sum_{i=1}^{n_{\mathrm{active}}} \sum_{j=1}^{m \times n} I_{\mathrm{d}i} \, I_{aj} \, \frac{\mathrm{d}M_{\mathrm{d}i, aj}}{\mathrm{d}z}} \tag{5.2}
$$

$$
\boxed{\frac{\mathrm{d}z}{\mathrm{d}t} = v} \tag{5.3}
$$

where $m_{\mathrm{a}}$ is the total accelerated mass (armature + payload).

### 5.3 Numerical Integration — Forward Difference Method

The system of circuit equations (Eq. 3.8) and motion equations (Eqs. 5.2–5.3) forms a coupled ODE system. The forward-difference (explicit Euler) method is employed for time-stepping [2]:

**Step 1 — Initialize** ($t = t_0$):
- $[I]_0 = 0$, $v_0 = 0$ (or given initial velocity)
- $[U]_0 = [U_{\mathrm{C}}(0), 0, \ldots, 0]^{\mathsf{T}}$
- Compute initial current derivatives:

$$
[\dot{I}]_0 = \bigl([L]_0 - [M_{\mathrm{I}}]_0\bigr)^{-1} \, [U]_0 \tag{5.4}
$$

**Step 2 — Time step $n$** ($t_n = t_{n-1} + \Delta t$):

Update capacitor voltage (**with crowbar diode check**, see §3.4):
$$
\begin{aligned}
&\text{If diode OFF } (s_i = 0): \quad U_{\mathrm{C},n}^{(i)} = U_{\mathrm{C},n-1}^{(i)} - \frac{\Delta t}{C_i} I_{\mathrm{d},n-1}^{(i)} \\
&\text{If } U_{\mathrm{C},n}^{(i)} < 0 \text{ (diode turns ON)}: \quad s_i \leftarrow 1, \; U_{\mathrm{C},n}^{(i)} \leftarrow 0 \\
&\text{If diode ON } (s_i = 1): \quad U_{\mathrm{C},n}^{(i)} = 0
\end{aligned} \tag{5.5}
$$
Set $[U_{\mathrm{eff}}]_n$: entries for diode-OFF stages use the capacitor voltage from (5.5); diode-ON stages use $0$.

Compute current derivatives with motional EMF [2]:
$$
\boxed{[\dot{I}]_n = \bigl([L]_n - [M_{\mathrm{I}}]_n\bigr)^{-1} \Bigl([U_{\mathrm{eff}}]_n + v_{n-1} [\dot{M}_{\mathrm{I}}]_n [I]_{n-1} - [R] [I]_{n-1} - [M] [I]_{n-1}\Bigr)} \tag{5.6}
$$

Update currents:
$$
[I]_n = [I]_{n-1} + \Delta t \, [\dot{I}]_{n-1} \tag{5.7}
$$

Compute force:
$$
F_n = \sum_{i=1}^{n_{\mathrm{active}}} \sum_{j=1}^{m \times n} I_{\mathrm{d}i,n} \, I_{aj,n} \, \left(\frac{\mathrm{d}M_{\mathrm{d}i, aj}}{\mathrm{d}z}\right)_n \tag{5.8}
$$

Compute acceleration, velocity, position:
$$
\begin{aligned}
a_n &= F_n / m_{\mathrm{a}} \\
v_n &= v_{n-1} + a_{n-1} \Delta t \\
z_n &= z_{n-1} + v_{n-1} \Delta t
\end{aligned} \tag{5.9}
$$

**Step 3 — Update geometry**: Recompute $[M_{\mathrm{I}}]_n$, $[\dot{M}_{\mathrm{I}}]_n$, and $[M]$ at the new position $z_n$.

**Step 4 — Trigger check**: If the armature reaches the trigger position of the next driving coil, activate it ($I_{\mathrm{d},n+1}$ begins conducting).

**Step 5 — Repeat** until $I_{\mathrm{d}}$ decays below a threshold or the armature exits the last coil.

### 5.4 Trigger Position

The optimal trigger position is where the armature center is at a distance such that the mutual inductance gradient is near maximum when the driving coil current reaches its peak [2,3]. Typically this is when the armature tail is near the driving coil's rear edge. The exact optimal position depends on the specific geometry and can be determined by parametric sweep simulation [3].

---

## 6. Temperature Rise Effects

### 6.1 Physical Model

Since the launch process is extremely brief (milliseconds), the armature heating is treated as **adiabatic** — heat conduction and radiation are neglected [1].

The ohmic heating power in each armature filament is $R I^2$:

$$
R I^2 = c_{\mathrm{p}}(T) \, m_{\mathrm{p}} \, \frac{\partial T}{\partial t} \tag{6.1}
$$

where $c_{\mathrm{p}}(T)$ is the temperature-dependent specific heat capacity (J/(kg·K)), $m_{\mathrm{p}}$ is the mass of the current filament (kg), and $T$ is absolute temperature (K).

### 6.2 Temperature-Dependent Specific Heat Capacity

From [1] Eq. 13, the specific heat capacities for aluminum and copper from room temperature to melting point follow exponential relationships:

**Aluminum** ($293 \,\mathrm{K} \leqslant T \leqslant 933 \,\mathrm{K}$, melting point):

$$
\boxed{c_{\mathrm{Al}}(T) = 0.819 \, \exp\bigl(3.07 \times 10^{-4} \, T\bigr) \quad \text{(kJ/(kg·K))}} \tag{6.2}
$$

**Copper** ($293 \,\mathrm{K} \leqslant T \leqslant 1356 \,\mathrm{K}$, melting point):

$$
\boxed{c_{\mathrm{Cu}}(T) = 0.333 \, \exp\bigl(3.917 \times 10^{-4} \, T\bigr) \quad \text{(kJ/(kg·K))}} \tag{6.3}
$$

> **Unit note**: Eqs. (6.2)–(6.3) give $c_{\mathrm{p}}$ in **kJ/(kg·K)**. When used in the temperature update Eq. (6.10) with SI quantities ($I^2 R \Delta t$ in joules, $m_{\mathrm{p}}$ in kg), multiply the specific heat by $1000$: $c_{\mathrm{p}}^{\mathrm{(SI)}} = 1000 \cdot c_{\mathrm{p}}$.

**Material densities** (for computing filament mass $m_{\mathrm{p}} = \rho_{\mathrm{m}} \cdot \pi (r_{\mathrm{out}}^2 - r_{\mathrm{in}}^2) \cdot \delta l$):
- Aluminum: $\rho_{\mathrm{Al}} = 2700 \; \mathrm{kg/m^3}$
- Copper: $\rho_{\mathrm{Cu}} = 8960 \; \mathrm{kg/m^3}$

### 6.3 Temperature-Dependent Resistivity and Filament Resistance

#### 6.3.1 Driving Coil Resistance

For a driving coil wound from flat copper wire of cross-sectional area $S_{\mathrm{wire}}$, with inner radius $r_{\mathrm{di}}$, outer radius $r_{\mathrm{de}}$, axial length $l_{\mathrm{d}}$, and fill factor $k$ (the ratio of total copper cross-section to the coil's rectangular cross-section, $0 < k < 1$), the total copper volume inside the coil is [2, §3.1.1]:

$$
V_{\mathrm{cu}} = k \cdot V_{\mathrm{coil}} = k \cdot \pi (r_{\mathrm{de}}^2 - r_{\mathrm{di}}^2) \, l_{\mathrm{d}} \tag{6.4}
$$

The total wire length is $l_{\mathrm{wire}} = V_{\mathrm{cu}} / S_{\mathrm{wire}}$, hence the DC resistance:

$$
\boxed{R_{\mathrm{d}} = \rho \cdot \frac{l_{\mathrm{wire}}}{S_{\mathrm{wire}}} = \rho \cdot \frac{k \pi (r_{\mathrm{de}}^2 - r_{\mathrm{di}}^2) \, l_{\mathrm{d}}}{S_{\mathrm{wire}}^2}} \tag{6.5}
$$

where $\rho$ is the wire resistivity. The turns density is $n_{\mathrm{c}} = N_{\mathrm{d}} / [(r_{\mathrm{de}} - r_{\mathrm{di}}) \cdot l_{\mathrm{d}}]$, related to $k$ through $N_{\mathrm{d}} S_{\mathrm{wire}} = k \cdot (r_{\mathrm{de}} - r_{\mathrm{di}}) \cdot l_{\mathrm{d}}$.

#### 6.3.2 Armature Filament Resistance

For an armature with inner radius $r_{\mathrm{ai}}$, outer radius $r_{\mathrm{ae}}$, length $l_{\mathrm{a}}$, discretized into $m$ axial slices and $n$ radial layers, each filament is a thin-walled circular ring [2, §3.1.1].

Define the radial increment $\delta r = (r_{\mathrm{ae}} - r_{\mathrm{ai}})/n$ and axial increment $\delta l = l_{\mathrm{a}} / m$.

For filament $(i,j)$ ($i = 1,\dots,m$ axial; $j = 1,\dots,n$ radial):

- **Cross-sectional area**: $A_{ij} = \delta l \cdot \delta r = \dfrac{l_{\mathrm{a}}}{m} \cdot \dfrac{r_{\mathrm{ae}} - r_{\mathrm{ai}}}{n}$
- **Outer boundary of radial layer $j$**: $r_j = r_{\mathrm{ai}} + j \cdot \delta r$
- **Inner boundary of radial layer $j$**: $r_{j-1} = r_{\mathrm{ai}} + (j-1) \cdot \delta r$
- **Mean radius of filament $(i,j)$**:

$$
\bar{r}_{ij} = \frac{r_j + r_{j-1}}{2} = r_{\mathrm{ai}} + \frac{r_{\mathrm{ae}} - r_{\mathrm{ai}}}{n} \cdot j - \frac{r_{\mathrm{ae}} - r_{\mathrm{ai}}}{2n} \tag{6.6}
$$

- **Current path length**: $l_{ij} = 2\pi \bar{r}_{ij}$ (circumference of the mean circular path)

Applying $R = \rho \cdot l / A$:

$$
\boxed{R_{ij} = 2\pi\rho \left[ \frac{m \cdot j}{l_{\mathrm{a}}} - \frac{m}{2l_{\mathrm{a}}} + \frac{m \cdot n \cdot r_{\mathrm{ai}}}{l_{\mathrm{a}}(r_{\mathrm{ae}} - r_{\mathrm{ai}})} \right]} \tag{6.7}
$$

This is the explicit formula from [2, Eq. 3-11]. Note that $R_{ij}$ depends on $j$ (radial index) but not on $i$ (axial index) — all filaments at the same radial layer have the same mean radius and therefore the same resistance.

#### 6.3.3 Temperature Dependence of Resistivity

#### 6.3.3 Temperature Dependence of Resistivity

Neglecting magnetoresistance effects, the primary factor affecting resistivity is temperature. The linear relationship is [1]:

$$
\boxed{\rho(T) = \rho_0 \, \bigl(1 + \beta \, \Delta T\bigr)} \tag{6.8}
$$

where $\beta$ is the temperature coefficient of resistivity:
- Copper: $\beta_{\mathrm{Cu}} = 4.1 \times 10^{-3} \, \mathrm{K}^{-1}$
- Aluminum: $\beta_{\mathrm{Al}} = 4.2 \times 10^{-3} \, \mathrm{K}^{-1}$

Substituting room-temperature reference values, the temperature-dependent resistivity functions are [1]:

$$
\begin{aligned}
\rho_{\mathrm{Cu}}(T) &= -3.5 \times 10^{-9} + 7.2 \times 10^{-11} \, T \quad (\Omega\cdot\mathrm{m}) \\
\rho_{\mathrm{Al}}(T) &= -6.57 \times 10^{-9} + 1.2 \times 10^{-10} \, T \quad (\Omega\cdot\mathrm{m})
\end{aligned} \tag{6.9}
$$

The filament resistance at temperature $T$ is obtained by substituting $\rho(T)$ into Eq. (6.5) for the driving coil or Eq. (6.7) for armature filaments.

### 6.4 Coupled Temperature-Circuit Iteration

At each time step $n$, for each armature filament $(i,j)$ [1]:

**1. Compute temperature increment** using the previous step's values:

$$
\boxed{T_{ij,n} = T_{ij,n-1} + \frac{I_{ij,n-1}^2 \, R_{ij,n-1} \, \Delta t}{m_{\mathrm{p}} \, c_{\mathrm{p}}(T_{ij,n-1})}} \tag{6.10}
$$

**2. Update material properties** at the new temperature:
- $c_{\mathrm{p}}(T_{ij,n})$ from Eq. (6.2) or (6.3)
- $\rho(T_{ij,n})$ from Eq. (6.9)
- $R_{ij,n}$ from Eq. (6.5) (driving coil) or Eq. (6.7) (armature filament)

**3. Proceed** to the next time step with updated resistances.

> **Figure 6.** Flowchart of armature temperature rise calculation. Source: [1] Fig. 3.

### 6.5 Key Simulation Results

From the 25-stage simulation in [1] (aluminum armature, mass = 1 kg, $C = 240 \, \mu\mathrm{F}$ per stage, $U = 4 \,\mathrm{kV}$, initial temperature $20^\circ\mathrm{C}$):

> **Table 3.** Simulation results with and without temperature rise consideration (from [1] Tables 1–2):

| Material | Temp. Rise Considered | Muzzle Velocity (m/s) | Max Temp. Rise (${}^{\circ}$C) | Efficiency (%) |
|:---:|:---:|:---:|:---:|:---:|
| Aluminum | Yes | 83.0 | 6.7 | 7.17 |
| Aluminum | No  | 87.5 | —   | 7.97 |
| Copper   | Yes | 103  | 14.0 | 11.05 |
| Copper   | No  | —    | —   | 12.15 |

**Key findings** [1]:

1. **Maximum temperature rise** occurs at the outermost surface of the armature tail (where eddy current density is highest due to the skin effect).
2. **Temperature decreases** axially from tail to front, with a secondary rise at the front tip.
3. **Copper vs. Aluminum at equal mass**: Copper has higher absolute temperature rise (due to larger induced eddy currents), but efficiency degrades less because $\beta_{\mathrm{Cu}} < \beta_{\mathrm{Al}}$ — copper's resistivity is less sensitive to temperature.
4. **Subdivision sensitivity**: Mesh density significantly affects temperature prediction but has modest influence on terminal velocity and efficiency (see [1] Tables 3–4).

> **Figure 7.** Temperature distribution in aluminum armature (25-stage simulation). Source: [1] Fig. 8.

> **Figure 8.** Velocity curves for copper and aluminum armatures. Source: [1] Fig. 12.

---

## 7. Efficiency Calculation

The energy conversion efficiency is defined as the ratio of the armature's muzzle kinetic energy to the total initial stored electrical energy [1]:

$$
\boxed{\eta = \frac{\frac{1}{2} m_{\mathrm{a}} v_{\mathrm{muzzle}}^2}{\frac{1}{2} n C U^2} = \frac{m_{\mathrm{a}} v_{\mathrm{muzzle}}^2}{n C U^2}} \tag{7.1}
$$

where $n$ is the total number of stages, $C$ is the capacitance per stage, and $U$ is the initial capacitor voltage per stage.

**Alternative efficiency form** (considering ohmic losses) [3]:

$$
\eta_{\mathrm{em}} = \frac{1}{1 + \frac{W_{\mathrm{da}}^{\mathrm{p} + \mathrm{b}}}{\Delta W_{\mathrm{k}}}} \tag{7.2}
$$

where $W_{\mathrm{da}}^{\mathrm{p} + \mathrm{b}}$ is the total ohmic loss (armature + barrel/driving coils) and $\Delta W_{\mathrm{k}}$ is the kinetic energy gain.

---

## 8. Complete Numerical Simulation Algorithm

### 8.1 Input Parameters

| Category | Parameters |
|:---|:---|
| **Driving coil(s)** | Inner radius $r_{\mathrm{di}}$, outer radius $r_{\mathrm{de}}$, length $l_{\mathrm{d}}$, turns $N_{\mathrm{d}}$, wire resistivity $\rho_{\mathrm{w}}$, wire cross-section $A_{\mathrm{w}}$, fill factor $k_{\mathrm{f}}$ |
| **Armature** | Inner radius $r_{\mathrm{ai}}$, outer radius $r_{\mathrm{ae}}$, length $l_{\mathrm{a}}$, mass $m_{\mathrm{a}}$, material (Cu or Al), initial position $x_0$, initial velocity $v_0$, axial divisions $m$, radial divisions $n$ |
| **Power supply** | Capacitance per stage $C$, initial voltage per stage $U_0$ |
| **Numerical** | Time step $\Delta t$, termination criteria |

> **Time step selection**: For forward Euler stability in an RLC-type system, $\Delta t$ should satisfy $\Delta t < 2 / \lambda_{\max}$ where $\lambda_{\max}$ is the largest-magnitude eigenvalue of the system matrix $([L]-[M_{\mathrm{I}}])^{-1}[R]$. A practical guideline is $\Delta t \lesssim \min(0.01 \times 2\pi\sqrt{L_{\mathrm{d}}C}, \; \delta l / v_{\max})$ — the smaller of 1% of the electrical oscillation period and the time for the armature to traverse one axial filament. Typical values range from $1\,\mu\mathrm{s}$ to $100\,\mu\mathrm{s}$. Verify convergence by halving $\Delta t$ and checking that results do not change significantly.

### 8.2 Precomputation (One-Time)

1. Driving coil resistance: $R_{\mathrm{d}} = \rho_{\mathrm{w}} \cdot k_{\mathrm{f}} \cdot \pi (r_{\mathrm{de}}^2 - r_{\mathrm{di}}^2) \cdot l_{\mathrm{d}} / A_{\mathrm{w}}^2$
2. Driving coil turns density: $n_{\mathrm{c}} = N_{\mathrm{d}} / [(r_{\mathrm{de}} - r_{\mathrm{di}}) \cdot l_{\mathrm{d}}]$
3. Driving coil self-inductance: $L_{\mathrm{d}} = 2\pi \mu_0 n_{\mathrm{c}}^2 r_{\mathrm{di}}^5 \, T(q, p)$ (Eq. 4.1)
4. For each armature filament $(i,j)$:
   - Inner/outer radius and axial position (Section 2.3)
   - Self-inductance $L_{ij}$ (Eq. 4.1 with appropriate parameters)
   - Initial resistance $R_{ij}$ (Eq. 6.6)
   - Mass $m_{\mathrm{p},ij} = \rho_{\mathrm{material}} \cdot \pi(r_{ij}^{\mathrm{(out)}2} - r_{ij}^{\mathrm{(in)}2}) \cdot \delta l$
5. Armature filament-to-filament mutual inductance matrix $[M]$ (Eq. 4.13, needed once since relative positions are fixed)
6. **Inter-coil mutual inductance** $M_{\mathrm{cc},ij}$ (for $i \neq j$): computed once at initialization using the same 4D integration (Eq. 4.13), applied to the cross-sections of driving coils $i$ and $j$. Since driving coils are fixed, these values are constant and populate the off-diagonal elements of $L_{\mathrm{c}}$ in Eq. (3.12).

### 8.3 Per-Timestep Computation

At each time step $n$:

1. **Update armature position**: $x_n = x_{n-1} + v_{n-1} \Delta t$
2. **Recompute** $[M_{\mathrm{I}}]_n$ and $[\dot{M}_{\mathrm{I}}]_n$ for all driving coil–armature filament pairs (Eqs. 4.13, 4.17)
3. **Solve circuit**: Compute $[\dot{I}]_n$ from Eq. (5.6), update $[I]_n$ (Eq. 5.7)
4. **Compute force**: $F_n$ from Eq. (5.8)
5. **Compute motion**: $a_n$, $v_n$, $z_n$ from Eq. (5.9)
6. **Compute temperature**: For each filament, update $T_{ij}$ from Eq. (6.7), then update $c_{\mathrm{p}}$, $\rho$, $R$ using Eqs. (6.2)–(6.6)
7. **Update capacitor voltages** (with crowbar diode per §3.4):
   - Diode-OFF stages: $U_{\mathrm{C},n}^{(i)} = U_{\mathrm{C},n-1}^{(i)} - (\Delta t / C_i) \cdot I_{\mathrm{d},n-1}^{(i)}$
   - If $U_{\mathrm{C},n}^{(i)} < 0$: set $U_{\mathrm{C},n}^{(i)} = 0$ and mark stage $i$ as diode-ON
   - Diode-ON stages: $U_{\mathrm{C},n}^{(i)} = 0$ (clamped)
8. **Trigger check**: If the armature center reaches the next coil's trigger position, activate the next stage (capacitor begins discharging). The optimal trigger position corresponds to the point where the mutual inductance gradient $\mathrm{d}M/\mathrm{d}x$ is near its maximum when the driving coil current peaks. In practice, this is approximately when the armature tail coincides with the driving coil's rear edge. The exact optimal position depends on armature velocity and circuit parameters — for a first implementation, place the trigger at $x_{\mathrm{trigger}}^{(i+1)} = x_{\mathrm{center}}^{(i+1)} - 0.5 l_{\mathrm{a}}$ (armature tail at driving coil center plane), then refine via parametric sweep.
9. **Termination check**: Stop the simulation when **either** (a) the current in all active driving coils drops below $I_{\mathrm{threshold}}$ (e.g., $10^{-3} \times I_{\mathrm{peak}}$) **and** the acceleration $|a| < a_{\mathrm{threshold}}$ (e.g., $10^{-3} \; \mathrm{m/s^2}$), **or** (b) the armature center position exceeds the barrel length (last coil center + coil length/2 + armature length/2). For multi-stage simulations, also terminate if no stages remain to be triggered and all active stage currents have decayed.

### 8.4 Computational Complexity

The dominant cost is the **4D integration** for $M$ and $\mathrm{d}M/\mathrm{d}x$:
- **Number of evaluations per timestep**: $n_{\mathrm{active}} \times (m \times n)$ driving coil–armature filament pairs.
- **Per evaluation**: $n_{\mathrm{quad}}^4$ calls to $M_{\text{loop}}$ (Eq. 4.8) or $\mathrm{d}M_{\text{loop}}/\mathrm{d}z$ (Eq. 4.15), each involving two elliptic integral evaluations.

**Optimization strategies**:
- **LRU caching**: Cache filament-level $M_{\text{loop}}$ and $\mathrm{d}M_{\text{loop}}/\mathrm{d}z$ results (same filament-pair geometries recur frequently).
- **Precompute $[M]$**: Inter-filament mutual inductance matrix is computed once — it does not change.
- **Exploit symmetry**: $M(a,b,d) = M(b,a,d)$ halves the number of unique filament-pair evaluations within the 4D integrand. For $[M_{\mathrm{I}}]$, the integration is over (driving coil $\times$ armature filament) which does not directly benefit from this symmetry in the outer loops, but the inner filament-pair calls benefit from caching.
- **Adaptive quadrature**: Use $\varepsilon_{\mathrm{abs}} = 10^{-6}$, $\varepsilon_{\mathrm{rel}} = 10^{-6}$ for 4D integrals, or use fixed-order Gauss-Legendre quadrature (Table 2) for deterministic cost.

---

## 9. Summary of Key Formulas

| Quantity | Formula | Eq. Ref |
|:---|:---|:---:|
| Axial force (pair) | $F = \frac{\mathrm{d}M}{\mathrm{d}x} I_{\mathrm{d}} I_{\mathrm{p}}$ | (1.2) |
| Axial force (all filaments) | $F = \sum_i \sum_j I_{\mathrm{d}i} I_{aj} \frac{\mathrm{d}M_{\mathrm{d}i,aj}}{\mathrm{d}z}$ | (5.1) |
| Motion | $m_{\mathrm{a}} \ddot{x} = F$, $\dot{x} = v$ | (5.2)–(5.3) |
| Circuit ODE | $[\dot{I}] = ([L]-[M_{\mathrm{I}}])^{-1}([U_{\mathrm{eff}}] + v[\dot{M}_{\mathrm{I}}][I] - [R][I] - [M][I])$ | (3.20) |
| Capacitor (diode OFF) | $U_{\mathrm{C}}(t) = U_0 - \frac{1}{C} \int I_{\mathrm{d}} \, \mathrm{d}t$ | (3.15) |
| Capacitor (diode ON) | $U_{\mathrm{C}} = 0$, RL decay | §3.4 |
| Self-inductance | $L = 2\pi\mu_0 n_{\mathrm{c}}^2 R_1^5 \, T(p,q)$ | (4.1) |
| $T(p,q)$ function | $T = \int_0^\infty U^2(1,p,x)(q x + e^{-q x} - 1)\mathrm{d}x$ | (4.3) |
| Filament M | $M_{\text{loop}} = \mu_0\sqrt{ab}\bigl[(\frac{2}{k}-k)K(k) - \frac{2}{k}E(k)\bigr]$ | (4.8) |
| Filament dM/dz | $\frac{\mathrm{d}M}{\mathrm{d}z} = \frac{\mu_0 k d}{4(1-k^2)\sqrt{ab}}[2(1-k^2)K(k) - (2-k^2)E(k)]$ | (4.15) |
| M symmetry | $M(a,b,d) = M(b,a,d) = M(a,b,-d)$ | (4.17) |
| dM/dz antisymmetry | $\frac{\mathrm{d}M}{\mathrm{d}z}(a,b,d) = -\frac{\mathrm{d}M}{\mathrm{d}z}(a,b,-d)$ | (4.18) |
| Coil M (4D) | $M = \frac{N_1 N_2}{16} \iiint_{-1}^1 M_{\text{loop}}' \, \mathrm{d}^4\mathbf{r}$ | (4.13) |
| Coil dM/dx (4D) | $\frac{\mathrm{d}M}{\mathrm{d}x} = \frac{N_1 N_2}{16} \iiint_{-1}^1 \frac{\mathrm{d}M_{\text{loop}}}{\mathrm{d}z}' \, \mathrm{d}^4\mathbf{r}$ | (4.16) |
| Driving coil resistance | $R_{\mathrm{d}} = \rho \cdot k\pi(r_{\mathrm{de}}^2 - r_{\mathrm{di}}^2) l_{\mathrm{d}} / S_{\mathrm{wire}}^2$ | (6.5) |
| Filament resistance | $R_{ij} = 2\pi\rho\bigl[\frac{mj}{l_{\mathrm{a}}} - \frac{m}{2l_{\mathrm{a}}} + \frac{mn r_{\mathrm{ai}}}{l_{\mathrm{a}}(r_{\mathrm{ae}}-r_{\mathrm{ai}})}\bigr]$ | (6.7) |
| Temperature update | $T_n = T_{n-1} + \frac{I_{n-1}^2 R_{n-1} \Delta t}{m_{\mathrm{p}} c_{\mathrm{p}}(T_{n-1})}$ | (6.10) |
| Specific heat (Al) | $c_{\mathrm{Al}} = 0.819 \exp(3.07\times 10^{-4} T)$ kJ/(kg·K) | (6.2) |
| Specific heat (Cu) | $c_{\mathrm{Cu}} = 0.333 \exp(3.917\times 10^{-4} T)$ kJ/(kg·K) | (6.3) |
| Resistivity | $\rho(T) = \rho_0(1 + \beta \Delta T)$, $\beta_{\mathrm{Cu}}=4.1\times10^{-3}$, $\beta_{\mathrm{Al}}=4.2\times10^{-3}$ K$^{-1}$ | (6.8) |
| Efficiency | $\eta = m_{\mathrm{a}} v^2 / (n C U^2)$ | (7.1) |
| Skin depth | $\delta = \sqrt{2/(\omega\sigma\mu)}$ | (2.1) |
| AGM for $K(k)$ | $K(k) = \pi/(2\cdot\text{AGM}(1,\sqrt{1-k^2}))$ | (4.21) |

---
