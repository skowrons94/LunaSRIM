# Physical models

lunasrim calculates stopping powers, ranges, energy loss and the transport of ions in amorphous targets with the models of SRIM-2013. The SRIM documentation describes these models in general terms [1, 2], but several details that change the numbers are only defined by the programs `TRIM.exe` and `SRModule.exe` themselves. We reconstructed them from the binaries and checked every routine against the original programs (see [validation.md](validation.md)). This document states the models as lunasrim implements them, with the approximations they contain. Section 1 fixes units and conventions. Sections 2 and 3 give the electronic and nuclear stopping powers, and Section 4 the stopping and range tables. Section 5 describes the TRIM Monte Carlo calculation. Sections 6 and 7 cover two analytic tools that are not part of SRIM: the energy loss and energy straggling of a beam crossing a layer, and Gaussian implantation profiles. Section 8 summarises what is needed to reproduce SRIM's numbers digit for digit.

Each model ends with an "Implementation" sentence that names the file and function in the `lunasrim` package.

## 1. Conventions and units

Ion energies are in keV at the interface and in eV inside TRIM. We write $E$ for the ion energy and $E_u = E/M_1$ for the energy per mass unit (keV/u), which sets the ion velocity. Lengths are in ångström (1 Å = 0.1 nm). Densities $\rho$ are in g/cm³. Masses are in atomic mass units (u). The ion mass $M_1$ defaults to the mass of the most abundant isotope, and the target masses $M_2$ default to the natural atomic weights, as in SRIM.

The stopping routines return the stopping cross section in eV/(10¹⁵ atoms/cm²). A target with stoichiometric numbers $s_i$ has atomic fractions $f_i = s_i / \sum_j s_j$, mean mass $\bar M_2 = \sum_i f_i M_{2,i}$ and atomic density

$$
N = \frac{N_A\,\rho}{\bar M_2}, \tag{1}
$$

where $N_A$ is Avogadro's number. TRIM uses $N_A = 0.6022$ in units of 10²⁴ mol⁻¹, which gives $N$ in atoms/Å³. The SR module uses $6.022\times10^{23}$ for the atomic density printed in the table header, $0.60223$ to convert stopping cross sections into stopping powers, and $0.6023$ in the range integration. We keep the three values because they change the last printed digits. A stopping cross section $S$ in eV/(10¹⁵ atoms/cm²) corresponds to the stopping power $10\,N S$ in eV/Å when $N$ is in atoms/Å³.

Reduced energies use the universal screening length of Ziegler, Biersack and Littmark [1],

$$
a = \frac{0.8853\,a_0}{Z_1^{0.23} + Z_2^{0.23}}, \tag{2}
$$

where $a_0 = 0.5292$ Å is the Bohr radius and $Z_1$, $Z_2$ are the atomic numbers of the ion and the target atom. The SR module writes the product as 0.4685 Å. The squared elementary charge is $e^2 = 14.4$ eV Å in most places; TRIM's pair constants use 14.41 eV Å. All elements from $Z = 1$ to $Z = 92$ are available as ions and as target atoms.

## 2. Electronic stopping

The electronic stopping cross section of an ion in an element is an empirical fit whose coefficients are read from the SRIM-2013 file `SCOEF03.DAT` [2]. The form of the fit depends on the ion: protons, helium and lithium ions have their own coefficient sets, and heavier ions are scaled from helium. Above the validity range of the fits the stopping follows the Bethe formula. Phase and compound corrections multiply the result.

### 2.1 Light ions (H, He, Li)

For protons, helium and lithium ions the stopping cross section between $E_u = 2$ keV/u and an upper limit $E_\mathrm{max}$ is

$$
S_\mathrm{fit}(E_u) = \frac{S_\mathrm{low}\,S_\mathrm{high}}{S_\mathrm{low} + S_\mathrm{high}}, \qquad
S_\mathrm{low} = c_1 E_u^{c_2} + c_3 E_u^{c_4} + c_9 E_u^{1/2}, \qquad
S_\mathrm{high} = \frac{c_5}{E_u^{c_6}} \left[ \ln\!\left( \frac{c_7}{E_u} + c_8 E_u \right) + \frac{c_{10}}{E_u^{0.2}} \right], \tag{3}
$$

where $c_1, \dots, c_{10}$ are the ten coefficients of the target element for the ion (columns 11 to 20 of `SCOEF03.DAT` for protons, 21 to 30 for helium and 31 to 40 for lithium) and $E_u$ is in keV/u. The limit is $E_\mathrm{max} = 1000$ keV/u for protons and 2000 keV/u for helium and lithium. Below 2 keV/u the stopping is taken proportional to the ion velocity, $S = S_\mathrm{fit}(2)\,(E_u/2)^{1/2}$. Between $E_\mathrm{max}$ and $1.2\,E_\mathrm{max}$ it is interpolated linearly between $S_\mathrm{fit}(E_\mathrm{max})$ and the Bethe value of Eq. (4) at $1.2\,E_\mathrm{max}$. Above $1.2\,E_\mathrm{max}$ it is the Bethe value.

### 2.2 High energies

At high velocity the stopping is the Bethe formula with shell, Barkas and Bloch terms:

$$
S_\mathrm{B} = \frac{4\pi e^4}{m_e c^2}\,\frac{Z_1^2 Z_2}{\beta^2}\,L, \qquad
L = \ln\frac{2 m_e c^2 \beta^2}{1-\beta^2} - \beta^2 - C_1 - C_2 - \ln I + L_1 + L_2, \tag{4}
$$

where $4\pi e^4/(m_e c^2) = 1/1961.16$ eV/(10¹⁵ atoms/cm²), $2 m_e c^2 = 1.022\times10^6$ eV, $\beta^2 = 1 - (1 + E_u/931494\ \mathrm{keV})^{-2}$, $C_1$ and $C_2$ are polynomials in $\ln(E_u / 1000\ \mathrm{keV})$ with 8 and 11 coefficients of the target element, $I$ is the mean excitation energy in eV, $L_1$ is the Barkas term and $L_2$ the Bloch term. For a gas target $C_1$ is set to zero and $I$ is replaced by the gas value of `SCOEF03.DAT` (one value for protons and one for all other ions). The Barkas term is

$$
L_1 = \frac{Z_1\,u\,q}{u + q}, \qquad u = \frac{E_u}{1000\ \mathrm{keV}}, \tag{5}
$$

where $q = 1.5\,E_u^{-0.4} + 45000/(Z_2 E_u^{1.6})$ for solids and $q = 1.35\,E_u^{-0.4}$ for gases, with $E_u$ in keV/u. The Bloch term is

$$
L_2 = -y \left[ 1.20206 - y \left( 1.042 - 0.8549\,y + 0.343\,y^2 \right) \right], \qquad y = \frac{(Z_1/137.036)^2}{\beta^2}, \tag{6}
$$

where 1.20206 is $\zeta(3)$, the leading coefficient of Bloch's correction.

### 2.3 Heavy ions

For ions with $Z_1 \ge 4$ the stopping is scaled from the helium stopping at the same velocity. Between 2 and 30000 keV/u,

$$
S(E_u) = S_\mathrm{He}(E_u)\,\frac{Z_1^2}{4}\,P_\mathrm{ion}(\ln E_u)\,P_\mathrm{tgt}(\ln E_u), \tag{7}
$$

where $S_\mathrm{He}$ is the helium stopping of Section 2.1 (including its Bethe part above 2400 keV/u), $P_\mathrm{ion}$ is a polynomial of degree 12 with coefficients of the ion (columns 41 to 53, described in the file layout as an effective-charge polynomial) and $P_\mathrm{tgt}$ is a polynomial of degree 12 with coefficients of the target element (columns 54 to 66). Below 2 keV/u the stopping is again proportional to velocity. Above 30000 keV/u it follows the Bethe formula with a high-energy effective charge,

$$
S(E_u) = S(30000)\,\frac{\gamma(E_u)\,S_\mathrm{B,He}(E_u)}{\gamma(30000)\,S_\mathrm{B,He}(30000)}, \qquad
\gamma(E_u) = h_1 + \frac{1 - h_1}{1 + \exp(-d/w)} + h_5 \exp\left(-\left|1.3 + h_2 - \log_{10} E_u\right|\right), \tag{8}
$$

where $S_\mathrm{B,He}$ is Eq. (4) for $Z_1 = 2$, $h_1, \dots, h_5$ are five parameters of the ion (columns 89 to 93), $d = \log_{10} E_u - h_2$, and $w = h_4$ for $d \ge 0$ and $w = h_3$ otherwise. Above $10^7$ keV/u the routine returns zero.

### 2.4 Phase and compound corrections

SRIM multiplies the stopping of each target element by a factor that is fully applied at low velocity and fades out at high velocity:

$$
F(E_u) = 1 + \frac{F_0 - 1}{1 + \exp\left[1.48\left(\sqrt{E_u/25} - c\right)\right]}, \tag{9}
$$

where $E_u$ is in keV/u, $c = 7$ for protons and $c = 4$ for all other ions, and $F_0 = f_\mathrm{phase}\,C$ is the product of a phase factor and the compound correction $C$ given by the user (a value of 0 or less is replaced by 1). The factor reaches half of its low-energy deviation at $E_u = 25c^2$, which is 1225 keV/u for protons and 400 keV/u for heavier ions. The phase factor is 1 for the noble gases and for elements in their natural phase. It is 0.92 for an element that SRIM classes as gaseous (H, N, O, F, Cl, Br, I, At) in a solid target, and $1/0.92$ for any other element in a gas target. The value 0.92 is SRIM's.

TRIM decides the phase of each element once: an element that appears in any gas layer is treated as a gas atom in every layer. lunasrim reproduces this.

The exponential in Eq. (9) is evaluated for every target element, also when $F_0 = 1$. Its argument exceeds the range of a double-precision number above $E_u = 5.92\times10^6$ keV/u for protons and $5.85\times10^6$ keV/u for heavier ions. `SRModule.exe` and `TRIM.exe` then stop with Visual Basic run-time error 6 (overflow). lunasrim raises `SrError` at the same energies.

### 2.5 Bragg's rule

The electronic stopping of a compound is the sum of the elemental stopping cross sections weighted with the atomic fractions (Bragg's rule),

$$
S_e(E) = \sum_i f_i\,F_i(E_u)\,S_{e,i}(E), \tag{10}
$$

where $S_{e,i}$ is the stopping in element $i$ from Sections 2.1 to 2.3 and $F_i$ the factor of Eq. (9) for that element. Bragg's rule neglects the effect of chemical bonds on the target electrons. The compound correction $C$ accounts for it on average. SRIM's compound dictionary (356 entries, `Compound.dat`) provides values of $C$ for many materials. In the list mode of the SR module (stopping powers at given energies) and in `lunasrim.stopping.compound_stopping`, `SRModule.exe` applies $C$ a second time: the sum of Eq. (10) is multiplied by Eq. (9) with $F_0 = C$ and $c = 7$ when $E < 10$ MeV. We reproduce this. The table mode of the SR module and TRIM apply $C$ once, through $F_i$.

Implementation: `lunasrim/stopping.py`, functions `electronic_stopping_element` (Eqs. 3 to 8), `fade_correction` and `phase_factor` (Eq. 9), `electronic_stopping` and `compound_stopping` (Eq. 10).

## 3. Nuclear stopping

The nuclear stopping describes the energy given to target atoms in elastic collisions. With the reduced energy

$$
\varepsilon = \frac{32.53\,M_2\,E}{Z_1 Z_2 (M_1 + M_2)\left(Z_1^{0.23} + Z_2^{0.23}\right)}, \tag{11}
$$

where $E$ is in keV, the reduced nuclear stopping is

$$
s_n(\varepsilon) = \frac{\ln(1 + b_1 \varepsilon)}{2\left(\varepsilon + b_2\,\varepsilon^{b_3} + b_4\,\varepsilon^{1/2}\right)} \quad (\varepsilon \le 30), \qquad
s_n(\varepsilon) = \frac{\ln \varepsilon}{2\varepsilon} \quad (\varepsilon > 30), \tag{12}
$$

and the nuclear stopping cross section is

$$
S_n = \frac{8.462\,Z_1 Z_2 M_1}{(M_1 + M_2)\left(Z_1^{0.23} + Z_2^{0.23}\right)}\,s_n(\varepsilon), \tag{13}
$$

in eV/(10¹⁵ atoms/cm²), where $b_1, \dots, b_4$ are fit coefficients. With $b = (1.1383, 0.01321, 0.21226, 0.19593)$, Eq. (12) is the universal nuclear stopping of Ziegler, Biersack and Littmark [1]. SRIM-2013 also ships coefficients for each of the 8464 ion-target pairs in `SNUC03.dat`. For 50 keV Ar in Si the pair fit gives 118.2 eV/(10¹⁵ atoms/cm²) and the universal form 110.7 eV/(10¹⁵ atoms/cm²), 6 % less.

The programs use the two forms in different places, and so does lunasrim. The list mode of the SR module, `lunasrim.stopping.nuclear_stopping` and the energy-loss calculation of Section 6 use the pair coefficients. The table mode of the SR module uses the universal form inside its range integration (Section 4), with $b_4 = 0.19594$ in the integration steps and 0.19593 at the starting point, as in `SRModule.exe`. TRIM does not use a nuclear stopping power: it follows the individual collisions (Section 5.3). For a compound, the nuclear stopping is summed with Bragg's rule without corrections.

Implementation: `lunasrim/stopping.py`, functions `reduced_energy` (Eq. 11), `nuclear_stopping` (pair coefficients) and `universal_nuclear_stopping`; `lunasrim/sr.py`, `StoppingTable._integrate`.

## 4. Stopping and range tables

The SR module tabulates electronic and nuclear stopping, the projected range and the longitudinal and lateral range straggling between two energies.

The energy grid has 26 points per decade with the mantissas 1.0, 1.1, ..., 1.8, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0 and 9.0, starting at 1 eV and ending at the first grid point above $E_\mathrm{max}$. The electronic stopping is evaluated on this grid with Eq. (10) and interpolated linearly in between.

The ranges are obtained with a PRAL-type integration of the first moments of the ion distribution, after Biersack, as in `SRModule.exe`. For each target element $k$ we define $A_k = M_{2,k}/M_1$, $\gamma_k = 4A_k/(1+A_k)^2$, the nuclear stopping $S_{n,k}(E)$ of Eq. (13) in the universal form, and a second nuclear moment

$$
S_{t,k}(E) = \frac{\pi a_k^2\,\gamma_k^2\,N}{\varepsilon_{f,k}^2}\;
\frac{1}{4 + 6.684\,\varepsilon^{-1.0494} + 0.197\,\varepsilon^{-1.6991}}, \tag{14}
$$

where $a_k$ is the screening length of Eq. (2), $\varepsilon_{f,k} = a_k A_k / [(1+A_k) Z_1 Z_{2,k} e^2]$ is the reduced energy per eV and $\varepsilon = \varepsilon_{f,k} E$. $S_{t,k}$ has the dimension of a mean square energy transfer per unit path length (eV²/Å). With the total stopping $S = S_e + S_n$ in eV/Å and the weighted sums $\nu = \sum_k f_k A_k S_{n,k}$, $\tau = \sum_k f_k A_k S_{t,k}$ and $\mu = \sum_k f_k (1 - A_k)^2 S_{t,k}$, each energy step $\Delta E$ updates the projected range $R_p$, the mean square total range $\langle R^2 \rangle$ and the mean square radial displacement $\Lambda$ as

$$
\Delta R_p = \frac{\Delta E}{E}\,\frac{4E^2 - (2E\nu + \tau)\,R_p}{4ES - 2\tau}, \qquad
\Delta\langle R^2\rangle = \frac{2 R_p\,\Delta E}{S}, \qquad
\Delta\Lambda = \frac{\Delta E}{E}\,\frac{\left(\langle R^2\rangle - 2\Lambda\right)(4E\nu - \mu) - 8E\tau R_p/S}{4ES - 8\tau}, \tag{15}
$$

where all quantities are evaluated at the end of the step. When $1 - 2\tau/(ES) < 0.01$ the last expression is replaced by $\Lambda = \langle R^2\rangle/2 - E^2 R_p/(2E\nu - \mu/2)$. Without nuclear stopping, $\nu = \tau = \mu = 0$ and Eq. (15) reduces to $dR_p/dE = 1/S$. The tabulated longitudinal straggling is $(\langle R^2\rangle - \Lambda - R_p^2)^{1/2}$ and the lateral straggling is $(\Lambda/2)^{1/2}$, the spread along one transverse axis. The integration starts at 1.1 eV from a low-energy solution in which $R_p$, $\langle R^2\rangle$ and $\Lambda$ are fixed by the ratio of nuclear and electronic stopping, and it divides every grid interval into 50 equal steps. The energy is accumulated step by step, so that the printed energies carry the rounding of 50 additions (for example `999.999 eV`). A row is written at the end of each grid interval above $E_\mathrm{min} - 0.2$ eV. A negative variance stops `SRModule.exe` with a run-time error; lunasrim raises `SrError` in that case.

Eq. (15) treats the ion distribution through its first two moments in a target of uniform composition. It does not describe the shape of the distribution beyond its width, and it neglects the energy straggling of the ion along its path. For a comparison with TRIM see Section 7.

The tables can be printed in the eight stopping units of SRIM: eV/Å, keV/µm, MeV/mm, keV/(µg/cm²), MeV/(mg/cm²), keV/(mg/cm²), eV/(10¹⁵ atoms/cm²) and reduced (LSS) units.

Implementation: `lunasrim/sr.py`, class `StoppingTable` (methods `calculate` and `_integrate`) for tables and `StoppingList` for lists of energies.

## 5. TRIM Monte Carlo

TRIM follows ions one by one through a target made of layers. Each layer is amorphous and homogeneous, with a composition, a density and a phase (solid or gas). The ion moves on straight free flight paths. At the end of each path it loses energy to the electrons continuously and collides elastically with one target atom. The calculation follows the method of Biersack and Haggmark [3] and of Ziegler, Biersack and Littmark [1] in the form of `TRIM.exe` of SRIM-2013. Crystal structure, channelling, temperature and target changes during the irradiation (dose effects) are not modelled.

### 5.1 Target constants

For every layer TRIM uses the mean atomic number $\bar Z_2$ and the mean mass $\bar M_2$ of the layer. $\bar Z_2$ is an integer: it is rounded to the nearest integer after the contribution of each element is added, as in `TRIM.exe`. The reduced energy of the ion in the layer is $\varepsilon = \varepsilon_f E$ with $\varepsilon_f = a \bar M_2 / [\bar Z_2 Z_1 e^2 (M_1 + \bar M_2)]$ and $a$ from Eq. (2) with $Z_2 = \bar Z_2$. The mean interatomic distance is $L_0 = N^{-1/3}$.

The electronic stopping power of each atom type (the ion and each target element as a recoil) is tabulated per layer with Eq. (10), multiplied by $10N$ to give eV/Å, on 1000 energies $E_k = \exp[G(k-1)]$ eV with $G = \ln(1000 E_0 + 1)/999$, where $E_0$ is the initial ion energy in keV. The table therefore extends from 1 eV to $E_0$ + 1 eV and is interpolated linearly.

### 5.2 Free flight path and impact parameter

The free flight path follows from a maximum impact parameter

$$
p_\mathrm{max} = \frac{a}{\tilde\varepsilon + \tilde\varepsilon^{1/2} + 0.125\,\tilde\varepsilon^{0.1}}, \qquad
\tilde\varepsilon = \left[\frac{5\,\varepsilon_f\,(1 + A)^2}{4A}\,\varepsilon\right]^{1/2}, \tag{16}
$$

where $A = M_1/\bar M_2$ and the factor 5 is in eV. The path length is

$$
L = \frac{1}{\pi p_\mathrm{max}^2 N}. \tag{17}
$$

The first path of an ion from the surface is $L_0$. In monolayer mode (calculation type 3) every path is $L_0$. With full damage cascades (type 2) every path is $L$ multiplied by a random number $R \in (0, 1)$ and is at least $L_0$. In quick Kinchin-Pease mode (type 1) the path is $L$. When a path would cross a layer boundary, the surface or the edge of the plot window, it is shortened by repeated multiplication with random numbers, with a lower limit of $L_0$ or of one depth bin. After such a cut, and after the first path, the next two paths are also multiplied by a random number and are at least $L_0$.

At the end of the path a target atom is chosen at random with the probabilities $f_i$ of the layer, and the impact parameter is

$$
p = \left(\frac{R}{\pi L N}\right)^{1/2}, \tag{18}
$$

where $R$ is a new random number. This samples the impact parameter uniformly over the area $1/(LN)$ that belongs to one target atom along the path.

### 5.3 Elastic scattering

The ion-atom interaction is the ZBL universal potential [1],

$$
V(r) = \frac{Z_1 Z_2 e^2}{r}\,\Phi(r/a), \qquad
\Phi(x) = 0.18175\,e^{-3.1998x} + 0.50986\,e^{-0.94229x} + 0.28022\,e^{-0.4029x} + 0.02817\,e^{-0.20162x}, \tag{19}
$$

where $a$ is the screening length of Eq. (2) for the pair. With the reduced energy $\varepsilon = a E / [Z_1 Z_2 e^2 (1 + M_1/M_2)]$ and the reduced impact parameter $b = p/a$, the centre-of-mass scattering angle $\theta$ is obtained as follows. For $\varepsilon > 10$ TRIM uses

$$
\sin^2\frac{\theta}{2} = \frac{1}{1 + \left[1 + b(1+b)\right](2\varepsilon b)^2}. \tag{20}
$$

For $\varepsilon \le 10$ it uses the MAGIC approximation of Biersack and Haggmark [3]. The reduced distance of closest approach $r_0$ is the root of $1 - \Phi(r)/(\varepsilon r) - b^2/r^2 = 0$, found by Newton iteration (at most 14 iterations, relative step below $3\times10^{-4}$). The scattering angle is then

$$
\cos\frac{\theta}{2} = \frac{b + \rho + \delta}{r_0 + \rho}, \qquad
\rho = -\frac{2\left[\varepsilon - V_r(r_0)\right]}{V_r'(r_0)}, \qquad
\delta = \frac{\Gamma\,\alpha\,(r_0 - b)}{1 + \Gamma}, \tag{21}
$$

where $V_r(r) = \Phi(r)/r$ is the reduced potential, $\alpha = 2\varepsilon\left(1 + 0.99229\,\varepsilon^{-1/2}\right) b^{\kappa}$ with $\kappa = (0.011615 + \varepsilon^{1/2})/(0.0071222 + \varepsilon^{1/2})$, and $\Gamma = \left[(1 + \alpha^2)^{1/2} - \alpha\right](9.3066 + \varepsilon)/(14.813 + \varepsilon)$. The five constants are the MAGIC coefficients for the universal potential [1].

The energy transferred to the target atom is $T = \gamma E \sin^2(\theta/2)$ with $\gamma = 4M_1M_2/(M_1+M_2)^2$. The ion is deflected by $\psi = \arctan[\sin\theta/(\cos\theta + M_1/M_2)]$ in the laboratory frame, the recoil leaves at $\arctan[\sin\theta/(1-\cos\theta)]$, and the azimuth is $2\pi R$. The ion loses $T$ after its electronic loss on the preceding path.

Implementation: `lunasrim/trim.py`, `TrimEngine.flight_path` (Eqs. 16, 17), `TrimEngine._run_ion` (Eq. 18), `TrimEngine.magic` (Eqs. 19 to 21) and `TrimEngine.new_direction`; constants in `lunasrim/trim_setup.py`, class `TrimTarget`.

### 5.4 Electronic energy loss and its straggling

On a path of length $L$ the ion loses

$$
\Delta E = L\,S_e(E) + \delta, \qquad \delta = \left(R - \tfrac12\right) w(E)\,L^{1/2}, \tag{22}
$$

where $S_e$ is the tabulated stopping power in eV/Å, $R$ a random number, and $|\delta|$ is limited to $L S_e$. The fluctuation $\delta$ is uniformly distributed with variance $w^2 L/12$. The width $w$ is

$$
w(E) = \begin{cases} kE, & kE \le w_\mathrm{B}, \\ w_\mathrm{B}\,\dfrac{1 - \beta^2/2}{1 - \beta^2}, & kE > w_\mathrm{B}, \end{cases} \qquad
w_\mathrm{B} = 177\,Z_1\left(\bar Z_2 N\right)^{1/2}, \qquad
k = \left[\frac{3.2 \times 12\,N\,(Z_1 + \bar Z_2)^{8/3}}{1000\,M_1}\right]^{1/2}, \tag{23}
$$

where $E$ is in eV, $N$ in atoms/Å³, $w$ in eV/Å$^{1/2}$, $\bar Z_2$ is the integer mean atomic number of the layer and $\beta^2 = 1 - [1 + E/(931.494\,M_1\ \mathrm{MeV})]^{-2}$. The variance per unit length, $w_\mathrm{B}^2/12 = 2610.75\,Z_1^2 \bar Z_2 N$ eV²/Å, is Bohr's straggling $4\pi e^4 Z_1^2 Z_2 N$ [4] (2605.63 $Z_1^2 Z_2 N$ eV²/Å with $e^2 = 14.39964$ eV Å) to 0.2 %. The low-energy branch $kE$ acts only at low energy: in silicon it applies below 92 eV for protons and below 3.8 keV for argon ions. Recoils use the same expression with their own $Z_1$, $M_1$ and stopping table. The ion stops when its energy falls to 1 eV.

### 5.5 Recoils and damage

A collision that transfers less than 1 eV gives its energy to phonons. Otherwise the struck atom becomes a recoil, which is treated in one of two ways.

Quick Kinchin-Pease calculation (type 1). TRIM does not follow the recoil. It splits the recoil energy $T$ into a damage energy $\nu$ and an ionisation part $T - \nu$ with the partition function used by Norgett, Robinson and Torrens [5],

$$
\nu = \frac{T}{1 + k_d\,g(\varepsilon_d)}, \qquad g(\varepsilon_d) = \varepsilon_d + 0.4\,\varepsilon_d^{3/4} + 3.4\,\varepsilon_d^{1/6}, \qquad
\varepsilon_d = \frac{0.01\,T}{\bar Z_2^{7/3}}, \qquad k_d = \frac{0.1334\,\bar Z_2^{2/3}}{\bar M_2^{1/2}}, \tag{24}
$$

where $T$ is in eV and $\bar Z_2$, $\bar M_2$ are the layer means. Norgett et al. write $\varepsilon_d = T/(86.93\,Z^{7/3})$, that is 0.0115 instead of 0.01. The number of vacancies is then

$$
N_v = \begin{cases} 0, & \nu < E_d, \\ 1, & E_d \le \nu < 2.5\,E_d, \\ 0.8\,\nu/(2E_d), & \nu \ge 2.5\,E_d, \end{cases} \tag{25}
$$

where $E_d$ is the displacement energy of the struck atom. This is the modified Kinchin-Pease estimate [5, 6]. The phonon tally receives $\nu - N_v E_L$, where $E_L$ is the lattice binding energy of the struck atom.

Full damage cascades (type 2) and monolayer steps (type 3). TRIM follows every recoil and every secondary recoil. An atom that receives $T > E_d$ leaves its site, which creates a vacancy, and starts with energy $T - E_L$ (secondary recoils within $2L_0$ of the surface keep the full $T$). Recoils move with flight paths $L_0$ and impact parameters $p = L_0\,(R/\pi)^{1/2}$, lose energy to electrons with Eqs. (22) and (23), and are followed down to 1 eV. Near the surface the impact parameters of light recoils ($Z < 21$, depth below $2L_0$) are modified, and atoms below $E_d$ within a random fraction of $9L_0$ of the surface are followed to allow sputtering. When a moving atom comes to rest and the last atom it struck was displaced and is of the same element, the event is counted as a replacement collision: the moving atom fills the vacancy and the vacancy count is reduced by one. The random number generator is saved before a cascade and restored afterwards, so that the ion trajectory does not depend on the cascade.

Types 1 and 2 nevertheless give different ion trajectories, because their flight paths differ (Section 5.2). The ion distributions agree within statistics: for 1000 ions of 50 keV Ar in Si the mean range is 586.4 Å with type 1 and 586.5 Å with type 2.

### 5.6 Sputtering, transmission and backscattering

A recoil that reaches the surface ($x < 0$) escapes if the energy of its motion normal to the surface exceeds the surface binding energy $E_s$, that is if $E\cos^2\alpha \ge E_s$, where $\alpha$ is the angle to the surface normal. Otherwise it is reflected back into the target. This is a planar surface barrier. The sputtering yield is the number of escaping atoms per ion. A recoil that leaves the back of the target is not followed further. An ion that leaves the front surface is backscattered and one that leaves the back is transmitted; their energies, positions and directions are written to `BACKSCAT.txt` and `TRANSMIT.txt`.

Missing displacement, lattice and surface energies of a TRIM.IN file take TRIM's defaults: 20 eV displacement energy for the first element; the others take its displacement energy, 2 eV lattice binding and 4 eV surface binding energy. The pysrim-compatible `Layer` class uses 25 eV, 3 eV and the heat of sublimation of `SCOEF03.DAT`.

Implementation: `lunasrim/trim.py`, `TrimEngine.eloss` (Eqs. 22, 23), `TrimEngine.recoil`, `on_kinchin_pease` (Eqs. 24, 25), `full_cascade`, `on_recoil_stop` and `replacement` (replacement collisions) and `on_recoil_surface` (sputtering).

### 5.7 Statistics and output distributions

The depth between the plot window limits $X_\mathrm{min}$ and $X_\mathrm{max}$ (by default the whole target) is divided into 100 bins. Every energy loss is distributed over the bins crossed by the path. TRIM writes the distributions of stopped ions and recoils (`RANGE.txt`), the ionisation (`IONIZ.txt`), the phonon energy (`PHONON.txt`), the energy given to recoils (`E2RECOIL.txt`), vacancies and replacement collisions (`VACANCY.txt`, `NOVAC.txt`) and the lateral spread (`LATERAL.txt`). The mean projected range is the mean depth $\langle x\rangle$ of the stopped ions, and the range straggling is $(\langle x^2\rangle - \langle x\rangle^2)^{1/2}$. The lateral range is $\langle|y|\rangle$ and the lateral straggling $\langle y^2\rangle^{1/2}$; the radial values use $r = (y^2 + z^2)^{1/2}$. The statistical uncertainty of these numbers scales as $n^{-1/2}$ with the number of ions $n$. For 50 keV Ar in Si with 1000 ions the standard error of the mean range is 7 Å, or 1.2 %.

Implementation: `lunasrim/trim.py`, `TrimEngine.tally_ion_path` and `ion_stats`; `lunasrim/tally.py`; `lunasrim/output.py`.

## 6. Energy loss and energy straggling in a layer

For an ion beam that crosses a layer, the mean exit energy and the energy spread can be calculated without a Monte Carlo. We integrate the mean energy $E(x)$ and the energy variance $\Omega^2(x)$ along the path,

$$
\frac{dE}{dx} = -S(E), \qquad \frac{d\Omega^2}{dx} = W(E) - 2\,S'(E)\,\Omega^2, \tag{26}
$$

where $S = S_e + S_n$ is the total stopping power of the SR module in list mode (Eqs. 10 and 13, with the pair coefficients of `SNUC03.dat` and the compound correction), $S' = dS/dE$ and $W$ is the variance of the electronic energy loss per unit path length. The second term of the variance equation describes how the slope of the stopping power changes an existing spread: where $S$ increases with decreasing energy ($S' < 0$) the slower ions lose more energy and the spread grows, and where $S' > 0$ it narrows.

We give $\Omega$ for two models of $W$. The default is TRIM's model of Eq. (23),

$$
W(E) = \frac{w(E)^2}{12}, \tag{27}
$$

which includes TRIM's relativistic factor and its integer mean atomic number $\bar Z_2$. At non-relativistic energies above the low-energy branch it is Bohr's straggling within 0.2 %. The second model is Bohr's formula [4],

$$
W_\mathrm{B} = 4\pi e^4 Z_1^2 Z_2 N, \tag{28}
$$

with the exact stoichiometric mean $Z_2 = \sum_i f_i Z_{2,i}$ and $e^2 = 14.39964$ eV Å. The two differ by the 0.2 % of Section 5.4, the relativistic factor, and the rounding of $\bar Z_2$ in compounds.

We tabulate $S_e$ and $S_n$ at 200 energies per decade between $\max(10^{-3}E_0, 1\ \mathrm{eV})$ and $1.02\,E_0$ and interpolate them linearly in $\ln S$ versus $\ln E$. The derivative $S'$ is a central difference over $\pm0.1$ % in energy. Eq. (26) is integrated with a fourth-order Runge-Kutta method in 500 equal steps. If the mean energy falls below 1 % of the initial energy, the integration stops and the result is flagged as stopped. The electronic and nuclear energy losses are separated with the ratio $S_e/S_n$ at the middle of each step. The FWHM is $2(2\ln 2)^{1/2}\,\Omega = 2.3548\,\Omega$. For a layer tilted by an angle $\phi$ to the beam, the path length is the thickness divided by $\cos\phi$.

Eq. (26) neglects three effects that a TRIM calculation includes or that real beams show.

- Nuclear energy-loss straggling. Nuclear collisions enter the mean energy through $S_n$ but not the spread. Single collisions with large energy transfer give a low-energy tail in the exit spectrum, which is not Gaussian and which a TRIM calculation shows (`TRANSMIT.txt`).
- Path-length straggling. Angular scattering makes the paths of individual ions longer than the layer thickness by different amounts. This adds to the spread for thick layers and heavy targets.
- Charge-state fluctuations and target-electron binding. $W$ contains no charge-exchange term, and Bohr's formula treats the target electrons as free, which overestimates the straggling at low velocity. TRIM makes the same approximation.

We compare the analytic $\Omega$ with the spread of lunasrim's TRIM for five ion-layer pairs. Both calculations use the defaults of `lunasrim.material` and `lunasrim.resolve_ion`: the ion masses 1.008 u (H), 4.003 u (He) and 12.000 u (C), and the densities 2.3212 g/cm³ (Si), 2.253 g/cm³ (C), 19.311 g/cm³ (Au), 0.001251 g/cm³ (N₂ gas) and 8.8955 g/cm³ (Ni). The TRIM runs are `lunasrim.trim(..., ions=4000, seed=7, calculation=1)`. The TRIM width is half the interval between the 16 % and 84 % quantiles of the energies of the transmitted ions, which excludes the nuclear tail.

| Ion, energy | Layer | Exit energy, analytic / TRIM median (keV) | $\Omega$, analytic / TRIM (keV) |
|---|---|---|---|
| H, 400 keV | Si, 1 µm | 330.5 / 331.6 | 4.48 / 4.55 |
| H, 100 keV | C, 300 nm | 50.1 / 50.1 | 2.22 / 2.38 |
| He, 2 MeV | Au, 300 nm | 1799.8 / 1799.8 | 12.27 / 12.62 |
| He, 300 keV | N₂ gas, 200 µm | 258.5 / 258.2 | 2.72 / 2.91 |
| C, 5 MeV | Ni, 500 nm | 3486.2 / 3485.1 | 33.58 / 34.95 |

The exit energies agree within 1.1 keV. The TRIM widths are larger by 1.6 %, 7.0 %, 2.8 %, 6.9 % and 4.1 %. For a normal distribution the standard error of this quantile estimate of the width is 1.5 % for 4000 ions, so the difference exceeds the statistical uncertainty in four of the five cases. This is consistent with the path-length and nuclear contributions that Eq. (26) omits.

The integer mean $\bar Z_2$ of TRIM's model matters for compounds of light elements. For 2 MeV He through 2 µm of Mylar at 45° (`lunasrim.material('Mylar')`, 1.397 g/cm³), the stepwise rounding gives $\bar Z_2 = 4$ instead of the stoichiometric mean 4.55, and the TRIM-model width is 11.71 keV against 12.47 keV with Bohr's formula.

Implementation: `lunasrim/straggling.py`, function `energy_loss` (Eqs. 26 to 28).

## 7. Implantation profiles

For ion implantation we approximate the depth distribution of the implanted atoms by a Gaussian,

$$
n(x) = \frac{D}{\sqrt{2\pi}\,\Delta R_p}\,\exp\left[-\frac{(x - R_p)^2}{2\,\Delta R_p^2}\right], \tag{29}
$$

where $n$ is the concentration in atoms/cm³, $D$ the dose in atoms/cm², $x$ the depth, $R_p$ the projected range and $\Delta R_p$ the longitudinal straggling of the SR table at the implantation energy (Section 4). The peak concentration is $D/(\sqrt{2\pi}\,\Delta R_p)$. Profiles of several energies are added.

Eq. (29) is a two-moment approximation and is poor where the distribution is skewed or has long tails. It also neglects sputtering of the surface, the change of the target composition with dose, channelling in crystals, and diffusion. Two TRIM runs with 1000 ions illustrate the size of these differences.

- 50 keV Ar in Si: TRIM gives a mean range of 586 ± 7 Å, a straggling of 227 Å, a skewness of 0.10 ± 0.08 and a kurtosis of 2.86. The SR table gives $R_p = 561$ Å and $\Delta R_p = 208$ Å. The shape is close to a Gaussian, and the table values are 4 % and 8 % smaller than those of TRIM.
- 400 keV H in Si: TRIM gives 45057 ± 75 Å and 2382 Å, and the SR table 44132 Å and 2351 Å (2 % and 1 % smaller). The TRIM distribution has a skewness of −8.9 and a kurtosis of 142, caused by a few ions that stop much closer to the surface after large-angle collisions. The Gaussian describes the peak but not this tail.

When the shape of the profile matters, a TRIM calculation (`RANGE.txt`, `RANGE_3D.txt`) should be used instead of Eq. (29).

Implementation: `lunasrim/web/core.py`, function `implant`.

## 8. Numerical reproduction of SRIM

The same formulas evaluated with a different arithmetic do not give the same output files. In a Monte Carlo calculation a difference in the last bit can change a decision, such as the number of Newton iterations in Eq. (21) or a path cut at a layer boundary, and from then on the trajectory differs. lunasrim therefore reproduces the arithmetic of the Visual Basic 5 runtime of SRIM-2013:

- Random numbers. TRIM uses L'Ecuyer's combined multiplicative congruential generator without the Bays-Durham shuffle of `ran2`, with Visual Basic's `CLng` (round half to even) in the integer division. The two seeds are the "Random Number Seed" and "Reminders" fields of TRIM.IN; zero seeds are replaced by 16381 and 12345 (`lunasrim/rng.py`).
- Variable types. TRIM's variables are Visual Basic Variants. The result type of each operation follows the OLE rules, so that, for example, an integer layer width divided by 100 gives a single-precision bin width, and the depth tallies accumulate in single precision. Recoil energies, positions and direction cosines are single-precision arrays.
- x87 arithmetic. The runtime runs the x87 FPU with control word 0x137F (64-bit precision). `Sqr` is therefore rounded to 64 bits and then to a double, which differs from a correctly rounded square root for many arguments just below 1. `Exp` and `^` are x87 instruction sequences of the Microsoft C runtime, and `Log`, `Atn`, `Sin` and `Cos` return the 80-bit result of one x87 instruction. lunasrim runs the same instruction sequences natively on x86-64 (`lunasrim/vbmath.py`, source in `tools/vbmath.S`).
- Number formatting. `Str$` and `CStr` use 15 significant digits (7 for single-precision values). `Format$` rounds a 17-digit conversion for scientific formats and a 15-digit conversion for fixed-point formats (`lunasrim/vbfmt.py`).

The x87 transcendental instructions are not correctly rounded and differ in the last bit between processor vendors. The reference outputs were made on an Intel processor, so lunasrim uses its native routines only when they reproduce the Intel values, which it checks when the module is imported. Elsewhere it uses Python's `math` module and an exact model of the 64-bit square root. [validation.md](validation.md) gives the consequences of this choice for the reference cases.

## References

1. J. F. Ziegler, J. P. Biersack and U. Littmark, *The Stopping and Range of Ions in Solids* (Pergamon, New York, 1985).
2. J. F. Ziegler, M. D. Ziegler and J. P. Biersack, SRIM – The stopping and range of ions in matter (2010), Nucl. Instr. Meth. B 268 (2010) 1818.
3. J. P. Biersack and L. G. Haggmark, A Monte Carlo computer program for the transport of energetic ions in amorphous targets, Nucl. Instr. Meth. 174 (1980) 257.
4. N. Bohr, On the decrease of velocity of swiftly moving electrified particles in passing through matter, Phil. Mag. 30 (1915) 581; N. Bohr, The penetration of atomic particles through matter, Mat. Fys. Medd. Dan. Vid. Selsk. 18, no. 8 (1948).
5. M. J. Norgett, M. T. Robinson and I. M. Torrens, A proposed method of calculating displacement dose rates, Nucl. Eng. Des. 33 (1975) 50.
6. G. H. Kinchin and R. S. Pease, The displacement of atoms in solids by radiation, Rep. Prog. Phys. 18 (1955) 1.
