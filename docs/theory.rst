.. This is a comment. Note how any initial comments are moved by
   transforms to after the document title, subtitle, and docinfo.

.. demo.rst from: http://docutils.sourceforge.net/docs/user/rst/demo.txt

.. |EXAMPLE| image:: _images/temp.png
   :width: 1em

**********************
Theory
**********************

..
  COMMENT: .. contents:: Table of Contents

Thermal Neutrons
=====================

This project aims to describe the ways in which low energy neutrons (with energy on the order of 1 eV or less) interact with material. Accurately describing these interactions is crucial for adequate modeling of thermal nuclear systems. A neutron at room temperature has an energy of approximately 0.025 eV, meaning that its de Broglie wavelength is about 1 angstrom which is close to typical interatomic spacing in materials. This can complicate neutron-target interactions and the wave-like behavior of neutrons must be carefully accounted for.


Overview of Scattering
=========================

The probability of a neutron with initial energy and solid angle :math:`(E,\Omega)` scattering to have some final energy and solid angle :math:`(E',\Omega')` is described using the double differential scattering cross section :math:`\sigma(E\rightarrow E', \Omega\rightarrow\Omega')`. This cross section describes both elastic and inelastic scattering, both of which have a coherent and incoherent contribution.

.. image:: _images/scatteringBreakdown.jpeg

In elastic scattering, total kinetic energy (i.e. sum of neutron and target kinetic energy) is conserved, which is not the case for inelastic scattering. Inelastic scattering thus requires for some excitation of the target to occur, which accounts for the difference between initial and final kinetic energy. 


Elastic Scattering
---------------------
.. math::
  \sigma(E\rightarrow E,\Omega\rightarrow\Omega') = \sigma_{coh}(E\rightarrow E,\Omega\rightarrow\Omega') + \sigma_{inc}(E\rightarrow E',\Omega\rightarrow\Omega') 




Inelastic Scattering
-----------------------

For inelastic scattering kinetic energy is not conserved, meaning that some excitation (or de-excitation) occurred. For high energy neutrons, this excitation is a *nuclear* excitation, where the target nucleus is brought to some excited state. This requires a significant amount of energy, and thus nuclear inelastic scattering is a *threshold* reaction, as seen below.

.. figure:: _images/U238_xs.png
    :width: 80%
    :align: center

    Elastic and nuclear inelastic scattering cross sections for U-238 (from NNDC). Note that nuclear inelastic scattering is a threshold reaction that does not appreciable contribute until incoming neutrons have an incoming energy of about 0.1 MeV.


For thermal (low energy) neutrons, this excitation is typically a *molecular* or *lattice* excitation, where vibrational modes of a multi-atom system are excited. Molecular excitations can be induced by neutrons with energy on the order of 1 eV and do not exhibit the same extreme threshold behavior as does nuclear excitations. Thermal inelastic scattering is thus focused on molecular excitations. 

.. math::
  \sigma(E\rightarrow E',\Omega\rightarrow\Omega') = \sigma_{coh}(E\rightarrow E',\Omega\rightarrow\Omega') + \sigma_{inc}(E\rightarrow E,\Omega\rightarrow\Omega') 




Pair Distribution Function
===============================

The scattering kernel :math:`\sigma_s(E\rightarrow E',\Omega\rightarrow\Omega')` is typically separated into a coherent and an incoherent contribution (Note that in separating the coherent and incoherent contributions, one ignores spin-correlation effects. These effects are of little to no importance for most applications, except for instances like thermal scattering in liquid hydrogen, which has correlated spins (CITE PARKS). Such materials are considered apart from this simple coherent/incoherent discussion), both of which can be defined in terms of so-called "van Hove pair distribution function" :math:`G(\boldsymbol{r},t)`, which contains information regarding the scattering material. The pair distribution function is split into two terms :math:`G(\boldsymbol{r},t)=G_s(\boldsymbol{r},t)+G_d(\boldsymbol{r},t)` which represent a "self" term and a "distinct" term, respectively. In a classical system, :math:`G(\boldsymbol{r},t)` can be interpreted as the probability that an atom will be at location :math:`\boldsymbol{r}` at time :math:`t`, given that an atom existed at the origin at time :math:`t=0`.

The first term, :math:`G_s(\boldsymbol{r},t)`, represents the probability the particle originally at the origin would later exist at position :math:`\boldsymbol{r}`. The latter term :math:`G_d(\boldsymbol{r},t)` assumes that the two particles observed were not the same (CITE pairDist,bell-glasstone).

Using these definitions of the pair distribution functions, the coherent and incoherent scattering kernels for a homogeneous system consisting of bound scatterers of a single nuclide can be described as

.. math::
      \sigma_{coh}(E\rightarrow E',\Omega\rightarrow\Omega') = \frac{\sigma_{coh}}{4\pi} \sqrt{\frac{E'}{E}} \frac{1}{2\pi} \int dt\int d\boldsymbol{r}~\mathrm{e}^{i(\boldsymbol{\kappa\cdot r}-\epsilon t/\hbar)} G(\boldsymbol{r},t)\label{eq:pairDistInXS_coh}


.. math::
      \sigma_{inc}(E\rightarrow E',\Omega\rightarrow\Omega') = \frac{\sigma_{inc}}{4\pi} \sqrt{\frac{E'}{E}} \frac{1}{2\pi} \int dt\int d\boldsymbol{r}~\mathrm{e}^{i(\boldsymbol{\kappa\cdot r}-\epsilon t/\hbar)} G_s(\boldsymbol{r},t)\label{eq:pairDistInXS_inc}

where :math:`\hbar\boldsymbol{\kappa}` is the change in neutron momentum, :math:`\epsilon` is the change in energy, and :math:`\sigma_{coh}` and :math:`\sigma_{inc}` are the bound coherent and incoherent scattering cross sections, respectively (CITE bell-glasstone). These bound cross sections can be defined in terms of the first and second moments of the scattering length :math:`b`,

.. math:: \sigma_{coh}=4\pi\langle a\rangle^2 
.. math:: \sigma_{inc}=4\pi\left(\langle a^2\rangle-\langle a\rangle^2 \right)

where :math:`\langle...\rangle` denotes the average (CITE sturm1993dynamic). These definitions assume that the spins of adjacent nuclei are randomly oriented. For reactor purposes this is a very good approximation, since the spins of neighboring nuclei are typically uncorrelated except at very low temperatures.


Incoherent Scattering (Elastic and Inelastic)
==============================================
LEAPR describes incoherent scattering by starting with a continuous (solid-type) distribution and considering additional feautres (e.g. the existence of very sharp peaks or diffusive behavior) if necessary. The continuous calculation *must always be performed*, and additional features that can also be considered include the existence of discrete oscillators (Einstein oscillators), translative/diffusive behavior, cold Hydrogen/Deuterium, and intermolecular interference. 

The continuous calculation will now be discussed, followed by each of the aforementioned features.


Continuous Treatment 
-------------------------

To calculate the incoherent contribution to the scattering law, the following equations must be solved,

.. math::
    S_{n.sym}(\alpha, \beta)=\frac{1}{2 \pi} \int_{-\infty}^{\infty} \mathrm{e}^{i \beta t} \mathrm{e}^{-\gamma(t)} d t

.. math::
    \gamma(t)=\alpha\lambda_s -\alpha \int_{-\infty}^\infty P(\beta')~\mathrm{e}^{-\beta'/2}~\mathrm{e}^{-i\beta' t}~d\beta'

.. math:: 
  P(\beta)=\frac{\rho(\beta)}{2\beta\sinh(\beta/2)},

where :math:`\rho(\beta)` is the phonon frequency distribution, and :math:`\lambda_s` is the Debye-Waller coefficient. By Taylor expanding :math:`\gamma(t)`, the above can be simplified to 

.. math:: 
    S(\alpha,\beta) = \mathrm{e}^{-\alpha\lambda_s}\sum_{n=0}^\infty \frac{\alpha^n}{n!} W_n(\beta)
.. math:: 
    W_1(\beta) = P(\beta)~\mathrm{e}^{-\beta/2}
.. math::
    W_n(\beta) = \int_{-\infty}^\infty W_1(\beta')~W_{n-1}(\beta-\beta')~d\beta'.


.. note::
  **What approximations are made?**
  In describing the scattering behavior with a continuous, solid-type spectrum, a number of approximations are made, which are summarized below.

  +-------------------+---------------------------------------------------+-----------------------------------+
  | Approximation     | Description                                       | Comments                          |
  |                   |                                                   |                                   |
  +===================+===================================================+===================================+
  | | Gaussian        | | In the definition of :math:`S(\alpha,\beta)`,   | | See Parks [REFERENCE] for more  | 
  | | approximation   | | :math:`\mathrm{exp}\Big(-\alpha\gamma(t)        | | discussion. Not typically       |
  |                   |   +\alpha^2\gamma_2(t)+\dots\Big)`                | | considered a significant source |
  |                   | | is approximated as :math:`\mathrm{exp}          | | of error                        | 
  |                   |   \big(-\alpha\gamma(t)\big)`.                    |                                   | 
  +-------------------+---------------------------------------------------+-----------------------------------+
  | | Incoherent      | | While this is not an approximation made         | | This is typically valid practice|
  | | approximation   | | in the above equations, it is important         | | for materials with either low   |
  |                   | | to note that LEAPR uses the incoherent          | | coherent cross sections or      |
  |                   | | equations to desribe both coherent and          | | randomly-oriented crystallites. |
  |                   | | incoherent scattering.                          |                                   |
  +-------------------+---------------------------------------------------+-----------------------------------+
  | | Validity of     | | The continuous calculation requires an          | | Before using a published phonon |
  | | phonon DOS      | | an input phonon spectrum, which are             | | spectrum, ensure that it        |
  |                   | | specific to material composition,               | | corresponds to the correct      |
  |                   | | crystalline structure, and temperature.         | | material, structure, and temp.  |
  +-------------------+---------------------------------------------------+-----------------------------------+

  ================================= =============
   Approximation                    Description 
  ================================= =============
  Gaussian approximation             We 
  Representation of coherence        False
  Validity of phonon distribution    A word of caution is to be said regarding phonon spectra use. 
  ================================= =============



Continuous Treatment 
-------------------------

To calculate the incoherent contribution to the scattering law, the following equations must be solved,

.. math::
    S_{n.sym}(\alpha, \beta)=\frac{1}{2 \pi} \int_{-\infty}^{\infty} \mathrm{e}^{i \beta t} \mathrm{e}^{-\alpha\gamma(t)} d t

.. math::
    \gamma(t)=\lambda_s -\int_{-\infty}^\infty P(\beta')~\mathrm{e}^{-\beta'/2}~\mathrm{e}^{-i\beta' t}~d\beta'

.. math:: 
  P(\beta)=\frac{\rho(\beta)}{2\beta\sinh(\beta/2)}


where :math:`\rho(\beta)` is the phonon distribution (also called a vibrational frequency or phonon density of states, DOS ). The DOS is provided by the user on an equally-spaced :math:`\beta` grid, where :math:`\beta` corresponds to unitless energy transfer.

The definition of :math:`\gamma(t)` also makes use of the *Debye-Waller factor* :math:`\lambda_s`, which is defined as 

.. math:: 
  \lambda_s=\int_{-\infty}^\infty P(\beta')~\mathrm{e}^{-\beta'/2}~d\beta'

For a more involved discussion on the phonon frequency spectrum and the Debye-Waller factor, please see [SECTION ON STUFF].


To facilitate calculation of the scattering law, the latter :math:`\gamma` exponential is expanded as a Taylor series,

.. math:: 
  \begin{align*}
    \mathrm{e}^{-\alpha\gamma(t)} &= \mathrm{e}^{-\alpha\lambda_s} \mathrm{exp}\left[\alpha \int_{-\infty}^\infty P(\beta')~\mathrm{e}^{-\beta'/2}~\mathrm{e}^{-i\beta' t}~d\beta'\right]\\
                            &= \mathrm{e}^{-\alpha\lambda_s} \sum_{n=0}^\infty\frac{1}{n!}\left[\alpha \int_{-\infty}^\infty P(\beta')~\mathrm{e}^{-\beta'/2}~\mathrm{e}^{-i\beta' t}~d\beta'\right]^n
   \end{align*}


This brings the scattering law to

.. math::
   \begin{align*} 
    S_{n.sym}(\alpha, \beta)&=\frac{1}{2 \pi} \int_{-\infty}^{\infty} \mathrm{e}^{i \beta t} \mathrm{e}^{-\alpha\gamma(t)} d t\\
    &=\frac{1}{2 \pi} \int_{-\infty}^{\infty}  \mathrm{e}^{i \beta t} \left(\mathrm{e}^{-\alpha\lambda_s} \sum_{n=0}^\infty\frac{1}{n!}\left[\alpha \int_{-\infty}^\infty P(\beta')~\mathrm{e}^{-\beta'/2}~\mathrm{e}^{-i\beta' t}~d\beta'\right]^n\right) d t\\
    &=\mathrm{e}^{-\alpha\lambda_s} \sum_{n=0}^\infty \frac{\alpha^n}{n!} \frac{1}{2\pi}\int_{-\infty}^{\infty}  \mathrm{e}^{i \beta t} \left[\int_{-\infty}^\infty P(\beta')~\mathrm{e}^{-\beta'/2}~\mathrm{e}^{-i\beta' t}~d\beta'\right]^nd t.
   \end{align*} 

By defining

.. math::
    W_n(\beta)=\frac{1}{2\pi}\int_{-\infty}^{\infty}  \mathrm{e}^{i \beta t} \left[\int_{-\infty}^\infty P(\beta')~\mathrm{e}^{-\beta'/2}~\mathrm{e}^{-i\beta' t}~d\beta'\right]^n

the scattering law is simplified to 

.. math:: 
    S(\alpha,\beta) = \mathrm{e}^{-\alpha\lambda_s}\sum_{n=0}^\infty \frac{\alpha^n}{n!} W_n(\beta)

where we see that 

.. math::
  \begin{align*}
    W_0(\beta) &= \frac{1}{2\pi}\int_{-\infty}^\infty\mathrm{e}^{i\beta t}~dt = \delta(\beta)\\
    W_1(\beta) &= \int_{-\infty}^\infty P(\beta')\mathrm{e}^{-\beta'/2}\left[\frac{1}{2\pi}\int_{-\infty}^\infty\mathrm{e}^{i(\beta-\beta')t}~dt\right]~d\beta' = P(\beta)~\mathrm{e}^{-\beta/2}
  \end{align*}

and it can even be shown that generation of successive :math:`W_n(\beta)` terms are obtained by convolving the previous term with the first term,

.. math::
  \begin{align*}
    W_n(\beta) &= \int_{-\infty}^\infty W_1(\beta')~W_{n-1}(\beta-\beta')~d\beta'.
  \end{align*}


These are the equations that LEAPR solves for while representing the scattering system with a continuous, solid-type frequency distribution. 


To summarize, the equations solved by leapr are 

.. math:: 
    S(\alpha,\beta) = \mathrm{e}^{-\alpha\lambda_s}\sum_{n=0}^\infty \frac{\alpha^n}{n!} W_n(\beta)
.. math:: 
    W_1(\beta) = P(\beta)~\mathrm{e}^{-\beta/2}
.. math::
    W_n(\beta) = \int_{-\infty}^\infty W_1(\beta')~W_{n-1}(\beta-\beta')~d\beta'.

where the most important user-provided input is the phonon frequency distribution :math:`\rho(\beta)`. 




Discrete Oscillators
-------------------------
Note that :math:`P(\beta)` is an even function, meaning that the definition of the Debye-Waller coefficient can be restated using hyperbolic cosine, 

.. math:: 
 \begin{align*}
  \lambda_s&=\int_{-\infty}^\infty P(\beta')~\mathrm{e}^{-\beta'/2}~d\beta'\\
           &=\int_{0}^\infty P(\beta')~2\cosh(\beta'/2)~d\beta.
 \end{align*}

Furthermore, when 


.. math::
   \begin{align*} 
   S_{n.sym}(\alpha, \beta)
    &=\frac{1}{2 \pi} \int_{-\infty}^{\infty}  \mathrm{e}^{i \beta t} \left(\mathrm{e}^{-\alpha\lambda_s} \sum_{n=0}^\infty\frac{1}{n!}\left[\alpha \int_{-\infty}^\infty P(\beta')~\mathrm{e}^{-\beta'/2}~\mathrm{e}^{-i\beta' t}~d\beta'\right]^n\right) d t\\
    &=\mathrm{e}^{-\alpha\lambda_s} \sum_{n=0}^\infty \frac{\alpha^n}{n!} \frac{1}{2\pi}\int_{-\infty}^{\infty}  \mathrm{e}^{i \beta t} \left[\int_{-\infty}^\infty P(\beta')~\mathrm{e}^{-\beta'/2}~\mathrm{e}^{-i\beta' t}~d\beta'\right]^nd t.
   \end{align*} 
..
   &=\frac{1}{2 \pi} \int_{-\infty}^{\infty} \mathrm{e}^{i \beta t} \mathrm{e}^{-\gamma(t)} d t\\




Translational and Diffusive Behavior
--------------------------------------

Free Gas
-----------

Cold Hydrogen and Deuterium
---------------------------



Coherent Scattering (Elastic Only)
==============================================

Hexagonal Lattices
-------------------------

Hexagonal Close Packed
-------------------------

Face Centered Cubic
--------------------------------------

Body Centered Cubic
---------------------







