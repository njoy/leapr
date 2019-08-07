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


Elastic vs. Inelastic Scattering
---------------------------------
Elastic scattering means that the total kinetic energy of the system (neutron plus target) is the same before and after the scattering collision. This is contrasted with inelastic scattering, where kinetic energy is not conserved. This change in energy is due to some excitation (or de-excitation) occurred. 

----------------------------------------------------------------------------

**Does an elastically scattered neutron have the same incoming/outgoing energy?**
If a neutron scatters off of a target of similar size (e.g. a hydrogen nucleus), the neutron can lose nearly all its energy. If the target is significantly more massive than the neutron, however, conservation of momentum prevents the neutron from losing a large amount of its energy. 

When fast / resonance range neutrons scatter off of media, their incoming energy is so much larger than the molecular bonds of the scattering material, which allows the neutron to effectively "see" the target as just a free atom. Thermal neutrons, however, do not have enough incoming energy to allow them to ignore molecular and lattice bonds. They thus cannot scatter off of a free atom, but rather scatter off of an atom that is part of a much larger aggregate system. This can make it harder for the neutron to lose a large fraction of its energy in an elastic scattering event. 

The incoming and outgoing energies of scattered thermal neutrons are not the same, but will have a tendency to be much closer than those of higher energy neutrons. 

----------------------------------------------------------------------------

**Isn't inelastic scattering a threshold reaction?**
For high energy neutrons, the excitation associated with an inelastic collision is a *nuclear* excitation, where the target nucleus is brought to some excited state. This requires a significant amount of energy, and thus nuclear inelastic scattering is a *threshold* reaction, as seen below.

.. figure:: _images/U238_xs.png
    :width: 60%
    :align: center

    Elastic and nuclear inelastic scattering cross sections for U-238 (from NNDC). Note that nuclear inelastic scattering is a threshold reaction that does not appreciable contribute until incoming neutrons have an incoming energy of about 0.1 MeV.


For thermal (low energy) neutrons, inelastic scattering is caused by some *molecular* or *lattice* excitation, where vibrational modes of a multi-atom system are excited. Molecular excitations can be induced by neutrons with energy on the order of 1 eV and do not exhibit the same extreme threshold behavior as does nuclear excitations. Thermal inelastic scattering is thus focused on molecular excitations. The availability of vibrational modes that could be excited in some lattice system is described by the vibrational frequency spectrum / phonon density of states / phonon frequency distribution. 

.. .. math::
  \sigma(E\rightarrow E',\Omega\rightarrow\Omega') = \sigma_{coh}(E\rightarrow E',\Omega\rightarrow\Omega') + \sigma_{inc}(E\rightarrow E,\Omega\rightarrow\Omega') 




Coherent vs. Incoherent 
--------------------------

Neutron scattering is the interaction of wavefunctions, where the incoming neuton wave interacts with a target and creates a scattered spherical wave. This is simply due to the fact that a large incoming wave hitting a (relatively) small target will result in a spherical scattered wave. Additionally, large thermal neutron wavelength means that the neutron can exist atop multiple atoms at once, creating simultaneous scattering sites. 

When these scattered spherical waves, which originate from different scattering sites interfere, they can do so either coherently (meaning periodic constructive growth or destructive cancellation) or incoherently (meaning that no large-scale periodic growth or cancellation occurs). 

Incoherent scattering is significantly easier to model, and LEAPR has the ability to desribe both elastic and inelastic incoherent scattering (which correspond to total kinetic energy conservation and change, respectively). Coherent scattering is harder to quantify, and LEAPR currently has the ability to describe only elastic coherent scattering for selected materials. 




Incoherent Scattering (Elastic and Inelastic)
==============================================
LEAPR describes incoherent scattering by starting with a continuous (solid-type) distribution and considering additional feautres (e.g. the existence of very sharp peaks or diffusive behavior) if necessary. The continuous calculation *must always be performed*, and additional features that can also be considered include the existence of discrete oscillators (Einstein oscillators), translative/diffusive behavior, cold Hydrogen/Deuterium, and intermolecular interference. 

The continuous calculation will now be discussed, followed by each of the aforementioned features.


Continuous Treatment 
-------------------------

To calculate the incoherent contribution to the scattering law, the following equations must be solved,

.. math::
    S^{(s)}_{n.sym}(\alpha, \beta)=\frac{1}{2 \pi} \int_{-\infty}^{\infty} \mathrm{e}^{i \beta t} \mathrm{e}^{-\gamma(t)} d t

.. math::
    \gamma(t)=\alpha\lambda_s -\alpha \int_{-\infty}^\infty P(\beta')~\mathrm{e}^{-\beta'/2}~\mathrm{e}^{-i\beta' t}~d\beta'

.. math:: 
  P(\beta)=\frac{\rho(\beta)}{2\beta\sinh(\beta/2)},

where :math:`S^{(s)}_{n.sym}` is the non-symmetric scattering law for solids, :math:`\rho(\beta)` is the phonon frequency distribution, and :math:`\lambda_s` is the Debye-Waller coefficient. For a introductory discussion on the phonon frequency spectrum and the Debye-Waller coefficient, please see [SECTION]. 

By Taylor expanding :math:`\gamma(t)`, the above can be simplified to 

.. math:: 
    S^{(s)}_{n.sym}(\alpha,\beta) = \mathrm{e}^{-\alpha\lambda_s}\sum_{n=0}^\infty \frac{\alpha^n}{n!} W_n(\beta)

where

.. math:: 
    W_1(\beta) = P(\beta)~\mathrm{e}^{-\beta/2}\qquad\mbox{and}\qquad W_n(\beta) = \int_{-\infty}^\infty W_1(\beta')~W_{n-1}(\beta-\beta')~d\beta'.


.. warning::
  **What approximations are made?**
  In describing the scattering behavior with a continuous, solid-type spectrum, a number of approximations are made, which are summarized below.

  +------------------+-------------------------------------------------+-----------------------------------+
  | Approximation    | Description                                     | Comments                          |
  |                  |                                                 |                                   |
  +==================+=================================================+===================================+
  | | Gaussian       | | In the definition of :math:`S(\alpha,\beta)`, | | See Parks [REFERENCE] for more  | 
  | | approximation  | | :math:`\mathrm{exp}\Big(-\alpha\gamma(t)      | | discussion. Not typically       |
  |                  |   +\alpha^2\gamma_2(t)+\dots\Big)`              | | considered a significant source |
  |                  | | is approximated as :math:`\mathrm{exp}        | | of error                        | 
  |                  |   \big(-\alpha\gamma(t)\big)`.                  |                                   | 
  +------------------+-------------------------------------------------+-----------------------------------+
  | | Incoherent     | | While this is not an approximation made       | | This is typically valid practice|
  | | approximation  | | in the above equations, it is important       | | for materials with either low   |
  |                  | | to note that LEAPR uses the incoherent        | | coherent cross sections or      |
  |                  | | equations to desribe both coherent and        | | randomly-oriented crystallites. |
  |                  | | incoherent scattering.                        |                                   |
  +------------------+-------------------------------------------------+-----------------------------------+
  | | Validity of    | | The continuous calculation requires an        | | Before using a published phonon |
  | | phonon DOS     | | an input phonon spectrum, which are           | | spectrum, ensure that it        |
  |                  | | specific to material composition,             | | corresponds to the correct      |
  |                  | | crystalline structure, and temperature.       | | material, structure, and temp.  |
  +------------------+-------------------------------------------------+-----------------------------------+
  | | Randomly       | | In order to separate scattering into          | | Spin-correlation is important   |
  | | oriented spins | | coherent and incoherent components, we        | | while considering materials like|
  |                  | | spin-correlation effects are ignored.         | | liquid hydrogen/deuterium, which|
  |                  | |                                               | | are considered separately later.|
  +------------------+-------------------------------------------------+-----------------------------------+
  | | Edge effects   | | Edge effects are not considered, meaning      | | Ignoring edge effects is a      |
  | | are ignored    | | that the scatterer is considered to be        | | standard approximations that is | 
  |                  | | infinitely large.                             | | is not typically considered to  |
  |                  |                                                 | | cause significant error.        |
  +------------------+-------------------------------------------------+-----------------------------------+

.. note::
  **Want more information?**

  +-------------------+---------------------------------------------+-----------------------------------+
  | Topic             | Internal resources                          | External resources                |
  +===================+=============================================+===================================+
  | | Phonon frequency| | please see [  ]                           |                                   |
  | | spectrum theory | |                                           |                                   |
  +-------------------+---------------------------------------------+-----------------------------------+
  | | Phonon frequency| | please see [  ]                           | | Materials project               |
  | | spectrum        |                                             | | can create your own             |
  | | availability    |                                             |                                   |
  +-------------------+---------------------------------------------+-----------------------------------+
  | | Derivation of   | | Please see [  ] this will show how we     |                                   |
  | | Eq. [   ]       | | got to the phonon expansion               |                                   |
  +-------------------+---------------------------------------------+-----------------------------------+
  | | Debye-Waller    | | please see [  ]                           |                                   |
  | | coefficient     | |                                           |                                   |
  +-------------------+---------------------------------------------+-----------------------------------+
  | |                 |                                             |                                   |
  +-------------------+---------------------------------------------+-----------------------------------+





Discrete Oscillators
-------------------------
The blue region in the below figure shows the vibrational frequency spectrum for hydrogen bound in water [inspired by [CITE DAMIAN]]. Note the two prominent peaks near 0.20 eV and 0.42 eV. If a user wanted to process scattering data for this material, they could provide this full spectrum to LEAPR and have it run the full continuous calculation.

.. figure:: _images/waterPhononDOS_hatch.png
    :width: 90%
    :align: center

    The vibrational frequency spectrum for H bound in water is shown above. 

Alternatively, LEAPR allows for users to represent these higher energy peaks as discrete oscillators, also known as "Einstein oscillators". These oscilltors are represented as weighted Dirac-:math:`\delta` functions in the frequency distribution, which brings the blue distribution in the above figure to become the red distribution. The lower energy, continuous distribution is still the same, but the two higher energy peaks are replaced with weighted :math:`\delta` functions (the weighting is not represented in the above figure).

As can be seen in the figure above, reducing the peaks to simple oscillators eliminates peak resolution and is **only recommended for validating and replicating existing data**. 

The scattering law contribution from a discrete oscillator is

.. math:: 
  S^{(i)}_{n.sym}(\alpha,\beta)=\mathrm{e}^{-\alpha\lambda_i}\sum_{n=-\infty}^\infty\delta(\beta-n\beta_i)~I_n\left[\frac{\alpha\omega_i}{\beta_i\sinh(\beta_i/2)}\right]~\mathrm{e}^{-n\beta_i/2}

which is a direct simplification of the scattering law from the continuous case (defined in Eq. [  ]).

To process the scattering law for a material described by discrete oscillators, the discrete ocsillator contribution :math:`S^{(i)}_{n.sym}(\alpha,\beta)` is calcululated for each :math:`i^{th}` oscillator. These individual contributions are convolved with the solid-type contribution :math:`S_{n.sym}^{(s)}(\alpha,\beta)` which, in the figure above corresponds with the lower-energy part of the red distribution.

.. warning::
  **What approximations are made?**
  The discrete oscillator formulation is a simplification of the continuous treatment, and thus adopts those along with additional approximations. Only the additional approximations are presented here. 

  +-------------------+----------------------------------------------+-----------------------------------+
  | Approximation     | Description                                  | Comments                          |
  |                   |                                              |                                   |
  +===================+==============================================+===================================+
  | | Einstein        | | The discrete oscillator approximation      | | See Parks [REFERENCE] for more  | 
  | | crystal approx. | | is the analytic solution for the           | | discussion. This has historical |
  |                   | | scattering law, when considering a         | | significance but is not         |
  |                   | | perfect cubic structure of atoms that      | | recommended for modern problems | 
  |                   | | all vibrate with the same frequency.       |                                   | 
  +-------------------+----------------------------------------------+-----------------------------------+

.. note::
  **Want more information?**

  +-----------------------+---------------------------------------------+-----------------------------------+
  | Topic                 | Internal resources                          | External resources                |
  +=======================+=============================================+===================================+
  | | Equivalence between | | please see [  ]                           |                                   |
  | | discrete oscillator | |                                           |                                   |
  | | and continuous      | |                                           |                                   |
  | | treatment           | |                                           |                                   | 
  +-----------------------+---------------------------------------------+-----------------------------------+
  | | Experimental support| | Please see [  ] this will show how we     |                                   |
  | | for validity of the | | got to the phonon expansion               |                                   |
  | | discrete oscillator | |                                           |                                   |
  | | treatment           | |                                           |                                   |
  +-----------------------+---------------------------------------------+-----------------------------------+





Translational and Diffusive Behavior
--------------------------------------
Thermal neutron scattering off of liquids can be described by solving a solid-type spectrum that is combined with a diffusive term. A popular diffusive model is the "Effective Width Model", which is defined as 

.. math:: 
  S_{n.sym}^{(t)}(\alpha,\beta) = \frac{2c\omega_t\alpha}{\pi}~\mathrm{exp}\left({2c^2\omega_t\alpha-\beta/2}\right)\sqrt{\frac{c^2+0.25}{\beta^2+4c^2\omega_t^2\alpha^2}}\mathrm{K}_1\left[\sqrt{c^2+0.25}\sqrt{\beta^2+4c^2\omega_t^2\alpha^2}\right]

with a corresponding frequency spectrum

.. math::
  \rho(\beta)=\omega_t\frac{4c}{\pi\beta}\sqrt{c^2+0.25}~\sinh(\beta/2)~\mbox{K}_1\left[\sqrt{c^2+0.25}~\beta\right]

where :math:`K_n(x)` is the modified Bessel function of the second kind with order :math:`n`.


An alternative to the effective width model is the free gas model, which is defined as 

.. math:: 
  S^{(f)}_{n.sym}(\alpha,\beta) = \frac{1}{\sqrt{4\pi\omega_t\alpha}}~\mathrm{exp}\left[-\frac{(\omega_t\alpha+\beta)^2}{4\omega_t\alpha}\right]


So to model a diffusive material, the solid-type solution obtained from the vibrational spectrum is convolved with some appropriate translative model (i.e. effective width model or free gas model). 




Coherent Scattering (Elastic Only)
==============================================

Coherent scattering is when periodic constructive growth or destructive cancellation of the scattered waves occur. This is a difficult phenomena to model, and thus LEAPR is currently limited to describing elastic coherent scattering for the following materials:
  +-----------------+------------------------------+
  | Materials       | Crystalline Structure        |
  +=================+==============================+
  | Graphite        | Hexagonal                    |
  +-----------------+------------------------------+
  | Beryllium Metal | Hexagonal Close-Packed (HCP) |
  +-----------------+------------------------------+
  | Beryllium Oxide | Inter-penetrating HCP        |
  +-----------------+------------------------------+
  | Aluminum        | Face-Centered Cubic (FCC)    |
  +-----------------+------------------------------+
  | Lead            | Face-Centered Cubic (FCC)    |
  +-----------------+------------------------------+
  | Iron            | Body-Centered Cubic (BCC)    |
  +-----------------+------------------------------+



The differential coherent scattering cross section is
.. math:: 
  \sigma_{coh}(E,\mu)=\frac{\sigma_c}{E}\sum_{E_i<E}f_i~\mathrm{e}^{-4W~E_i}~\delta(\mu-\mu_i)

where :math:`W` is the effevtive Debye-Waller coefficient, :math:`\sigma_c` is the bound coherent scattering cross section. :math:`E_i` are Bragg Edges, defined in term


.. math::
  E_i = \frac{\hbar^2\tau_i^2}{8m}

:math:`f_i` are defined in terms of the crystallographic structure factors :math:`F`.



  

Hexagonal Lattices
-------------------------
LEAPR's treatment of hexagonal lattices is heavily influenced from the HEXSCAT code. Summary of the theory will be presented here.

Hexagonal Close Packed
-------------------------

Face Centered Cubic
--------------------------------------

Body Centered Cubic
---------------------



Special Cases and Misc. Functions
==============================================

Cold Hydrogen and Deuterium 
-----------------------------


Short-time collision approximation 
------------------------------------

Skold and Vineyard
-----------------------------




