.. This is a comment. Note how any initial comments are moved by
   transforms to after the document title, subtitle, and docinfo.

.. demo.rst from: http://docutils.sourceforge.net/docs/user/rst/demo.txt

.. |EXAMPLE| image:: _images/temp.png
   :width: 1em

**********************
Overview
**********************

..
  COMMENT: .. contents:: Table of Contents

What is LEAPR?
=====================
LEAPR is one of 24 modules that comprise the NJOY nuclear data preaparation code. Two modules (LEAPR and THERMR) handle data corresponding the thermal (low energy) neutron scattering. It was originally written and is currently maintained at Los ALamos National Laboratory. For more information regarding the NJOY code, please visit https://www.njoy21.io/NJOY2016/.

.. LEAPR prepares the **scattering law** :math:`S(\alpha,\beta)`, and THERMR writes the scattering law in a convenient way for use in simulations, etc. The scattering law can be used to calculate the scattering cross sections, where :math:`\alpha` and :math:`\beta` are unitless momentum and energy change, respectively.




Thermal Neutrons
=====================

LEAPR aims to describe the ways in which low energy neutrons (with energy on the order of 1 eV or less) interact with material. Accurately describing these interactions is crucial for adequate modeling of thermal nuclear systems. A neutron at room temperature has an energy of approximately 0.025 eV, meaning that its de Broglie wavelength is about 1 angstrom which is close to typical interatomic spacing in materials. This can complicate neutron-target interactions, and thus describing thermal scattering must account for the wave-like behavior of neutrons. 


Types of Scattering
=========================

The probability of a neutron with initial energy and solid angle :math:`(E,\Omega)` scattering to have some final energy and solid angle :math:`(E',\Omega')` is described using the double differential scattering cross section :math:`\sigma(E\rightarrow E', \Omega\rightarrow\Omega')`. This cross section describes both elastic and inelastic scattering, both of which have a coherent and incoherent contribution.

.. image:: _images/scatteringBreakdown.jpeg

In elastic scattering, total kinetic energy (i.e. sum of neutron and target kinetic energy) is conserved, which is not the case for inelastic scattering. Inelastic scattering thus requires for some excitation of the target to occur, which accounts for the difference between initial and final kinetic energy. 


Elastic vs. Inelastic 
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

For reactor systems, incoherent scattering primarily dominates (which facilitates data preparation, as incoherent scattering is simpler to process). There are some instances, however, where neglecting coherent scattering could result in significant error. For a brief discussion detailing *when* coherent scattering is important, please see [SOURCE].

Overview of Objectives
--------------------------
As seen above, there are multiple types of thermal scattering that occur in a reactor system. LEAPR aims to prepare data describing these different types of interactions, namely the **scattering law**, :math:`S(\alpha,\beta)`. The scattering law is a function of dimensionless momentum and energy exchange (:math:`\alpha` and :math:`\beta`, respectively), 

.. math:: 
  \alpha = \frac{E'+E-2\mu\sqrt{EE'}}{Ak_bT}\qquad~\qquad\beta=\frac{E'-E}{k_bT}.

Once obtained, the scattering law can be used to calculate the double differential scattering cross section,

.. math::
  \frac{d^2\sigma}{dE'~d\Omega }=\frac{\sigma_b}{2k_bT}\sqrt{\frac{E'}{E}}S(\alpha,\beta)

[IS THIS IN THE INCOHERENT APPROXIMATION I THINK SO BUT I"M NOT SURE CHECK YOUR BOOK]

Thus, the goal of LEAPR is to calculate this scattering law for some user-provided :math:`\alpha,\beta` grid. Doing so requires many approximation that will be described in the coming sections. Incoherent scattering is significantly easier to describe, and LEAPR calculates this contribution to the scattering law by use of a user-provided vibrational frequency spectrum (sometimes known as frequency distribution, phonon distribution, phonon density of states, etc.). Cohernt scattering is significantly more difficult to describe, and so LEAPR's ability to calculate coherent scattering contributions is much more limited. Coherent elastic scattering capabilities are available for certain materials, and coherent inelastic scattering can be approximated if additional material data (namely, the static structure factor) is provided.

