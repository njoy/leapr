.. This is a comment. Note how any initial comments are moved by
   transforms to after the document title, subtitle, and docinfo.

.. demo.rst from: http://docutils.sourceforge.net/docs/user/rst/demo.txt

.. |EXAMPLE| image:: _images/temp.png
   :width: 1em

**********************
Coding Details
**********************

..
  COMMENT: .. contents:: Table of Contents

The main ``leapr`` function consists of a loop over temperatures that is performed for the primary scatterer, and then for a secondary scatterer if one exists. There are not many differences in how the primary and secondary scatterer are processed: they each have their own phonon distribution ``rho`` and phonon grid spacing ``delta``, and the :math:`\alpha` grid is scaled differently for the secondary scatterer. For the primary scatterer, the :math:`\alpha` and :math:`\beta` grids are scaled by :math:`0.0253/\mathrm{k_bT}` if ``lat = 1``, and are left alone otherwise. For the secondary scatterer, both grids are scaled by :math:`0.0253/\mathrm{k_bT}` if ``lat = 1``, but the :math:`\alpha` grid is *always* scaled by the mass ratio of the secondary scatterer to the primary scatterer (regardless of the value of ``lat``). 


Incoherent Scattering (Elastic and Inelastic)
==============================================

The incoherent scattering treatment is comprised of three functions (``contin``, ``trans``, and ``discre``), which calculate the continuous, the translational, and discrete oscillator contributions to the scattering law.

``contin`` is called first, and takes basic parameters (most imporantly, the :math:`\alpha` and :math:`\beta` grids and the phonon distribution) and returns the scattering law :math:`S(\alpha,\beta)`. If translational behavior is requested, the scattering law and translational-related inputs are given to the ``trans`` function, where it calculates the translational contribution and convolves it with the scattering law. This process is repeated in ``discre``, where the contribution of each oscillator is calculated and convolved with the existing scattering law. 


Continuous Treatment 
-------------------------
The continuous treatment occurs in the ``contin`` function, whose main purpose is to perform the phonon expansion calculation that is described in :ref:`Theory of Continuous Distribution Calculations<theory_incoherent_contin>`. ``contin`` uses the ``start`` function to obtain the first term in the convolution integral :math:`\mathcal{T}_1(\beta)` as well as the Debye-Waller factor :math:`\lambda_s` and the effective temperature :math:`T_{eff}`. In obtaining these values, ``start`` normlizes the :math:`\rho(\beta)` and :math:`P(\beta)` functions. Note that since the phonon distribution is known to behave parabolically :math:`\big(\rho(\beta)\sim c\beta^2\big)` for small values of :math:`\beta`, the first value of :math:`P(\beta)`, denoted as :math:`P_0`, is approximated as :math:`P_0\approx\rho_0\beta_1^{-2}`.

``contin`` then begins the phonon expansion loop, where successive terms :math:`\mathcal{T}_n(\beta)` are computed and convolved with the preceeding terms. The convolutions are performed using the ``convol`` function. Simplicity of ``convol`` function requires that the :math:`\mathcal{T}(\beta)` terms be on an equally spaced :math:`\beta` grid, which is why the input phonon distribution :math:`\rho(\epsilon)` is requested on an equally spaced grid. Trapezoidal integration is used in ``convol``.

At any one time, only :math:`\mathcal{T}_1,\mathcal{T}_{last},` and :math:`\mathcal{T}_{next}` need be stored. Note that due to the nature of the convolution integrals that are used to calculate the :math:`\mathcal{T}_n` functions, :math:`\mathcal{T}_{last}` and :math:`\mathcal{T}_{next}` will grow in size with each subsequent step. The size of nonzero values in these vectors is tracked by the ``nLast`` and ``nNext`` variables. Note that these :math:`\mathcal{T}_{last}` and :math:`\mathcal{T}_{next}` arrays contain only the :math:`-\beta` part of the non-symmetric :math:`\mathcal{T}` functions, since they tend to be more reasonably sized than those in the :math:`+\beta` side of the non-symmetric and those in the symmetric :math:`\mathcal{T}` functions. 

:math:`\mathcal{T}_1,\mathcal{T}_{last},`

.. math::
    S^{(s)}_{n.sym}(\alpha, \beta)=\frac{1}{2 \pi} \int_{-\infty}^{\infty} \mathrm{e}^{i \beta t} \mathrm{e}^{-\gamma(t)} d t

.. math::
    \gamma(t)=\alpha\lambda_s -\alpha \int_{-\infty}^\infty P(\beta')~\mathrm{e}^{-\beta'/2}~\mathrm{e}^{-i\beta' t}~d\beta'



Translational Behavior
--------------------------------------

The translational component is solved for in ``trans``, which is prepared to handle either a diffusive law or a free-gas law. Which law gets invoked depends on whether the diffusion coefficient provided is equal to zero. If it is equal to zero, a free gas scattering law is calculated and convolved with the existing :math:`S(\alpha,\beta)` that is output from ``contin``. Otherwise, the effective width model is solved and convolved with the scattering law. 

The first step done in ``contin`` is to compute the free gas or diffusive shape on an appropriate :math:`\beta` grid, which is done using ``getFreeGas`` or ``getDiffusion``, respectively. Then, ``sbfill`` is used to remap the crrent scattering law onto that same :math:`\beta` grid. The convolution integral is them performed. 

Consideration of translational behavior does not change the Debye-Waller coefficient, but it does change the effective temperature. 


Discrete Oscillators
-------------------------
The discrete oscillator treatment in ``discre`` begins by calling ``prepareParams``, which sets up vectors will be later used as arguments in the discrete oscillator equations. 


.. math:: 
  S^{(i)}_{n.sym}(\alpha,\beta)=\mathrm{e}^{-\alpha\lambda_i}\sum_{n=-\infty}^\infty\delta(\beta-n\beta_i)~I_n\left[\frac{\alpha\omega_i}{\beta_i\sinh(\beta_i/2)}\right]~\mathrm{e}^{-n\beta_i/2}


``discre`` loops through all :math:`\alpha` values as well as through all oscillators. Note that the sum has :math:`n` go from :math:`-\infty\rightarrow\infty`, and both the positive and negative terms are handled in the ``posNegTerms`` function. 


Coherent Scattering (Elastic)
==============================================


Coherent Scattering (Inelastic) Approximations
================================================

Skold and Vineyard
------------------------


Cold Hydrogen and Deuterium 
-------------------------------


