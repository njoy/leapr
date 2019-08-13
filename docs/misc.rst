


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




