.. This is a comment. Note how any initial comments are moved by
   transforms to after the document title, subtitle, and docinfo.

.. demo.rst from: http://docutils.sourceforge.net/docs/user/rst/demo.txt

.. |EXAMPLE| image:: _images/temp.png
   :width: 1em

**********************
User's Guide
**********************

..
  COMMENT: .. contents:: Table of Contents

Overview of LEAPR
=====================
LEAPR is an NJOY module that prepares the thermal neutron scattering law :math:`S(\alpha,\beta)`. The scattering law is used to calculate the scattering cross sections, where :math:`\alpha` and :math:`\beta` are unitless momentum and energy change, respectively.

LEAPR is able to calculate the coherent elastic, incoherent elastic, and incoherent inelastic contributions to the scattering law. The incoherent contributions are either approximated using the free gas approximation [INTRODUCED IN THIS SECIONT] or calculated using equations from [THIS SECTION]. Coherent elastic contributions are calculated using equations from [THIS SECTION]. 



What choices does a user have?
===============================
LEAPR has the ability to represent coherent and incoherent thermal neutron scattering, both of which have an elastic and an inelastic component. If you are not familiar with these distinctions, consider visiting the "Theory" section. 

Below is a flowchart illustrating options available for data processing in LEAPR. The options are split into four main chunks, corresponding to the coherent/incoherent and elastic/inelastic options. Note that the only processing step that is strictly necessary is the continuous processing in the beginning, which describes both incoherent elastic and incoherent inelastic. The different options will be briefly discussed below.


.. figure:: _images/LEAPR_flowchart.jpg
    :width: 100%
    :align: center

    Flowchart showing options available for users. 









