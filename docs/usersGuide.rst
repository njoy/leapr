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


What choices does a user have?
===============================
LEAPR has the ability to represent coherent and incoherent thermal neutron scattering, both of which have an elastic and an inelastic component. If you are not familiar with these distinctions, consider visiting :ref:`overview_types_of_scattering`. 

Below is a flowchart illustrating options available for data processing in LEAPR. The options are split into four main chunks, corresponding to the coherent/incoherent and elastic/inelastic options. Note that the only processing step that is strictly necessary is the continuous treatment in the beginning, which describes both incoherent elastic and incoherent inelastic. All other processes are optional, can be used to either supplement the incoherent treatment or to represent coherent scattering. 


.. figure:: _images/LEAPR_flowchart.jpg
    :width: 100%
    :align: center

.. Flowchart showing options available for users. 


In this discussion, thermal scattering treatment is broken into three main sections: incoherent (which includes elastic and inelastic), coherent elastic, and coherent inelastic. For those interested in the equations being solved and approximations made in each of these cases, please visit :ref:`theory`. This guide primarily serves to introduce the necessary inputs and where these inputs may be found. 


Incoherent (Elastic and Inelastic)
====================================
Every LEAPR run requires that a continuous, solid-type spectrum be processed. Translational/diffusive and discrete oscillator treatment can also be applied, if desired. These three options and their necessary inputs will be discussed below.


General Information
-----------------------

+--------------------+------------+----------------------------+--------------+---------+
| Parameter Name     | Symbol     |  Description               | Restriction  | Card    |
+====================+============+============================+==============+=========+
| | Number of        | ``nalpha`` | | Length of                |              |         | 
|   :math:`\alpha`   |            |   :math:`\alpha` vector    | :math:`>0`   | 7.a     | 
|   values           |            |                            |              |         |
+--------------------+------------+----------------------------+--------------+---------+
| | Number of        |  ``nbeta`` | | Length of                |              |         | 
|   :math:`\beta`    |            |   :math:`\alpha` vector    |              | 7.b     | 
|   values           |            |                            |  :math:`>0`  |         |
+--------------------+------------+----------------------------+--------------+---------+
| | :math:`\alpha`   | ``lat``    | | If ``lat`` is set to     | 0 or 1       |         | 
|   and :math:`\beta`|            |   1, all :math:`\alpha`    |              |         |
|   scaling flag     |            |   and :math:`\beta`        |              | 7.c     | 
|                    |            | | values are scaled        |              |         |
|                    |            |   by                       |              |         |
|                    |            |   0.0253/                  |              |         |
|                    |            |  :math:`\mathrm{k_bT}`     |              |         |
|                    |            | | *(default value of 0)*   |              |         | 
+--------------------+------------+----------------------------+--------------+---------+
| | :math:`\alpha`   | ``alpha``  | | :math:`\alpha` values    | :math:`\geq  |         |
|   values           |            |   given in increasing order| 0`           | 8       |
|                    |            | | A total of ``nalpha``    |              |         |
|                    |            |   values needed            |              |         |
+--------------------+------------+----------------------------+--------------+---------+
| | :math:`\beta`    | ``beta``   | | :math:`\beta` values     | :math:`\geq  |         |
|   values           |            |   given in increasing order| 0`           | 9       |
|                    |            | | A total of ``nbeta``     |              |         |
|                    |            |   values needed            |              |         |
+--------------------+------------+----------------------------+--------------+---------+
| | Temperature      | ``temp``   | | Temperature in Kelvin    | :math:`>0`   | 10      |
+--------------------+------------+----------------------------+--------------+---------+







Continuous Treatment
-----------------------
The continuous treatment takes in a vibrational frequency spectrum (also called a phonon distribution) and computes a scattering law :math:`S(\alpha,\beta)` via the *phonon expansion method*. This approach makes numerous assumptions, which are outlined in :ref:`Theory of Incoherent Scattering<theory_incoherent>`. The following input parameters are necessary to perform a continuous treatment calculation.


+--------------------+------------+-----------------------------------+---------------+---------+
| Parameter Name     | Symbol     |  Description                      | Restriction   | Card    |
+====================+============+===================================+===============+=========+
| | Order of phonon  | ``nphon``  | | Number of terms used in the     |               | 3.c     |
| | expansion        |            | | phonon expansion sum            |               |         |
|                    |            | | *(default value 100)*           |  :math:`>0`   |         |
+--------------------+------------+-----------------------------------+---------------+---------+
| | Phonon grid      | ``delta``  | | The phonon distribution will be |               |         |
| | spacing          |            | | provided on a uniform energy    |               | 11.a    |
|                    |            | | grid, starting at 0. This is the|  :math:`>0`   |         |
|                    |            | | energy spacing in eV            |               |         |
+--------------------+------------+-----------------------------------+---------------+---------+
| | Number of points | ``ni``     | | Number of values in the phonon  |  :math:`>0`   |         |
| | in phonon grid   |            | | distribution that will be given |               | 11.b    |
+--------------------+------------+-----------------------------------+---------------+---------+
| | Phonon           | ``rho``    | | Phonon distribution values,     | | All ``rho`` |         |
| | distribution     |            |   given                           |   values      | 12      |
|                    |            | | on an equally-spaced grid of    | | must be     |         |
|                    |            | | length ``ni`` with spacing      |   :math:`>0`  |         |
|                    |            |   ``delta``                       |               |         |
+--------------------+------------+-----------------------------------+---------------+---------+
| | Normalization    |            | | Continuous, translational, and  |               |         |
| | for continuous   |            | | discrete oscillators all have   | :math:`0<`    |         |
| | component        | ``tbeta``  |   weights                         | ``tbeta``     | 13.a    |
|                    |            | | which sum to 1. This is the     | :math:`\leq1` |         |
|                    |            |   weight                          |               |         |
|                    |            | | for the continuous spectrum     |               |         |
+--------------------+------------+-----------------------------------+---------------+---------+



                            




Translational/Diffusive
-------------------------

+--------------------+------------+-----------------------------------+--------------+---------+
| Parameter Name     | Symbol     |  Description                      | Restriction  | Card    |
+====================+============+===================================+==============+=========+
| | Diffusion        | ``c``      | | The translational term can be   |              |         |
| | constant         |            | | either a free-gas law           | :math:`\geq0`| 13.b    |
|                    |            |   (``c``:math:`=0.0`)             |              |         |
|                    |            | | or a                            |              |         |
|                    |            |   diffusive law ``c``             |              |         |
|                    |            |   (:math:`>0.0`)                  |              |         |
+--------------------+------------+-----------------------------------+--------------+---------+
| | Normalization    |            | | Continuous, translational, and  |              |         |
| | for translational|            | | discrete oscillators all have   | :math:`0     |         |
| | component        | ``twt``    |   weights                         | \leq`        | 13.c    |
|                    |            | | which sum to 1. This is the     | ``twt``      |         |
|                    |            |   weight                          | :math:`<1`   |         |
|                    |            | | for the translational spectrum  |              |         |
+--------------------+------------+-----------------------------------+--------------+---------+




Discrete Oscillators
-----------------------

+--------------------+------------+-----------------------------------------+--------------+---------+
| Parameter Name     | Symbol     |  Description                            | Restriction  | Card    |
+====================+============+=========================================+==============+=========+
| | Number of        | ``nd``     | | The number of oscillators to be       |              |         |
| | oscillators      |            | | convolved with the existing scattering| :math:`\geq0`| 14      |
|                    |            | | law. If ``nd``:math:`=0` then cards 15|              |         |
|                    |            |   and 16                                |              |         |
|                    |            | | will not be read                      |              |         |
+--------------------+------------+-----------------------------------------+--------------+---------+
| | Oscillator       | ``bdel``   | | Energy locations of the oscillators   |              |         |
| | energies         |            | | are given here, in eV. There must be  | :math:`\geq0`| 15      |
|                    |            | | ``nd`` values provided                |              |         |
+--------------------+------------+-----------------------------------------+--------------+---------+
| | Oscillator       | ``adel``   | | Continuous, translational, and        |              |         |
| | weights          |            | | discrete oscillators all have         | :math:`\geq0`| 16      |
|                    |            |   weights                               |              |         |
|                    |            | | which sum to 1. These are the         |              |         |
|                    |            | | weights of each oscillator.           |              |         |
|                    |            | | There must be ``nd`` values           |              |         |
|                    |            |   given                                 |              |         |
+--------------------+------------+-----------------------------------------+--------------+---------+






Coherent Elastic 
==================


+--------------------+------------+-----------------------------------------+--------------+---------+
| Parameter Name     | Symbol     |  Description                            | Restriction  | Card    |
+====================+============+=========================================+==============+=========+
| | Number of        | ``nd``     | | The number of oscillators to be       |              |         |
| | oscillators      |            | | convolved with the existing scattering| :math:`\geq0`| 14      |
|                    |            | | law. If ``nd``:math:`=0` then cards 15|              |         |
|                    |            |   and 16                                |              |         |
|                    |            | | will not be read                      |              |         |
+--------------------+------------+-----------------------------------------+--------------+---------+
| | Oscillator       | ``bdel``   | | Energy locations of the oscillators   |              |         |
| | energies         |            | | are given here, in eV. There must be  | :math:`\geq0`| 15      |
|                    |            | | ``nd`` values provided                |              |         |
+--------------------+------------+-----------------------------------------+--------------+---------+
| | Oscillator       | ``adel``   | | Continuous, translational, and        |              |         |
| | weights          |            | | discrete oscillators all have         | :math:`\geq0`| 16      |
|                    |            |   weights                               |              |         |
|                    |            | | which sum to 1. These are the         |              |         |
|                    |            | | weights of each oscillator.           |              |         |
|                    |            | | There must be ``nd`` values           |              |         |
|                    |            |   given                                 |              |         |
+--------------------+------------+-----------------------------------------+--------------+---------+











Coherent Inelastic
===================





