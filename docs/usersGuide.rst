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

+--------------------+------------+----------------------------+--------------+---------+---------+
| Parameter Name     | Symbol     |  Description               | Restriction  | Default | Card    |
+====================+============+============================+==============+=========+=========+
| | Number of        | ``nalpha`` | | Length of                |              |         |         | 
|   :math:`\alpha`   |            |   :math:`\alpha` vector    |              |         | 7.a     | 
| | values           |            |   that                     |  :math:`>0`  |         |         |
|                    |            | | is later given           |              |         |         |
+--------------------+------------+----------------------------+--------------+---------+---------+
| | Number of        |  ``nbeta`` | | Length of                |              |         |         | 
|   :math:`\beta`    |            |   :math:`\alpha` vector    |              |         | 7.b     | 
| | values           |            |   that                     |  :math:`>0`  |         |         |
|                    |            | | is later given           |              |         |         |
+--------------------+------------+----------------------------+--------------+---------+---------+
| | :math:`\alpha,   | ``lat``    | | If ``lat`` is set to     | 0 or 1       |         |         | 
|   \beta` scaling   |            |   1 all :math:`\alpha`     |              |         |         |
| | flag             |            | | and :math:`\beta`        |              |  0      | 7.c     | 
|                    |            |   values are scaled        |              |         |         |
|                    |            | | by                       |              |         |         |
|                    |            |   0.0253/:math:`k_bT`      |              |         |         |
+--------------------+------------+----------------------------+--------------+---------+---------+
| | :math:`\alpha`   | ``alpha``  | | :math:`\alpha`           | :math:`\geq  |         |         |
|   values           |            |   values provided in       | 0`           |         | 8       |
|                    |            | | increasing order         |              |         |         |
+--------------------+------------+----------------------------+--------------+---------+---------+
| | :math:`\beta`    | ``beta``   | | :math:`\beta`            | :math:`\geq  |         |         |
|   values           |            |   values provided in       | 0`           |         | 9       |
|                    |            | | increasing order         |              |         |         |
+--------------------+------------+----------------------------+--------------+---------+---------+
| | Temperature      | ``temp``   | | Temperature in Kelvin    | :math:`>0`   |         | 10      |
+--------------------+------------+----------------------------+--------------+---------+---------+







Continuous Treatment
-----------------------
The continuous treatment takes in a vibrational frequency spectrum (also called a phonon distribution) and computes a scattering law :math:`S(\alpha,\beta)` via the *phonon expansion method*. This approach makes numerous assumptions, which are outlined in :ref:`Theory of Incoherent Scattering<theory_incoherent>`. The following input parameters are necessary to perform a continuous treatment calculation.


+--------------------+------------+----------------------------+--------------+---------+---------+
| Parameter Name     | Symbol     |  Description               | Restriction  | Default | Card    |
+====================+============+============================+==============+=========+=========+
| | Order of phonon  | ``nphon``  | | Number of terms used     |              |         | 3.c     |
| | expansion        |            | | in the phonon            |              | 100     |         |
|                    |            | | expansion sum            |  :math:`>0`  |         |         |
+--------------------+------------+----------------------------+--------------+---------+---------+
| | Number of        | ``nalpha`` | | Length of                |              |         |         | 
|   :math:`\alpha`   |            |   :math:`\alpha` vector    |              |         | 7.a     | 
| | values           |            |   that                     |  :math:`>0`  |         |         |
|                    |            | | is later given           |              |         |         |
+--------------------+------------+----------------------------+--------------+---------+---------+
| | Number of        |  ``nbeta`` | | Length of                |              |         |         | 
|   :math:`\beta`    |            |   :math:`\alpha` vector    |              |         | 7.b     | 
| | values           |            |   that                     |  :math:`>0`  |         |         |
|                    |            | | is later given           |              |         |         |
+--------------------+------------+----------------------------+--------------+---------+---------+
| | :math:`\alpha,   | ``lat``    | | If ``lat`` is set to     | 0 or 1       |         |         | 
|   \beta` scaling   |            |   1 all :math:`\alpha`     |              |         |         |
| | flag             |            | | and :math:`\beta`        |              |  0      | 7.c     | 
|                    |            |   values are scaled        |              |         |         |
|                    |            | | by                       |              |         |         |
|                    |            |   0.0253/:math:`k_bT`      |              |         |         |
+--------------------+------------+----------------------------+--------------+---------+---------+
| | :math:`\alpha`   | ``alpha``  | | :math:`\alpha`           | :math:`\geq  |         |         |
|   values           |            |   values provided in       | 0`           |         | 8       |
|                    |            | | increasing order         |              |         |         |
+--------------------+------------+----------------------------+--------------+---------+---------+
| | :math:`\beta`    | ``beta``   | | :math:`\beta`            | :math:`\geq  |         |         |
|   values           |            |   values provided in       | 0`           |         | 9       |
|                    |            | | increasing order         |              |         |         |
+--------------------+------------+----------------------------+--------------+---------+---------+
| | Temperature      | ``temp``   | | Temperature in Kelvin    | :math:`>0`   |         | 10      |
+--------------------+------------+----------------------------+--------------+---------+---------+
| | Phonon grid      | ``delta``  | | The phonon distribution  |              |         |         |
| | spacing          |            | | will be provided on a    |              |         | 11.a    |
|                    |            | | uniform energy grid that |              |         |         |
|                    |            | | starts at 0. This is the |              |         |         |
|                    |            | | energy grid spacing in eV|              |         |         |
+--------------------+------------+----------------------------+--------------+---------+---------+
| | Number of points | ``ni``     | | Number of values that    |              |         |         |
| | in phonon grid   |            | | will be provided in the  |              |         | 11.b    |
|                    |            | | phonon distribution      |  :math:`>0`  |         |         |
+--------------------+------------+----------------------------+--------------+---------+---------+
| | Phonon           | ``rho``    | | Phonon distribution,     | | All values |         |         |
| | distribution     |            |   given                    | | must be    |         | 12      |
|                    |            | | on an equally-spaced grid| | :math:`>0` |         |         |
|                    |            | | of length ``ni`` with    |              |         |         |
|                    |            | | spacing ``delta``        |              |         |         |
+--------------------+------------+----------------------------+--------------+---------+---------+
| | Normalization    |            | | Continuous,              |              |         |         |
| | for continuous   |            |   translational,           | | (0,1]      |         |         |
| | component        | ``tbeta``  | | and discrete treatments  | | (cannot be |         | 13.a    |
|                    |            | | all have weights, which  | | 0 but can  |         |         |
|                    |            | | must sum to 1. This is   | | be 1)      |         |         |
|                    |            | | the weighting for the    |              |         |         |
|                    |            | | continuous, solid-type   |              |         |         |
|                    |            | | spectrum                 |              |         |         |
+--------------------+------------+----------------------------+--------------+---------+---------+



                            




Translational/Diffusive
-----------------------

+--------------------+------------+----------------------------+--------------+---------+---------+
| Parameter Name     | Symbol     |  Description               | Restriction  | Default | Card    |
+====================+============+============================+==============+=========+=========+
| | Number of        | ``nalpha`` | | Length of                |              |         |         | 
|   :math:`\alpha`   |            |   :math:`\alpha` vector    |              |         | 7.a     | 
| | values           |            |   that                     |  :math:`>0`  |         |         |
|                    |            | | is later given           |              |         |         |
+--------------------+------------+----------------------------+--------------+---------+---------+
| | Number of        |  ``nbeta`` | | Length of                |              |         |         | 
|   :math:`\beta`    |            |   :math:`\alpha` vector    |              |         | 7.b     | 
| | values           |            |   that                     |  :math:`>0`  |         |         |
|                    |            | | is later given           |              |         |         |
+--------------------+------------+----------------------------+--------------+---------+---------+
| | :math:`\alpha,   | ``lat``    | | If ``lat`` is set to     | 0 or 1       |         |         | 
|   \beta` scaling   |            |   1 all :math:`\alpha`     |              |         |         |
| | flag             |            | | and :math:`\beta`        |              |  0      | 7.c     | 
|                    |            |   values are scaled        |              |         |         |
|                    |            | | by                       |              |         |         |
|                    |            |   0.0253/:math:`k_bT`      |              |         |         |
+--------------------+------------+----------------------------+--------------+---------+---------+
| | :math:`\alpha`   | ``alpha``  | | :math:`\alpha`           | :math:`\geq  |         |         |
|   values           |            |   values provided in       | 0`           |         | 8       |
|                    |            | | increasing order         |              |         |         |
+--------------------+------------+----------------------------+--------------+---------+---------+
| | :math:`\beta`    | ``beta``   | | :math:`\beta`            | :math:`\geq  |         |         |
|   values           |            |   values provided in       | 0`           |         | 9       |
|                    |            | | increasing order         |              |         |         |
+--------------------+------------+----------------------------+--------------+---------+---------+
| | Temperature      | ``temp``   | | Temperature in Kelvin    | :math:`>0`   |         | 10      |
+--------------------+------------+----------------------------+--------------+---------+---------+
| | Phonon grid      | ``delta``  | | The phonon distribution  |              |         |         |
| | spacing          |            | | will be provided on a    |              |         | 11.a    |
|                    |            | | uniform energy grid that |              |         |         |
|                    |            | | starts at 0. This is the |              |         |         |
|                    |            | | energy grid spacing in eV|              |         |         |
+--------------------+------------+----------------------------+--------------+---------+---------+
| | Number of points | ``ni``     | | Number of values that    |              |         |         |
| | in phonon grid   |            | | will be provided in the  |              |         | 11.b    |
|                    |            | | phonon distribution      |  :math:`>0`  |         |         |
+--------------------+------------+----------------------------+--------------+---------+---------+
| | Phonon           | ``rho``    | | Phonon distribution,     | | All values |         |         |
| | distribution     |            |   given                    | | must be    |         | 12      |
|                    |            | | on an equally-spaced grid| | :math:`>0` |         |         |
|                    |            | | of length ``ni`` with    |              |         |         |
|                    |            | | spacing ``delta``        |              |         |         |
+--------------------+------------+----------------------------+--------------+---------+---------+
| | Normalization    |            | | Continuous,              |              |         |         |
| | for continuous   |            |   translational,           | | (0,1]      |         |         |
| | component        | ``tbeta``  | | and discrete treatments  | | (cannot be |         | 13.a    |
|                    |            | | all have weights, which  | | 0 but can  |         |         |
|                    |            | | must sum to 1. This is   | | be 1)      |         |         |
|                    |            | | the weighting for the    |              |         |         |
|                    |            | | continuous, solid-type   |              |         |         |
|                    |            | | spectrum                 |              |         |         |
+--------------------+------------+----------------------------+--------------+---------+---------+




Discrete Oscillators
-----------------------

Coherent Elastic 
==================

Coherent Inelastic
===================





