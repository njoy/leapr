.. This is a comment. Note how any initial comments are moved by
   transforms to after the document title, subtitle, and docinfo.

.. demo.rst from: http://docutils.sourceforge.net/docs/user/rst/demo.txt

.. |EXAMPLE| image:: _images/temp.png
   :width: 1em

**********************
QuickStart
**********************



NJOY input files are able to call multiple modules in sequence, which can be used together to perform complicated manipulations of the nuclear data. Broadly speaking, an NJOY input file will be structured as follows.

.. literalinclude:: exampleInputs/genericInput
  :language: html

Here, the first module name is specified, and the inputs are divided into "tapes" (i.e. lines) of inputs. Each module has its own number of possible tapes, and each tape has its own number of inputs. 

Simple :math:`\mbox{H}` in :math:`\mbox{H}_2\mbox{O}`
========================================================
This example serves as a *very* simplified input that would allow a user to prepare the thermal scattering law :math:`S(\alpha,\beta)` for describing hydrogen bound in water. These values are for instructional purposes only. Note that parentheses following the forward slash will indicate the tape number.

..
  COMMENT: .. contents:: Table of Contents
   :emphasize-lines: 1


|


.. literalinclude:: exampleInputs/simple_H_H2O
   :language: html
   :lineno-start: 0
   :lines: 1-2

The first line defines the module, which is LEAPR. Many modules are designed to follow another module, but LEAPR can be called alone. Typically, the first tape is used to specify input and output files. LEAPR does not read in any auxiliary files, so the first tape only specifies the output file, which in this case will be located in the (unless otherwise specified) bin directory as "tape24".

-------------------------------------------------------------------------------

.. literalinclude:: exampleInputs/simple_H_H2O
   :language: html
   :lineno-start: 0
   :lines: 1-3

Tape 2 contains the descriptive title that will be used in the final LEAPR output. 

-------------------------------------------------------------------------------

.. literalinclude:: exampleInputs/simple_H_H2O
   :language: html
   :lineno-start: 0
   :lines: 1-4

Here, we want to process the material at a single temperature, so the number of temperatures is set to 1. In addition to the output file (for us, tape24), NJOY also returns a summary output file, which is traditionally labeled "output". To control the amount of detail used when preparing this document, the second input of tape 3 is the print control option. This value can be (0,1,2) which indicate increasing level of detail. The final input value is the phonon expansion number, which indicates the number of terms that should be computed when approximating the sum in [EQUATION SECTION ETC]. 

.. math:: 
    S^{(s)}_{n.sym}(\alpha,\beta) = \mathrm{e}^{-\alpha\lambda_s}\sum_{n=0}^\infty \frac{\alpha^n}{n!} W_n(\beta)





-------------------------------------------------------------------------------

.. literalinclude:: exampleInputs/simple_H_H2O
   :language: html
   :lineno-start: 0
   :lines: 1-5














