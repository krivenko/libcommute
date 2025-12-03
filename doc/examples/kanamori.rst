.. _hubbard_kanamori:

Kanamori interaction Hamiltonian
================================

The Kanamori multiorbital Hamiltonian describes Coulomb repulsion of
:math:`t_{2g}` electronic states as relevant to a transition-metal ion
in a cubic crystal field with an octahedral environment [GdMM2013]_.
In its rotationally invariant form, the Hamiltonian reads

.. math::

  \begin{multline}
  \hat H_K = U \sum_{m} n_{m\uparrow} n_{m\downarrow} +
             (U-2J)\sum_{m\neq m'} n_{m\uparrow} n_{m'\downarrow} +\\+
             (U-3J)\sum_{m<m',\sigma} n_{m\sigma} n_{m'\sigma}
             -J\sum_{m\neq m'} c^\dagger_{m\uparrow} c_{m\downarrow}
                               c^\dagger_{m'\uparrow} c_{m'\uparrow}
             +J\sum_{m\neq m'} c^\dagger_{m\uparrow} c^\dagger_{m\downarrow}
                               c_{m'\downarrow} c_{m'\uparrow}.
  \end{multline}

The orbital indices :math:`m`, :math:`m'` run over the :math:`t_{2g}` triplet
:math:`p_x, p_y, p_z`. :math:`U` and :math:`J` are Coulomb integrals
(independent parameters of the model). The rotationally invariant Hamiltonian
can be expressed in terms of a few integrals of motion (Eq. (7) of [GdMM2013]_),

.. math::

  \hat H_{t_{2g}} = (U-3J)\frac{\hat N (\hat N-1)}{2}
                  - 2J\mathbf{S}^2 - \frac{J}{2}\mathbf{L}^2
                  + \frac{5}{2}J \hat N,

where :math:`\hat N` is the total number of electrons on the shell,

.. math::

  \hat N = \sum_{m\sigma} n_{m\sigma},

:math:`\mathbf{S}` is the total spin vector,

.. math::

  \mathbf{S} = \frac{1}{2} \sum_m \sum_{\sigma\sigma'}
    c^\dagger_{m\sigma} \boldsymbol{\tau}_{\sigma\sigma'} c_{m\sigma'},

and :math:`\mathbf{L}` is the total orbital isospin,

.. math::

  L_m = i \sum_{mm'} \sum_{\sigma} \epsilon_{mm'm''}
    c^\dagger_{m'\sigma} c_{m''\sigma}.

The following program shows that :math:`\hat H_K = \hat H_{t_{2g}}` indeed.

.. literalinclude:: ../../examples/kanamori.cpp
  :language: cpp
  :lines: 24-
  :linenos:

.. [GdMM2013] "Strong Correlations from Hund's Coupling",
   A. Georges, L. de' Medici and J. Mravlje,
   Annu. Rev. Condens. Matter Phys. 2013. 4:137–78,
   https://doi.org/10.1146/annurev-conmatphys-020911-125045
