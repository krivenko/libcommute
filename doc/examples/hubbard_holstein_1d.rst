.. _hubbard_holstein_1d:

Partial diagonalization of a one-dimensional Hubbard-Holstein model
===================================================================

The model describes behavior of strongly correlated electrons coupled to
localized phonons [PBK74]_, [TC03]_. The Hamiltonian of the model defined on
a chain with :math:`L` sites reads

.. math::

    \hat H = -t \sum_{\sigma}\sum_{i=0}^{L-2}
        (c^\dagger_{i,\sigma} c_{i+1,\sigma} + h.c.)
    -\mu \sum_{i=0}^{L-1} n_{i,\sigma}
    +U \sum_{i=0}^{L-1} n_{i,\uparrow} n_{i,\downarrow} +\\+
    \omega \sum_{i=0}^{L-1} a^\dagger_i a_i +
    g \sum_{i=0}^{L-1} (n_{i,\uparrow} + n_{i,\downarrow})(a^\dagger_i + a_i),

where :math:`c^\dagger_{i,\sigma}`, :math:`c_{i,\sigma}` and
:math:`n_{i,\sigma}` are electron creation, annihilation and occupation number
operators respectively. Operators :math:`a^\dagger_i`, :math:`a_i` create and
destroy phonons localized at the chain site :math:`i`. :math:`t` is the hopping
constant of the electrons, :math:`\mu` is the chemical potential, and :math:`U`
is the on-site Coulomb repulsion constant. The localized phonons have frequency
:math:`\omega` and are coupled to the electrons with strength :math:`g`.
Each localized phonon is allowed to occupy only two states, :math:`|0\rangle`
and :math:`|1\rangle`. Diagonalization is performed within the sector with
:math:`N_{el} = 2` electrons by means of Eigen's ``SelfAdjointEigenSolver``.

.. literalinclude:: ../../examples/hubbard_holstein_1d.cpp
  :language: cpp
  :lines: 29-
  :linenos:

.. note::

  A Krylov subspace iteration eigensolver used instead of Eigen's
  ``SelfAdjointEigenSolver`` would require only repeated evaluation of
  :math:`\hat H|\psi\rangle`. Therefore, it would be unnecessary to store and
  fill the whole matrix ``Hmat`` at once.

.. [PBK74] "Low-temperature properties of the one-dimensional polaron band.
           I. Extreme-band-narrowing regime",
           G. Beni, P. Pincus and J. Kanamori,
           Phys. Rev. B 10, pp. 1896-1901 (1974),
           https://doi.org/10.1103/PhysRevB.10.1896

.. [TC03] "Possibility of a metallic phase in the
          charge-density-wave--spin-density-wave crossover region in the
          one-dimensional Hubbard-Holstein model at half filling",
          Y. Takada and A. Chatterjee,
          Phys. Rev. B 8, p. 081102 (2003),
          https://doi.org/10.1103/PhysRevB.67.081102
