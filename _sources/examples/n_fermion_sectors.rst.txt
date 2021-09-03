.. _n_fermion_sectors:

N-fermion sector views and exact diagonalization with Eigen 3
=============================================================

In this example we show how to partially diagonalize a large model of lattice
fermions. More specifically, we are going to use the single band Fermi-Hubbard
model on a 2D square lattice with periodic boundary conditions.

.. math::

  \hat H = -t\sum_{\langle i,j \rangle, \sigma}
            c^\dagger_{i,\sigma} c_{j,\sigma}
            -\mu \sum_{i, \sigma} n_{i,\sigma}
            + U \sum_i n_{i,\uparrow} n_{i,\downarrow}.

The number of atoms in the considered lattice (cluster) is set to 16, which
makes for an intractably large Hilbert space of dimension :math:`2^{32}`.
Since the Hamiltonian of the Hubbard model preserves total numbers of
spin-up and spin-down electrons, one can diagonalize the Hamiltonian within a
single :math:`N`-fermion sector or within a
:math:`(N_\uparrow, N_\downarrow)`-multisector.

In *libcommute*'s terms, an :math:`N`-fermion
:cpp:class:`sector <libcommute::n_fermion_sector_view>` is a subspace of a full
Hilbert space, which is spanned by all basis states with a fixed total
occupation of fermionic degrees of freedom (FDOF). Similarly, a
:cpp:class:`multisector <libcommute::n_fermion_multisector_view>` is
spanned by all basis states with a fixed occupation :math:`N_1` of a subset of
the FDOF :math:`\{S_1\}`, occupation :math:`N_2` of another subset
:math:`\{S_2\}` and so on. There can be any number of pairs
:math:`(\{S_i\}, N_i)` (sectors contributing to the multisector) as long as
all the subsets :math:`\{S_i\}` are disjoint.

In our example we consider a moderately sized sector with :math:`N = 2`
(:math:`\dim = {32 \choose 2} = 496`)
and a multisector with :math:`N_\uparrow = 1, N_\downarrow = 1`
(:math:`\dim = {16 \choose 1}{16 \choose 1} = 256`).

.. note::

  In general, the N-fermion (multi)sector views and functions do not require
  a purely fermionic system. The definitions of a sector and a
  multisector stand in presence of bosons, spins and any other non-fermionic
  degrees of freedom.

We use `Eigen <https://eigen.tuxfamily.org/>`_'s real vector container to store
elements of state vectors within the (multi)sector, and the
``SelfAdjointEigenSolver`` class to diagonalize Hamiltonian matrices.

.. literalinclude:: ../../examples/n_fermion_sectors.cpp
  :language: cpp
  :lines: 24-
  :linenos:
