.. _jaynes_cummings:

Jaynes-Cummings ladder
======================

The Jaynes-Cummings model was introduced in quantum optics to describe
interaction of a two-level atom with a single quantised mode of an optical
cavity [JC63]_. Phonons in the mode are created/destroyed by bosonic ladder
operators :math:`a^\dagger`/:math:`a`. Atom's Hilbert space is spanned by the
ground state :math:`|g\rangle` and the excited state :math:`|e\rangle`, and is
equivalent to that of a spin-half. The operators acting in this isospin space
are raising and lowering operators :math:`S_+ = |e\rangle\langle g|`,
:math:`S_- = |g\rangle\langle e|` and the atomic inversion operator
:math:`2S_z = |e\rangle\langle e| - |g\rangle\langle g|`.

Under the assumption that the photon frequency :math:`\omega_c` and the atomic
transition frequency :math:`\omega_a` satisfy
:math:`|\omega_c - \omega_a| \ll \omega_c + \omega_a` (small detuning
:math:`\delta = \omega_a - \omega_c`), one can write down the Jaynes-Cummings
Hamiltonian in the rotating wave approximation as follows,

.. math::

  \hat H = \hbar \omega_c a^\dagger a + \hbar \omega_a S_z +
           \hbar \frac{\Omega}{2}(a S_+ + a^\dagger S_-),

where :math:`\Omega` is the atom-field interaction strength. The Jaynes-Cummings
Hamiltonian commutes with operator :math:`\hat N = a^\dagger a + S_z`, which
means :math:`\hat H` takes on a block-diagonal structure. The ground state
:math:`|n = 0, g\rangle` is within its own one-dimensional block, while the
rest of the blocks are spanned by a pair :math:`|n-1, e\rangle`,
:math:`|n, g\rangle` each (:math:`n = \overline{1,\infty}`). Accordingly,
the ground state energy is :math:`E_g = -(1/2)\hbar\omega_a`, while all excited
states form a series of doublets :math:`E_\pm(n)` -- the so called
Jaynes-Cummings ladder.

.. math::

  E_\pm(n) = \hbar\omega_c \left(n - \frac{1}{2}\right) \pm
             \frac{\hbar}{2} \sqrt{\delta^2 + n\Omega^2}.

The respective eigenstates in each block (dressed states) are

.. math::

    |n, +\rangle &= \cos\left(\frac{\alpha_n}{2}\right) |n - 1, e\rangle +
                    \sin\left(\frac{\alpha_n}{2}\right) |n, g\rangle,\\
    |n, -\rangle &= \sin\left(\frac{\alpha_n}{2}\right) |n - 1, e\rangle -
                    \cos\left(\frac{\alpha_n}{2}\right) |n, g\rangle\\

with mixing angle

.. math::

    \alpha_n = \tan^{-1}\left(\frac{\Omega\sqrt{n}}{\delta}\right).

The program below verifies these analytical results numerically.

.. literalinclude:: ../../examples/jaynes_cummings.cpp
  :language: cpp
  :lines: 24-
  :linenos:

.. [JC63] "Comparison of quantum and semiclassical radiation theories
   with application to the beam maser",
   E. T. Jaynes and F.W. Cummings,
   Proc. IEEE, Vol. 51, issue 1, pp. 89-109 (1963),
   https://doi.org/10.1109/PROC.1963.1664
