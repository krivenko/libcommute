.. _spin1_heisenberg:

Spin-1 Heisenberg ring: Compressed storage of state vectors
===========================================================

In this example we explore the concept of a
:ref:`sparse Hilbert space <sparse_hilbert_space>`, and illustrate how using the
:cpp:class:`compressed_state_view <libcommute::compressed_state_view>` objects
can lead to a dramatic reduction in memory footprint. To this end, we consider
a ring of :math:`N` spins :math:`S=1` interacting via the antiferromagnetic
Heisenberg exchange coupling,

.. math::

  \hat H = \sum_{i=0}^{N-1} \mathbf{S}_i \cdot \mathbf{S}_{(i+1)\bmod N} =
           \sum_{\substack{i=0,\\j = (i+1)\bmod N}}^{N-1} S_{z,i} S_{z,j} +
              \frac{1}{2} \left[ S_{+,i} S_{-,j} + S_{-,i} S_{+,j} \right].

The state space of each spin is spanned by the :math:`S_z`-triplet
:math:`|0\rangle, |\pm1\rangle`, and the dimension of the model's Hilbert space
is thus :math:`3^N`. From the computing standpoint, however, it takes 2 bits to
specify one of the three basis states, with the bit combination ``0b11`` not
representing any physical state. According to the
:ref:`convention <hilbert_space>` adopted in *libcommute*, one then needs a bit
string of length :math:`2N` to encode a basis state of the entire system.
Although not all such strings correspond to physical basis states (most of the
strings are, in fact, unphysical for :math:`N\geq 3`), naively one would have to
allocate an array of size :math:`4^N` to store a state vector compatible with
:ref:`linear operators <loperator>` acting in this Hilbert space. This situation
is an example of a *sparse Hilbert space*, i.e. of a space whose basis state
indices do not form a continuous integer range. Sparse Hilbert spaces are
characterized by non-power-of-two dimensions.

It is a valid (and faster) option to allocate ``std::vector``
(``Eigen::Vector``, etc) objects of the size :math:`4^N` to store state
vectors, even though some of their elements will never be referenced. Another
option offered by *libcommute* is to allocate containers of the exponentially
smaller size :math:`3^N`, and to employ the
:cpp:class:`compressed_state_view <libcommute::compressed_state_view>` adaptor
that -- at the cost of a certain speed loss -- translates indices of
the physical basis states onto the continuous range :math:`[0;3^N-1]`.

.. literalinclude:: ../../examples/spin1_heisenberg.cpp
  :language: cpp
  :lines: 24-
  :linenos:
