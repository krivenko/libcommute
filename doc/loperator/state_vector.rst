.. _state_vector:

``StateVector`` concept
=======================

.. default-domain:: cpp

.. namespace:: libcommute

Let us say we want to make type :type:`SV` a *libcommute*-compatible state
vector type so that :ref:`linear operators <loperator>` can act on instances
of :type:`SV`. For this we have to make :type:`SV` model the ``StateVector``
concept.
In a nutshell, a ``StateVector`` type is a one-dimensional array allowing
integer indexing and implementing a certain interface. Elements of the array
do not have to be stored contiguously.
Acceptable index values must be at least 64-bit wide unsigned integers since
*libcommute* uses the following type for basis state indexing.

.. type:: sv_index_type = std::uint64_t

  *Defined in <libcommute/loperator/state_vector.hpp>*

  Type of basis state indices.

The table below shows the interface (a set of free functions and a
metafunction) that needs be implemented for an object :expr:`sv` of type
:type:`SV`.

*libcommute* provides a default implementation of the ``StateVector`` concept
(see *<libcommute/loperator/state_vector.hpp>*), which works without extra
effort for :type:`std::vector` and similar container types having methods
``size()`` and ``operator[]()``. Even if the default implementation works for
your type, it might be beneficial to override some of the functions from
the table to achieve optimal performance.

.. list-table::
  :header-rows: 1

  * - Function/metafunction
    - Description
    - Default implementation
  * - :expr:`element_type<SV>::type`
    - The type of the elements.
    - :expr:`decltype(std::declval<SV>()[0])`

  * - :expr:`get_element(sv, n)`
    - Return the :expr:`n`-th element of :expr:`sv`.
    - :expr:`sv[n]`

  * - :expr:`update_add_element(sv, n, value)`
    - Add a value of some type :type:`T` to the :expr:`n`-th element of
      :expr:`sv`.
    - :expr:`sv[n] += value` or :expr:`sv[n] = sv[n] + value`

      The compound-assignment from type :expr:`T` will be used
      whenever :expr:`sv`'s elements support it. Otherwise, the implementation
      will fall back to the simple addition.

  * - :expr:`set_zeros(sv)`
    - Fill :expr:`sv` with zeros.
    - :expr:`sv[n] = zero` in a for-loop.

      The zero value is created by
      :expr:`make_const(0)` as described in ":ref:`custom_scalar_type`".

  * - :expr:`zeros_like(sv)`
    - Return an object of the same type and size as :expr:`sv` but filled with
      zeros.
    - Creates a new object as :expr:`res = SV(sv.size())` and calls
      :expr:`set_zeros(res)`,

  * - :expr:`foreach(sv, f)`
    - Apply a function-like object :expr:`f` to all basis state index/non-zero
      element pairs :expr:`(n, a)` in :expr:`sv`.
    - In a for-loop, calls :expr:`f(n, a)` for all non-zero elements :expr:`a`
      as detected by :expr:`is_zero()` (see ":ref:`custom_scalar_type`").
