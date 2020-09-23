.. _factories:

.. default-domain:: cpp

Factory functions for creation/annihilation/spin operators
==========================================================

Factory functions create the most basic :ref:`expressions <expression>`:
Single creation/annihilation/number operators of fermions and bosons, as well as
components of spin operators. These can be further combined to build all sorts
of many-body Hamiltonians and observables. The expression index types are
automatically deduced from the types of arguments passed to the factories.

There are six sets of factory functions - each defined in a separate namespace -
differing in the scalar type of the created expression
(real, complex or arbitrary) and its index types being static/dynamic.

.. namespace:: libcommute

.. list-table::
  :header-rows: 1

  * - Namespace
    - Expression type
  * - :expr:`libcommute::static_indices`
    - :expr:`expression<ScalarType, IndexTypes...>`
  * - :expr:`libcommute::static_indices::real`
    - :expr:`expression<double, IndexTypes...>`
  * - :expr:`libcommute::static_indices::complex`
    - :expr:`expression<std::complex<double>, IndexTypes...>`
  * - :expr:`libcommute::dynamic_indices`
    - :expr:`expression<ScalarType, dynamic_indices::dyn_indices>`
  * - :expr:`libcommute::dynamic_indices::real`
    - :expr:`expression<double, dynamic_indices::dyn_indices>`
  * - :expr:`libcommute::dynamic_indices::complex`
    - :expr:`expression<std::complex<double>, dynamic_indices::dyn_indices>`

.. note:: Factory functions for spin component operators :math:`S_x` and
          :math:`S_y` are defined only for the complex scalar type.

.. _factories_static:

Statically typed indices
------------------------

*Defined in <libcommute/expression/factories.hpp>*

.. rubric:: Namespace :expr:`libcommute::static_indices`

.. function:: template<typename ScalarType, typename... IndexTypes> \
              expression<ScalarType, IndexTypes...> \
              static_indices::c_dag(IndexTypes&&... indices)

  Make a fermionic creation operator :math:`c^\dagger` with given indices and
  an arbitrary scalar type.

.. function:: template<typename ScalarType, typename... IndexTypes> \
              expression<ScalarType, IndexTypes...> \
              static_indices::c(IndexTypes&&... indices)

  Make a fermionic annihilation operator :math:`c` with given indices and
  an arbitrary scalar type.

.. function:: template<typename ScalarType, typename... IndexTypes> \
              expression<ScalarType, IndexTypes...> \
              static_indices::n(IndexTypes&&... indices)

  Make a fermionic number operator :math:`n = c^\dagger c` with given indices
  and an arbitrary scalar type.

.. function:: template<typename ScalarType, typename... IndexTypes> \
              expression<ScalarType, IndexTypes...> \
              static_indices::a_dag(IndexTypes&&... indices)

  Make a bosonic creation operator :math:`a^\dagger` with given indices and
  an arbitrary scalar type.

.. function:: template<typename ScalarType, typename... IndexTypes> \
              expression<ScalarType, IndexTypes...> \
              static_indices::a(IndexTypes&&... indices)

  Make a bosonic annihilation operator :math:`a` with given indices and
  an arbitrary scalar type.

.. function:: template<typename ScalarType, typename... IndexTypes> \
              expression<ScalarType, IndexTypes...> \
              static_indices::S_p(IndexTypes&&... indices)

  Make a spin :math:`S=1/2` raising operator :math:`S_+` with given indices and
  an arbitrary scalar type.

.. function:: template<typename ScalarType, typename... IndexTypes> \
              expression<ScalarType, IndexTypes...> \
              static_indices::S_m(IndexTypes&&... indices)

  Make a spin :math:`S=1/2` lowering operator :math:`S_-` with given indices and
  an arbitrary scalar type.

.. function:: template<typename ScalarType, typename... IndexTypes> \
              expression<ScalarType, IndexTypes...> \
              static_indices::S_z(IndexTypes&&... indices)

  Make a spin :math:`S=1/2` z-projection operator :math:`S_z` with given indices
  and an arbitrary scalar type.

.. function:: template<int Multiplicity, \
              typename ScalarType, typename... IndexTypes> \
              expression<ScalarType, IndexTypes...> \
              static_indices::S_p(IndexTypes&&... indices)

  Make a general spin raising operator :math:`S_+` with given indices and
  an arbitrary scalar type. The :expr:`Multiplicity` template parameter must
  equal :math:`2S+1`, where :math:`S` is the wanted spin.

.. function:: template<int Multiplicity, \
              typename ScalarType, typename... IndexTypes> \
              expression<ScalarType, IndexTypes...> \
              static_indices::S_m(IndexTypes&&... indices)

  Make a general spin lowering operator :math:`S_-` with given indices and
  an arbitrary scalar type. The :expr:`Multiplicity` template parameter must
  equal :math:`2S+1`, where :math:`S` is the wanted spin.

.. function:: template<int Multiplicity, \
              typename ScalarType, typename... IndexTypes> \
              expression<ScalarType, IndexTypes...> \
              static_indices::S_z(IndexTypes&&... indices)

  Make a general spin z-projection operator :math:`S_z` with given indices
  and an arbitrary scalar type. The :expr:`Multiplicity` template parameter must
  equal :math:`2S+1`, where :math:`S` is the wanted spin.

.. rubric:: Namespace :expr:`libcommute::static_indices::real`

.. function:: template<typename... IndexTypes> \
              expression<double, IndexTypes...> \
              static_indices::real::c_dag(IndexTypes&&... indices)
              template<typename... IndexTypes> \
              expression<double, IndexTypes...> \
              static_indices::real::c(IndexTypes&&... indices)
              template<typename... IndexTypes> \
              expression<double, IndexTypes...> \
              static_indices::real::n(IndexTypes&&... indices)
              template<typename... IndexTypes> \
              expression<double, IndexTypes...> \
              static_indices::real::a_dag(IndexTypes&&... indices)
              template<typename... IndexTypes> \
              expression<double, IndexTypes...> \
              static_indices::real::a(IndexTypes&&... indices)
              template<typename... IndexTypes> \
              expression<double, IndexTypes...> \
              static_indices::real::S_p(IndexTypes&&... indices)
              template<typename... IndexTypes> \
              expression<double, IndexTypes...> \
              static_indices::real::S_m(IndexTypes&&... indices)
              template<typename... IndexTypes> \
              expression<double, IndexTypes...> \
              static_indices::real::S_z(IndexTypes&&... indices)
              template<int Multiplicity, \
              typename... IndexTypes> \
              expression<double, IndexTypes...> \
              static_indices::real::S_p(IndexTypes&&... indices)
              template<int Multiplicity, \
              typename... IndexTypes> \
              expression<double, IndexTypes...> \
              static_indices::real::S_m(IndexTypes&&... indices)
              template<int Multiplicity, \
              typename... IndexTypes> \
              expression<double, IndexTypes...> \
              static_indices::real::S_z(IndexTypes&&... indices)

  Specializations of the factory functions from
  :expr:`libcommute::static_indices` for the real expressions
  (scalar type :expr:`double`).

.. rubric:: Namespace :expr:`libcommute::static_indices::complex`

.. function:: template<typename... IndexTypes> \
              expression<std::complex<double>, IndexTypes...> \
              static_indices::complex::c_dag(IndexTypes&&... indices)
              template<typename... IndexTypes> \
              expression<std::complex<double>, IndexTypes...> \
              static_indices::complex::c(IndexTypes&&... indices)
              template<typename... IndexTypes> \
              expression<std::complex<double>, IndexTypes...> \
              static_indices::complex::n(IndexTypes&&... indices)
              template<typename... IndexTypes> \
              expression<std::complex<double>, IndexTypes...> \
              static_indices::complex::a_dag(IndexTypes&&... indices)
              template<typename... IndexTypes> \
              expression<std::complex<double>, IndexTypes...> \
              static_indices::complex::a(IndexTypes&&... indices)
              template<typename... IndexTypes> \
              expression<std::complex<double>, IndexTypes...> \
              static_indices::complex::S_p(IndexTypes&&... indices)
              template<typename... IndexTypes> \
              expression<std::complex<double>, IndexTypes...> \
              static_indices::complex::S_m(IndexTypes&&... indices)
              template<typename... IndexTypes> \
              expression<std::complex<double>, IndexTypes...> \
              static_indices::complex::S_z(IndexTypes&&... indices)
              template<int Multiplicity, \
              typename... IndexTypes> \
              expression<std::complex<double>, IndexTypes...> \
              static_indices::complex::S_p(IndexTypes&&... indices)
              template<int Multiplicity, \
              typename... IndexTypes> \
              expression<std::complex<double>, IndexTypes...> \
              static_indices::complex::S_m(IndexTypes&&... indices)
              template<int Multiplicity, \
              typename... IndexTypes> \
              expression<std::complex<double>, IndexTypes...> \
              static_indices::complex::S_z(IndexTypes&&... indices)

  Specializations of the factory functions from
  :expr:`libcommute::static_indices` for the complex expressions
  (scalar type :expr:`std::complex<double>`).

.. function:: template<typename... IndexTypes> \
              expression<std::complex<double>, IndexTypes...> \
              static_indices::complex::S_x(IndexTypes&&... indices)

  Make a spin :math:`S=1/2` x-projection operator :math:`S_x` with given indices
  and the complex scalar type.

.. function:: template<typename... IndexTypes> \
              expression<std::complex<double>, IndexTypes...> \
              static_indices::complex::S_y(IndexTypes&&... indices)

  Make a spin :math:`S=1/2` y-projection operator :math:`S_y` with given indices
  and the complex scalar type.

.. function:: template<int Multiplicity, \
              typename... IndexTypes> \
              expression<std::complex<double>, IndexTypes...> \
              static_indices::complex::S_x(IndexTypes&&... indices)

  Make a general spin x-projection operator :math:`S_x` with given indices
  and the complex scalar type. The :expr:`Multiplicity` template parameter must
  equal :math:`2S+1`, where :math:`S` is the wanted spin.

.. function:: template<int Multiplicity, \
              typename... IndexTypes> \
              expression<std::complex<double>, IndexTypes...> \
              static_indices::complex::S_y(IndexTypes&&... indices)

  Make a general spin y-projection operator :math:`S_y` with given indices
  and the complex scalar type. The :expr:`Multiplicity` template parameter must
  equal :math:`2S+1`, where :math:`S` is the wanted spin.

.. note:: Passing a C string literal as an index argument to a factory function
          will result in an expression with the corresponding index being a
          :type:`std::string`.

.. _factories_dyn:

[C++17] Dynamically typed index sequences
-----------------------------------------

*Defined in <libcommute/expression/factories_dyn.hpp>*

.. rubric:: Namespace :expr:`libcommute::dynamic_indices`

.. function:: template<typename ScalarType, typename... IndexTypes> \
              expression<ScalarType, dyn_indices> \
              dynamic_indices::c_dag(IndexTypes&&... indices)

  Make a fermionic creation operator :math:`c^\dagger` with given dynamically
  typed indices and an arbitrary scalar type.

.. function:: template<typename ScalarType, typename... IndexTypes> \
              expression<ScalarType, dyn_indices> \
              dynamic_indices::c(IndexTypes&&... indices)

  Make a fermionic annihilation operator :math:`c` with given dynamically typed
  indices and an arbitrary scalar type.

.. function:: template<typename ScalarType, typename... IndexTypes> \
              expression<ScalarType, dyn_indices> \
              dynamic_indices::n(IndexTypes&&... indices)

  Make a fermionic number operator :math:`n = c^\dagger c` with given
  dynamically indices typed and an arbitrary scalar type.

.. function:: template<typename ScalarType, typename... IndexTypes> \
              expression<ScalarType, dyn_indices> \
              dynamic_indices::a_dag(IndexTypes&&... indices)

  Make a bosonic creation operator :math:`a^\dagger` with given dynamically
  typed indices and an arbitrary scalar type.

.. function:: template<typename ScalarType, typename... IndexTypes> \
              expression<ScalarType, dyn_indices> \
              dynamic_indices::a(IndexTypes&&... indices)

  Make a bosonic annihilation operator :math:`a` with given dynamically typed
  indices and an arbitrary scalar type.

.. function:: template<typename ScalarType, typename... IndexTypes> \
              expression<ScalarType, dyn_indices> \
              dynamic_indices::S_p(IndexTypes&&... indices)

  Make a spin :math:`S=1/2` raising operator :math:`S_+` with given dynamically
  typed indices and an arbitrary scalar type.

.. function:: template<typename ScalarType, typename... IndexTypes> \
              expression<ScalarType, dyn_indices> \
              dynamic_indices::S_m(IndexTypes&&... indices)

  Make a spin :math:`S=1/2` lowering operator :math:`S_-` with given dynamically
  typed indices and an arbitrary scalar type.

.. function:: template<typename ScalarType, typename... IndexTypes> \
              expression<ScalarType, dyn_indices> \
              dynamic_indices::S_z(IndexTypes&&... indices)

  Make a spin :math:`S=1/2` z-projection operator :math:`S_z` with given
  dynamically typed indices and an arbitrary scalar type.

.. function:: template<int Multiplicity, \
              typename ScalarType, typename... IndexTypes> \
              expression<ScalarType, dyn_indices> \
              dynamic_indices::S_p(IndexTypes&&... indices)

  Make a general spin raising operator :math:`S_+` with given dynamically typed
  indices and an arbitrary scalar type. The :expr:`Multiplicity` template
  parameter must equal :math:`2S+1`, where :math:`S` is the wanted spin.

.. function:: template<int Multiplicity, \
              typename ScalarType, typename... IndexTypes> \
              expression<ScalarType, dyn_indices> \
              dynamic_indices::S_m(IndexTypes&&... indices)

  Make a general spin lowering operator :math:`S_-` with given dynamically typed
  indices and an arbitrary scalar type. The :expr:`Multiplicity` template
  parameter must equal :math:`2S+1`, where :math:`S` is the wanted spin.

.. function:: template<int Multiplicity, \
              typename ScalarType, typename... IndexTypes> \
              expression<ScalarType, dyn_indices> \
              dynamic_indices::S_z(IndexTypes&&... indices)

  Make a general spin z-projection operator :math:`S_z` with given dynamically
  typed indices and an arbitrary scalar type. The :expr:`Multiplicity` template
  parameter must equal :math:`2S+1`, where :math:`S` is the wanted spin.

.. rubric:: Namespace :expr:`libcommute::dynamic_indices::real`

.. function:: template<typename... IndexTypes> \
              expression<double, dyn_indices> \
              dynamic_indices::real::c_dag(IndexTypes&&... indices)
              template<typename... IndexTypes> \
              expression<double, dyn_indices> \
              dynamic_indices::real::c(IndexTypes&&... indices)
              template<typename... IndexTypes> \
              expression<double, dyn_indices> \
              dynamic_indices::real::n(IndexTypes&&... indices)
              template<typename... IndexTypes> \
              expression<double, dyn_indices> \
              dynamic_indices::real::a_dag(IndexTypes&&... indices)
              template<typename... IndexTypes> \
              expression<double, dyn_indices> \
              dynamic_indices::real::a(IndexTypes&&... indices)
              template<typename... IndexTypes> \
              expression<double, dyn_indices> \
              dynamic_indices::real::S_p(IndexTypes&&... indices)
              template<typename... IndexTypes> \
              expression<double, dyn_indices> \
              dynamic_indices::real::S_m(IndexTypes&&... indices)
              template<typename... IndexTypes> \
              expression<double, dyn_indices> \
              dynamic_indices::real::S_z(IndexTypes&&... indices)
              template<int Multiplicity, \
              typename... IndexTypes> \
              expression<double, dyn_indices> \
              dynamic_indices::real::S_p(IndexTypes&&... indices)
              template<int Multiplicity, \
              typename... IndexTypes> \
              expression<double, dyn_indices> \
              dynamic_indices::real::S_m(IndexTypes&&... indices)
              template<int Multiplicity, \
              typename... IndexTypes> \
              expression<double, dyn_indices> \
              dynamic_indices::real::S_z(IndexTypes&&... indices)

  Specializations of the factory functions from
  :expr:`libcommute::dynamic_indices` for the real expressions (scalar type
  :expr:`double`).

.. rubric:: Namespace :expr:`libcommute::dynamic_indices::complex`

.. function:: template<typename... IndexTypes> \
              expression<std::complex<double>, dyn_indices> \
              dynamic_indices::complex::c_dag(IndexTypes&&... indices)
              template<typename... IndexTypes> \
              expression<std::complex<double>, dyn_indices> \
              dynamic_indices::complex::c(IndexTypes&&... indices)
              template<typename... IndexTypes> \
              expression<std::complex<double>, dyn_indices> \
              dynamic_indices::complex::n(IndexTypes&&... indices)
              template<typename... IndexTypes> \
              expression<std::complex<double>, dyn_indices> \
              dynamic_indices::complex::a_dag(IndexTypes&&... indices)
              template<typename... IndexTypes> \
              expression<std::complex<double>, dyn_indices> \
              dynamic_indices::complex::a(IndexTypes&&... indices)
              template<typename... IndexTypes> \
              expression<std::complex<double>, dyn_indices> \
              dynamic_indices::complex::S_p(IndexTypes&&... indices)
              template<typename... IndexTypes> \
              expression<std::complex<double>, dyn_indices> \
              dynamic_indices::complex::S_m(IndexTypes&&... indices)
              template<typename... IndexTypes> \
              expression<std::complex<double>, dyn_indices> \
              dynamic_indices::complex::S_z(IndexTypes&&... indices)
              template<int Multiplicity, \
              typename... IndexTypes> \
              expression<std::complex<double>, dyn_indices> \
              dynamic_indices::complex::S_p(IndexTypes&&... indices)
              template<int Multiplicity, \
              typename... IndexTypes> \
              expression<std::complex<double>, dyn_indices> \
              dynamic_indices::complex::S_m(IndexTypes&&... indices)
              template<int Multiplicity, \
              typename... IndexTypes> \
              expression<std::complex<double>, dyn_indices> \
              dynamic_indices::complex::S_z(IndexTypes&&... indices)

  Specializations of the factory functions from
  :expr:`libcommute::dynamic_indices` for the complex expressions (scalar type
  :expr:`std::complex<double>`).

.. function:: template<typename... IndexTypes> \
              expression<std::complex<double>, dyn_indices> \
              dynamic_indices::complex::S_x(IndexTypes&&... indices)

  Make a spin :math:`S=1/2` x-projection operator :math:`S_x` with given
  dynamically typed indices and the complex scalar type.

.. function:: template<typename... IndexTypes> \
              expression<std::complex<double>, dyn_indices> \
              dynamic_indices::complex::S_y(IndexTypes&&... indices)

  Make a spin :math:`S=1/2` y-projection operator :math:`S_y` with given
  dynamically typed indices and the complex scalar type.

.. function:: template<int Multiplicity, \
              typename... IndexTypes> \
              expression<std::complex<double>, dyn_indices> \
              dynamic_indices::complex::S_x(IndexTypes&&... indices)

  Make a general spin x-projection operator :math:`S_x` with given dynamically
  typed indices and the complex scalar type. The :expr:`Multiplicity` template
  parameter must equal :math:`2S+1`, where :math:`S` is the wanted spin.

.. function:: template<int Multiplicity, \
              typename... IndexTypes> \
              expression<std::complex<double>, dyn_indices> \
              dynamic_indices::complex::S_y(IndexTypes&&... indices)

  Make a general spin y-projection operator :math:`S_y` with given dynamically
  typed indices and the complex scalar type. The :expr:`Multiplicity` template
  parameter must equal :math:`2S+1`, where :math:`S` is the wanted spin.

