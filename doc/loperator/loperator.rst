.. _loperator:

Linear operators
================

.. default-domain:: cpp

.. namespace:: libcommute

*libcommute*'s linear operators implement action of
:ref:`polynomial expressions <expression>` on
:ref:`state vectors <state_vector>` in
:ref:`finite-dimensional Hilbert spaces <hilbert_space>`. Linear operator
classes are core devices used to build exact diagonalization routines based on
*libcommute*.

The linear operators come in two flavors, simple (:class:`loperator`)
and parametric (:class:`parametric_loperator`).

.. _simple_loperator:

Simple linear operator
----------------------

Simple linear operators are instances of class template :class:`loperator`.

Let us assume that :expr:`expr` is a real expression involving only the
predefined algebras (fermions, bosons and spins), and :expr:`hs` is a Hilbert
space object compatible with :expr:`expr`. Then creating a simple linear
operator is as easy as

.. code-block:: cpp

  auto L = libcommute::make_loperator(expr, hs);

Furthermore, let us declare a couple of objects that model the
:ref:`state vectors <state_vector>` concept (in this particular case, they are
standard vectors)

.. code-block:: cpp

  std::vector<double> psi(hs.vec_size());
  std::vector<double> phi(hs.vec_size());

The following code fragment shows three (nearly) equivalent ways to act with our
operator on a state, :math:`|\phi\rangle = \hat L |\psi\rangle`.

.. code-block:: cpp

  // Function call syntax
  phi = L(psi);

  // Multiplication syntax
  phi = L * psi;

  // 'In-place' action syntax
  L(psi, phi);

Although all three forms are semantically equivalent, the last one is faster as
it eliminates the need for a temporary object to store
:math:`\hat L|\psi\rangle`.

.. note::

  An important note must be made about various algebra support in
  :class:`loperator`. Expressions can dynamically accommodate new algebras via
  inheritance from the polymorphic base :class:`generator`. Unlike them, linear
  operators represent action of a certain set of algebras that is fixed at
  compile time. This design decision allows to remove the virtual function call
  overhead from the hottest parts of codes, where linear operators are
  repeatedly applied to state vectors. Template parameter pack
  :type:`loperator::AlgebraIDs` (sorted list of integers) determines what
  algebras a given linear operator object supports. :func:`make_loperator`
  returns operators covering all three predefined algebras. For the sake of
  optimization, it may make sense to manually create :class:`loperator` objects
  with a more restricted set of IDs if some of the predefined algebras will
  appear in expressions.

.. class:: template<typename ScalarType, int... AlgebraIDs> loperator

  *Defined in <libcommute/loperator/loperator.hpp>*

  Linear operator representation of a polynomial expression acting in
  a finite-dimensional Hilbert space.

  :type:`ScalarType` - coefficient type of the corresponding
  polynomial expression.

  :type:`AlgebraIDs` - IDs of algebras supported by this object.

  .. rubric:: Constructor

  .. function:: template<typename... IndexTypes> \
                loperator(expression<scalar_type, IndexTypes...> const& expr, \
                          hilbert_space<IndexTypes...> const& hs)

    Construct from an expression and a Hilbert space.

  .. rubric:: Copy/move-constructors and assignments

  .. function:: loperator(loperator const&) = default
  .. function:: loperator(loperator&&) noexcept = default
  .. function:: loperator& operator=(loperator const&) = default
  .. function:: loperator& operator=(loperator&&) noexcept = default

  .. rubric:: Action on a state vector

  .. function:: template<typename StateVector> \
                StateVector operator()(StateVector const& psi) const
                template<typename StateVector> \
                StateVector operator*(StateVector const& psi) const

    Act on a state vector :expr:`psi` and return the result
    :math:`\hat L|\psi\rangle`.

  .. function:: template<typename SrcStateVector, typename DstStateVector> \
                void operator()(SrcStateVector && psi, \
                                DstStateVector && phi) const

    Act on a state vector :expr:`psi` and write the result into :expr:`phi`,
    :math:`|\phi\rangle = \hat L |\psi\rangle`. This method is faster than the
    previous two because the result is written directly into :expr:`phi` without
    making a temporary object.

.. function:: template<typename ScalarType, typename... IndexTypes> \
              loperator<ScalarType, fermion, boson, spin> \
              make_loperator(expression<ScalarType, IndexTypes...> const& expr,\
              hilbert_space<IndexTypes...> const& hs)

  *Defined in <libcommute/loperator/loperator.hpp>*

  A helper factory function that constructs an :class:`loperator` instance.
  This function is a more convenient equivalent of :class:`loperator`'s
  constructor.

.. _param_loperator:

Parametric linear operator
--------------------------

Parametric linear operators are similar to
the :ref:`simple ones <simple_loperator>` with one important distinction. They
represent polynomial expressions with coefficients that support function-like
call syntax. This generalization can be of relevance for expressions that depend
on various parameters (hence the name). A rather common example here is when
the Hamiltonian of a quantum system explicitly depends on time or some external
fields.

When a parametric operator acts on a state, it takes a list of extra arguments.
These extra arguments are passed to the callable coefficients of the
corresponding expression and their return values are used to eventually form the
output state vector. In other words, the extra arguments are effectively
'substituted' into the expression.

.. code-block:: cpp

  // Coefficients of 'expr' depend on one real parameter, i.e. its
  // ScalarType implements the call operator taking one 'double' argument
  // and returning a 'double'.

  // Create a parametric linear operator
  auto PL = libcommute::make_param_loperator(expr, hs);

  // Input and output state vectors
  std::vector<double> psi(hs.vec_size());
  std::vector<double> phi(hs.vec_size());

  //
  // Substitute 2.0 into the expression and act on 'psi'
  //

  // Function call syntax
  phi = PL(psi, 2.0);

  // 'In-place' action syntax
  PL(psi, phi, 2.0);

.. note:: Omitting the extra arguments will result in the coefficients being
          called without arguments.

There is also an option to optimize the substitution process by preallocating
storage for the coefficient return values.

.. code-block:: cpp

  std::vector<double> evaluated_coeffs;

  // Only the first call will cause memory allocation
  PL.act_and_store_coeffs(psi, phi, evaluated_coeffs, .0);
  PL.act_and_store_coeffs(psi, phi, evaluated_coeffs, 1.0);
  PL.act_and_store_coeffs(psi, phi, evaluated_coeffs, 2.0);
  PL.act_and_store_coeffs(psi, phi, evaluated_coeffs, 3.0);

Finally, one can permanently transform a parametric operator into the
non-parametric form.

.. code-block:: cpp

  // Substitute 4.0 into PL and save the result as a simple linear operator
  auto L = PL.at(4.0);

For an extensive example of :class:`parametric_loperator`'s use have a look at
":ref:`parametric_loperator`".

.. class:: template<typename ScalarType, int... AlgebraIDs> \
           parametric_loperator

  *Defined in <libcommute/loperator/loperator.hpp>*

  Linear operator representation of a polynomial expression acting in
  a finite-dimensional Hilbert space. This class supports parameter substitution
  upon acting on a state.

  .. rubric:: Constructor

  .. function:: template<typename... IndexTypes> \
                parametric_loperator( \
                  expression<scalar_type, IndexTypes...> const& expr, \
                  hilbert_space<IndexTypes...> const& hs)

    Construct from an expression and a Hilbert space.

  .. rubric:: Copy/move-constructors and assignments

  .. function:: parametric_loperator(parametric_loperator const&) = default
  .. function:: parametric_loperator(parametric_loperator&&) noexcept = default
  .. function:: parametric_loperator& operator=(parametric_loperator const&) \
                  = default
  .. function:: parametric_loperator& operator=(parametric_loperator&&) \
                  noexcept = default

  .. rubric:: Action on a state vector

  .. function:: template<typename StateVector, typename... CoeffArgs> \
                StateVector \
                operator()(StateVector const& psi, CoeffArgs&&... args) const

        Act on a state vector :expr:`psi` and return the result
        :math:`\hat L(\text{args}\ldots)|\psi\rangle`.

  .. function:: template<typename StateVector, typename... CoeffArgs> \
                void operator()(StateVector const& psi, \
                                StateVector & phi, \
                                CoeffArgs&&... args) const

    Act on a state vector :expr:`psi` and write the result into :expr:`phi`,
    :math:`|\phi\rangle = \hat L(\text{args}\ldots) |\psi\rangle`.
    This method is faster than the previous one because the result is written
    directly into :expr:`phi` without making a temporary object.

  .. function:: template<typename StateVector, typename... CoeffArgs> \
                void act_and_store_coeffs( \
                StateVector const& psi, \
                StateVector & phi, \
                std::vector<evaluated_coeff_t<CoeffArgs...>>& evaluated_coeffs,\
                CoeffArgs&&... args) const

  Similar to the previous method, but using the external vector
  :expr:`evaluated_coeffs` to store coefficient values after parameter
  substitution. When called multiple times, only the first invocation will
  resize :expr:`evaluated_coeffs` and allocate memory.

  .. function:: template<typename... CoeffArgs> \
                loperator<evaluated_coeff_t<CoeffArgs...>, AlgebraIDs...> \
                at(CoeffArgs&&... args) const

    Transform this parametric linear operator into the non-parametric form
    by substituting parameters :expr:`args` into it.

.. function:: template<typename ScalarType, typename... IndexTypes> \
              parametric_loperator<ScalarType, fermion, boson, spin> \
              make_param_loperator( \
                expression<ScalarType, IndexTypes...> const& expr, \
                hilbert_space<IndexTypes...> const& hs)

  *Defined in <libcommute/loperator/loperator.hpp>*

  A helper factory function that constructs a :class:`parametric_loperator`
  instance. This function is a more convenient equivalent of
  :class:`parametric_loperator`'s constructor.

