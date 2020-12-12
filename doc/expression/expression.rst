.. _expression:

Polynomial expressions
======================

This page describes how *libcommute* represent polynomial expressions built
out of non-commuting quantum-mechanical operators and what features such
expressions support.

.. _expr_classes:

Class definitions
-----------------

.. default-domain:: cpp

.. namespace:: libcommute

.. class:: template<typename ScalarType, typename... IndexTypes> \
           expression

  *Defined in <libcommute/expression/expression.hpp>*

  A polynomial expression, the main building block of *libcommute*'s DSL.

  Expressions are finite ordered sums of monomials
  :math:`M^{(n)}_{i_1 i_2 \ldots i_n}` accompanied by their respective
  coefficients :math:`C_{i_1 i_2 \ldots i_n}`,

  .. math::

    E = C M^{(0)} + \sum_i C_i M^{(1)}_i + \sum_{ij} C_{ij} M^{(2)}_{ij} +
        \ldots

  The :ref:`monomials <monomial>` are in turn canonically ordered products of
  :ref:`algebra generators <generator>` (operators
  :math:`c^\dagger`/:math:`c`, :math:`a^\dagger`/:math:`a`, etc).

  :type:`ScalarType` is the type of coefficients
  :math:`C_{i_1 i_2 \ldots i_n}`. The most common choices of :type:`ScalarType`
  are :expr:`double` for expressions with real coefficients and
  :expr:`std::complex<double>` for the complex expressions. There are two
  convenience type aliases in the :expr:`libcommute::static_indices` namespace
  for these particular scalar types, :type:`static_indices::expr_real` and
  :type:`static_indices::expr_complex`.
  Use of custom scalar types -- subject to some
  requirements -- is also allowed. See Section :ref:`custom_scalar_type`
  for more details.

  Indices :math:`i`, :math:`j`, etc in the definition above are compound
  indices with a fixed number of components. Types of the components are
  given by the template parameter pack :type:`IndexTypes`. **All index types
  must be less-comparable and form strictly ordered sets.**

  .. code-block:: cpp
    :caption: Example: Expression types

    // Polynomial expressions with real coefficients. Monomials in 'expr1'
    // and 'expr2' are products of operators with one integer index.
    libcommute::expression<double, int> expr1;
    libcommute::static_indices::expr_real<int> expr2;

    // Polynomial expressions with complex coefficients. Monomials in 'expr3'
    // and 'expr4' are products of operators with one integer
    // and one string index.
    libcommute::expression<std::complex<double>, int, std::string> expr3;
    libcommute::static_indices::expr_complex<int, std::string> expr4;

  .. rubric:: Member type aliases

  .. type:: scalar_type = ScalarType
  .. type:: index_types = std::tuple<IndexTypes...>
  .. type:: monomial_t = monomial<IndexTypes...>
  .. type:: monomials_map_t = std::map<monomial_t, ScalarType>

  .. rubric:: Constructors

  .. function:: expression() = default

   Construct a trivial expression, i.e. an expression containing no monomials.

  .. function:: template<typename S> \
                expression(expression<S, IndexTypes...> const& x)

  Construct a copy of an expression :expr:`x` with a different
  scalar type :expr:`S` by converting its coefficients to :type:`ScalarType`.

  .. function:: template<typename S> explicit expression(S const& x)

  Construct a constant expression :math:`E = x M^{(0)}`.

  .. function:: template<typename S> \
                explicit expression(S const& x, monomial_t const& monomial)

  Construct an expression made of exactly one monomial,
  :math:`E = x M^{(n)}_{i_1 i_2 \ldots i_n}`.

  The constructors listed here have limited functionality. One is supposed to
  use the :ref:`factory functions <factories>` to build expressions most
  of the time.

  .. rubric:: Copy/move-constructors and assignments

  .. function:: expression(expression const&) = default
  .. function:: expression(expression&&) noexcept = default
  .. function:: expression& operator=(expression const&) = default
  .. function:: expression& operator=(expression&&) noexcept = default
  .. function:: template<typename S> \
                expression& operator=(expression<S, IndexTypes...> const& x)

  .. rubric:: Arithmetic operations

  Two expression objects can be added, subtracted and multiplied as long as
  their :type:`IndexTypes` agree and the corresponding arithmetic operation
  is defined for their scalar types. For example, it is possible to mix real
  and complex expressions as operands of ``+``, ``-`` and ``*``.
  The unary minus is defined for an expression type iff it is defined for
  expression's scalar type. Another possibility is to mix expression objects
  with objects of other arbitrary types in binary operations. Object of the
  non-expression types are treated as constants, i.e. contributions to
  :math:`C M^{(0)}`.

  The following table shows how result types of arithmetic operations
  are calculated.

  .. list-table::
    :header-rows: 1

    * - Arithmetic operation
      - Result type
    * - :expr:`expression<S1, IT...>{} + expression<S2, IT...>{}`
      - :expr:`expression<decltype(S1{} + S2{}), IT...>`
    * - :expr:`expression<S1, IT...>{} - expression<S2, IT...>{}`
      - :expr:`expression<decltype(S1{} - S2{}), IT...>`
    * - :expr:`expression<S1, IT...>{} * expression<S2, IT...>{}`
      - :expr:`expression<decltype(S1{} * S2{}), IT...>`
    * - :expr:`-expression<S, IT...>{}`
      - :expr:`expression<decltype(-S{}), IT...>`
    * - :expr:`expression<S1, IT...>{} + S2{}`
      - :expr:`expression<decltype(S1{} + S2{}), IT...>`
    * - :expr:`expression<S1, IT...>{} - S2{}`
      - :expr:`expression<decltype(S1{} - S2{}), IT...>`
    * - :expr:`expression<S1, IT...>{} * S2{}`
      - :expr:`expression<decltype(S1{} * S2{}), IT...>`
    * - :expr:`S1{} + expression<S2, IT...>{}`
      - :expr:`expression<decltype(S1{} + S2{}), IT...>`
    * - :expr:`S1{} - expression<S2, IT...>{}`
      - :expr:`expression<decltype(S1{} - S2{}), IT...>`
    * - :expr:`S1{} * expression<S2, IT...>{}`
      - :expr:`expression<decltype(S1{} * S2{}), IT...>`

  Compound assignments ``+=``, ``-=``, ``*=`` are available under the same
  scalar type compatibility conditions between LHS and RHS. If the RHS is
  of a non-expression type :expr:`S`, *libcommute* will attempt to select
  the optimized compound operator :expr:`ScalarType::operator+=(S const& x)`
  first (similarly for ``-=``, ``*=``). If it fails,
  :expr:`ScalarType::operator+(S const& x)` and the regular assignment will
  be called instead.


  .. rubric:: :ref:`Iteration interface and transformations <expr_iteration>`

  .. class:: const_iterator

    Constant bidirectional iterator over monomial-coefficient pairs
    :math:`(M,C)` in a polynomial expression. Given an iterator :expr:`it`,
    :expr:`it->monomial` returns a constant reference to the :type:`monomial`
    object :math:`M`, and :expr:`it->coeff` is a constant reference to
    the respective coefficient :math:`C`.

  .. function:: const_iterator begin() const noexcept
                const_iterator cbegin() const noexcept

    Constant iterator to the first monomial-coefficient pair.

  .. function:: const_iterator end() const noexcept
                const_iterator cend() const noexcept

    Constant past-the-end iterator.

  .. function:: template<typename F, typename NewScalarType>  \
                friend expression<NewScalarType, IndexTypes...> \
                transform(expression const& expr, F&& f)

    Apply functor ``f`` to each monomial-coefficient pair in ``expr``.
    Return a new expression obtained by replacing coefficients in ``expr`` with
    respective values returned by ``f``. The expected signature of ``f`` is

    .. code-block:: cpp

      NewScalarType f(monomial_t const&, ScalarType const&)

    The transformed scalar type ``NewScalarType`` is automatically deduced
    from ``f``'s return type. When ``f`` returns a zero for a certain monomial,
    that monomial is excluded from the resulting expression.

  .. rubric:: Other methods and friend functions

  .. function:: std::map<monomial_t, ScalarType> const& get_monomials() const

    Direct read-only access to the list of monomials.

  .. function:: size_t size() const

    Number of monomials in this polynomial expression.

  .. function:: void clear()

    Set expression to zero by removing all monomials.

  .. function:: friend bool \
                operator==(expression const& e1, expression const& e2)

    Check if ``e1`` and ``e2`` contain identical lists of monomials.

  .. function:: friend bool \
                operator!=(expression const& e1, expression const& e2)

    Check if ``e1`` and ``e2`` contain different lists of monomials.

  .. function:: friend expression conj(expression const& expr)

    Return Hermitian conjugate of ``expr``.

  .. function:: friend std::ostream& operator<< \
                (std::ostream& os, expression const& expr)

    Output stream insertion operator.

.. type:: template<typename... IndexTypes> \
          static_indices::expr_real = \
          expression<double, IndexTypes...>;

  *Declared in <libcommute/expression/expression.hpp>*

  Shorthand type for expressions with real coefficients and statically typed
  indices.

.. type:: template<typename... IndexTypes> \
          static_indices::expr_complex = \
          expression<std::complex<double>, IndexTypes...>;

  *Declared in <libcommute/expression/expression.hpp>*

  Shorthand type for expressions with complex coefficients and statically typed
  indices.

.. _custom_scalar_type:

Custom scalar types
-------------------

Choosing :expr:`double` or :expr:`std::complex<double>` as the scalar type of
expressions covers the vast majority of practically important cases.
Nonetheless, sometimes it may be desirable to go beyond and pass a
user-defined type as :type:`ScalarType <libcommute::expression::ScalarType>`.
For instance, one may want to use types from an arbitrary-precision arithmetic
library to represent expansion coefficients of a quantum-mechanical operator.
Another potential use - making coefficients depend on a parameter, such as time
or an external field. The ``polynomial`` class from the Boost Math Toolkit or
a similar type would allow to represent the functional dependence on
the parameter while also implementing addition, subtraction and multiplication
operations.

Mathematically speaking, instances of
:type:`ScalarType <libcommute::expression::ScalarType>` are assumed to form
a ring without a multiplicative identity
(`a.k.a. rng <https://en.wikipedia.org/wiki/Ring_(mathematics)#Rng>`_). More
specifically, the set of :type:`ScalarType <libcommute::expression::ScalarType>`
instances must be

- An abelian group under addition (with binary ``operator+`` and ``operator-``).
  In particular, there must be a well-defined zero element.
- A semigroup under multiplication (``operator*``).
- Multiplication must be distributive with respect to addition.

Let us say we have a type :expr:`S` with the required algebraic properties.
Before using it as a scalar type, we must define a specialization of structure
:expr:`scalar_traits` in the namespace :expr:`libcommute` to teach *libcommute*
how to deal with the new type.

.. code-block:: cpp

  namespace libcommute {

  template<> struct scalar_traits<S> {

    // Test whether 's' is the zero element
    static bool is_zero(S const& s) { ... }

    // Make a constant of type 'S' from a double value 'x'
    static S make_const(double x) { ... }

    // OPTIONAL: Complex conjugate of 's'
    static S conj(S const& s) { ... }
  };

  }

The static member :expr:`scalar_traits<S>::conj()` is optional and will only be
called by the Hermitian conjugation function
:func:`conj() <libcommute::expression::conj()>`.

.. note::

  For the built-in floating-point types, the zero-value test method is
  implemented as

  .. code-block:: cpp

    static bool is_zero(S const& s) {
      return std::abs(s) < 100 * std::numeric_limits<S>::epsilon();
    }

  One can adjust the test and change the constant 100 to something else by
  defining a special macro :expr:`LIBCOMMUTE_FLOATING_POINT_TOL_EPS`.

.. _dyn_indices:

[C++17] Dynamically typed index sequences
-----------------------------------------

On the most basic level, index sequences of all generators found in a single
expression must agree in types with the :type:`IndexTypes
<libcommute::expression::IndexTypes>` template parameter pack. In many
situations, however, it is more natural to have generators with different
numbers/types of indices mixed in one expression. This is where the dynamically
typed index sequences step in. They are instantiations of the
:expr:`dyn_indices_generic` class template defined in a special nested namespace
:expr:`libcommute::dynamic_indices`.

.. class:: template<typename... IndexTypes> \
           dynamic_indices::dyn_indices_generic

  A wrapper around :expr:`std::vector<std::variant<IndexTypes...>>`.

  .. type:: indices_t = std::vector<std::variant<IndexTypes...>>

    Underlying dynamically typed sequence of indices.

  .. function:: dyn_indices_generic() = default

    Construct an empty index sequence.

  .. function:: dyn_indices_generic(indices_t indices)

    Construct an index sequence from a vector of indices.

  .. function:: dyn_indices_generic(dyn_indices_generic const&) = default
                dyn_indices_generic(dyn_indices_generic&&) noexcept = default
                dyn_indices_generic& \
                operator=(dyn_indices_generic const&) = default
                dyn_indices_generic& \
                operator=(dyn_indices_generic&&) noexcept = default

    Copy/move-constructors and assignments

  .. function:: size_t size() const

    Number of indices in the sequence.

  .. function:: explicit operator indices_t const& () const

    Explicit cast to the underlying index sequence type.

  .. function:: friend bool operator==(dyn_indices_generic const& ind1, \
                                       dyn_indices_generic const& ind2)
                friend bool operator!=(dyn_indices_generic const& ind1, \
                                       dyn_indices_generic const& ind2)
                friend bool operator<(dyn_indices_generic const& ind1, \
                                      dyn_indices_generic const& ind2)
                friend bool operator>(dyn_indices_generic const& ind1, \
                                      dyn_indices_generic const& ind2)

    Compare two dynamic index sequences :expr:`ind1` and :expr:`ind2`.
    These operators compare sequences' lengths first, and in the case of equal
    lengths call the `corresponding methods of std::vector
    <https://en.cppreference.com/w/cpp/container/vector/operator_cmp>`_.

  .. function:: friend std::ostream & operator<< \
                (std::ostream & os, dyn_indices_generic const& ind)

    Output stream insertion operator.

.. code-block:: cpp
  :caption: Dynamic index sequence example

  #include <libcommute/expression/dyn_indices.hpp>

  // A user-defined index type
  enum spin {up, down};

  // Our own dynamically typed index sequence (mind the namespace!)
  using my_dyn_indices =
    libcommute::dynamic_indices::dyn_indices_generic<int, std::string, spin>;

  // An expression with dynamically typed indices
  libcommute::expression<double, my_dyn_indices> dyn_expr;

:expr:`dyn_expr` is an expression with dynamically typed indices. It can
contain generators that effectively carry a variable number of indices,
each of type :expr:`int`, :expr:`std::string` or of the user-defined enumeration
type :expr:`spin`.

There is also a special set of :ref:`factory functions <factories_dyn>` defined
in namespaces nested under :expr:`libcommute::dynamic_indices`. Those return
commonly used QM operators with the dynamically typed indices.
Some related type aliases are declared in the same namespace for the sake of
convenience.

.. type:: dynamic_indices::dyn_indices = \
          dyn_indices_generic<int, std::string>

  *Declared in <libcommute/expression/dyn_indices.hpp>*

  Dynamic mixture of integer and string indices.

.. type:: dynamic_indices::expr_real \
          = expression<double, dyn_indices>
          dynamic_indices::expr_complex \
          = expression<std::complex<double>, dyn_indices>

  *Declared in <libcommute/expression/expression.hpp>*

  Real/complex expression shorthand types.

.. _expr_iteration:

Iteration interface and transformations
---------------------------------------

Both expressions and their constituent monomials can be easily iterated over
using STL-compatible :class:`iterators <expression::const_iterator>` or
range-based ``for``-loops.
This allows for writing complex expression analysis algorithms.

.. code-block:: cpp
  :caption: Expression iteration example

    using namespace libcommute;

    // We are going to analyse the structure of this expression
    expression<double, int> E;

    //
    // Fill expression 'E' ...
    //

    // Iterate over all monomial-coefficient pairs in 'E'
    for(auto const& mc : E) {
      std::cout << "Coefficient: " << mc.coeff << "\n";
      std::cout << "Monomial: ";

      // Iterate over algebra generators in current monomial
      for(auto const& g : mc.monomial) {

        if(is_fermion(g)) { // Print information about fermionic operators
          auto const& f = dynamic_cast<generator_fermion<int> const&>(g);
          std::cout << (f.dagger() ? "  c^+(" : "  c(");
          std::cout << std::get<0>(g.indices());
          std::cout << ")";
        }

        if(is_fermion(g)) { // Print information about bosonic operators
          auto const& a = dynamic_cast<generator_boson<int> const&>(g);
          std::cout << (a.dagger() ? "  a^+(" : "  a(");
          std::cout << std::get<0>(g.indices());
          std::cout << ")";
        }

        if(is_spin(g)) { // Print information about spin operators
          auto const& s = dynamic_cast<generator_spin<int> const&>(g);
          switch(s.component()) {
            case plus:
              std::cout << "S_+("; break;
            case minus:
              std::cout << "S_-("; break;
            case z:
              std::cout << "S_z("; break;
          }
          std::cout << std::get<0>(g.indices());
          std::cout << ")";
        }

      }
      std::cout << std::endl;
    }

.. note:: Only the constant iterators are implemented by :type:`expression` and
          :type:`monomial`.

Another common task is building a new expression out of an existing one by
changing some of monomial coefficients. Removing some monomials (for instance,
those corresponding to particle interaction terms) is a special case that
amounts to setting coefficients of the unwanted monomials to zero. In the
following example we show how to use function
:func:`libcommute::transform() <libcommute::expression::transform>` to change
Hamiltonian of a finite atomic chain

.. math::

  \hat H = v \sum_{a=1}^{N-1} (c^\dagger_a c_{a+1} + c^\dagger_{a+1} c_a)

into the Su-Schrieffer-Heeger (SSH) model, where the hopping constant is
taken to be different on odd and even chain links.

.. code-block:: cpp
  :caption: :func:`transform() <libcommute::expression::transform>` example

    using namespace libcommute;

    // Hamiltonian of the atomic chain
    expression<double, int> H;

    using static_indices::c_dag;
    using static_indices::c;

    // Length of the chain
    const int N = 10;
    // Hopping matrix element
    const double v = 1.0;

    // Add hopping terms
    for(int a = 0; a < N - 1; ++a) {
      H += v * (c_dag(a) * c(a + 1) + c_dag(a + 1) * c(a));
    }

    std::cout << "H = " << H << std::endl;

    // Construct the Su-Schrieffer-Heeger (SSH) model by changing hopping
    // constants on all even links of the chain.

    // Hopping matrix element on the even links
    const double w = 0.5;

    // Transformation operation
    auto update_hopping_element = [v, w](decltype(H)::monomial_t const& m,
                                  double coeff) -> double {
      // The monomial 'm' we are expecting here is a product c^\dagger(a1) c(a2)

      // Site index of c^\dagger
      int a1 = std::get<0>(m[0].indices());
      // Site index of c
      int a2 = std::get<0>(m[1].indices());
      // Index of link connecting a1 and a2
      int link = std::min(a1, a2);

      // Return the updated matrix element
      return (link % 2 == 1 ? w : v);
    };

    auto H_SSH = transform(H, update_hopping_element);

    std::cout << "H_SSH = " << H_SSH << std::endl;

.. _hc:

:math:`\pm H.c.` notation
-------------------------

It is common to use the :math:`\pm H.c.` ("plus/minus Hermitian conjugate")
notation to abbreviate analytical expressions when writing down Hamiltonians
and other Hermitian operators. The same trick works with *libcommute*'s
expressions, thanks to the ``hc`` singleton object.

Here is a simple usage example that constructs two operators

.. math::

  \hat H_A = c^\dagger_1 c_2 + c^\dagger_2 c_1 = c^\dagger_1 c_2 + H.c.,\\
  \hat H_B = i(c^\dagger_1 c_2 - c^\dagger_2 c_1) = ic^\dagger_1 c_2 - H.c.

.. code-block:: cpp

  #include <libcommute/expression/hc.hpp>

  auto H_A = c_dag(1) * c(2) + hc;

  std::complex<double> I(0,1);
  auto H_B = I * c_dag(1) * c(2) - hc;
