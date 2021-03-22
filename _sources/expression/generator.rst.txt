.. _generator:

Algebra generators
==================

.. default-domain:: cpp

.. namespace:: libcommute

Algebra generators, such as creation/annihilation operators :math:`c^\dagger`/
:math:`c` in `fermionic and bosonic algebras`__, are atomic structural units of
any expression. Within *libcommute*'s framework, algebra generators are classes
derived from the abstract base :type:`libcommute::generator`. Every generator
carries an index sequence (possibly of zero length), and the respective
index types must be passed as template parameters to
:type:`libcommute::generator`. The index types must be less-comparable and
form strictly ordered sets so that tuples of the indices (index sequences) are
also less-comparable and form a strictly ordered set.

The few methods that derived classes have to override serve a multitude
of purposes.

.. _CCR_and_CAR: https://en.wikipedia.org/wiki/CCR_and_CAR_algebras
__ CCR_and_CAR_

- Assign a numerical ID shared by all generators of a particular algebra
  (fermions, bosons, spin operators or a user-defined algebra) and unique
  to that algebra.

- Establish a canonical ordering of generators belonging to the same algebra.

- Describe commutation/anti-commutation rules that can be used to simplify a
  product of two generators and to put them into the canonical order.

- Describe further simplification rules applicable to products of 3 and more
  generators (currently, only one simplification of this kind is implemented).

Generators of different algebras always commute. In a canonically ordered
product, generators are placed in the non-decreasing order of their algebra IDs.

The following table summarizes information about predefined generators.

.. list-table::
  :header-rows: 1

  * - Algebra
    - Generator type
    - Algebra ID
  * - Fermions :math:`c^\dagger_i`/:math:`c_i`
    - :expr:`libcommute::generator_fermion`
    - :expr:`libcommute::fermion`
  * - Bosons :math:`a^\dagger_i`/:math:`a_i`
    - :expr:`libcommute::generator_boson`
    - :expr:`libcommute::boson`
  * - Spins :math:`S_\pm`/:math:`S_z`
    - :expr:`libcommute::generator_spin`
    - :expr:`libcommute::spin`
  * - User-defined algebra
    - A class derived from :type:`libcommute::generator`
    - >= :expr:`libcommute::LIBCOMMUTE_MIN_USER_DEFINED_ALGEBRA_ID`

Integer constants :expr:`fermion`, :expr:`boson`, :expr:`spin` and
:expr:`LIBCOMMUTE_MIN_USER_DEFINED_ALGEBRA_ID`
mentioned in the 3rd column are defined in *<libcommute/algebra_ids.hpp>*.

.. _algebra_ids:

.. var:: static constexpr int fermion = -3
.. var:: static constexpr int boson = -2
.. var:: static constexpr int spin = -1
.. var:: static constexpr int LIBCOMMUTE_MIN_USER_DEFINED_ALGEBRA_ID = 0

.. _gen_base:

``generator``: abstract base class for algebra generators
---------------------------------------------------------

.. class:: template<typename... IndexTypes> generator

  *Defined in <libcommute/expression/generator.hpp>*

  The abstract base class for algebra generator types.

  :type:`IndexTypes` - types of indices carried by this generator.

  .. rubric:: Member type aliases

  .. type:: index_types = std::tuple<IndexTypes...>

    Index tuple type.

  .. type:: linear_function_t = linear_function<std::unique_ptr<generator>>

    Linear combination of generators. This type is used by various methods
    dealing with transformations of generator products.

  .. rubric:: Constructor

  .. function:: template<typename... Args> generator(Args&&... indices)

    Construct generator with given indices.

  .. rubric:: Copy/move-constructors, assignments and destructor

  .. function:: generator(generator const&) = default
  .. function:: generator(generator&&) noexcept = default
  .. function:: generator& operator=(generator const&) = default
  .. function:: generator& operator=(generator&&) noexcept = default
  .. function:: virtual ~generator()
  .. function:: virtual std::unique_ptr<generator> clone() const = 0

    Virtual copy-constructor. Makes a copy of this generator managed by a
    unique pointer.

  .. rubric:: Algebra ID

  .. function:: virtual int algebra_id() const = 0

    Get the ID of the algebra this generator belongs to.

  .. rubric:: Index sequence

  .. function:: index_types const& indices() const

    Read-only access to the index tuple carried by the generator.

  .. rubric:: Canonical ordering

  .. function:: protected virtual bool equal(generator const& g) const
                protected virtual bool less(generator const& g) const
                protected virtual bool greater(generator const& g) const

    These methods can be overridden by the derived classes to establish
    the canonical order of :expr:`g` w.r.t. :expr:`*this` assuming both
    generators belong to the same algebra. The default implementation compares
    index tuples of :expr:`*this` and :expr:`g`.

  .. function:: friend bool operator==(generator const& g1, generator const& g2)
                friend bool operator!=(generator const& g1, generator const& g2)
                friend bool operator<(generator const& g1, generator const& g2)
                friend bool operator>(generator const& g1, generator const& g2)

    Comparison operators for a pair of generators. First, they compare algebra
    IDs of :expr:`g1` and :expr:`g2`. If those are equal, :expr:`g1.equal(g2)`,
    :expr:`g1.less(g2)` or :expr:`g1.greater(g2)` is called.

  .. rubric:: Product simplification/transformation

  .. function:: virtual double swap_with\
                (generator const& g2, linear_function_t & f) const = 0

    Given a pair of generators :math:`g_1` (:expr:`*this`) and :math:`g_2`
    such that :math:`g_1 > g_2`, :expr:`swap_with()` must signal what
    transformation :math:`g_1 g_2 \mapsto c g_2  g_1 + f(g)` should
    be applied to the product :math:`g_1 g_2` in order to put it into the
    canonical order. :expr:`swap_with()` returns the constant :math:`c` and
    writes the linear function of generators :math:`f(g)` into its second
    argument. :math:`c` is allowed to be zero.

  .. function:: virtual bool \
                simplify_prod(generator const& g2, linear_function_t & f) const

    Given a pair of generators :math:`g_1` (:expr:`*this`) and :math:`g_2` such
    that :math:`g_1 g_2` is in the canonical order (:math:`g_1 \leq g_2`),
    optionally apply a simplifying transformation :math:`g_1 g_2 \mapsto f(g)`.
    If a simplification is actually possible, :expr:`simplified_prod()` must
    return :expr:`true` and write the linear function :math:`f(g)` into its
    second argument. Otherwise return :expr:`false`.

    The default implementation always returns :expr:`false`.

  .. function:: virtual bool \
                reduce_power(int power, linear_function_t & f) const

    Given a generator :math:`g_1` (:expr:`*this`) and a power :math:`p > 2`
    (:expr:`power`), optionally apply a simplifying transformation
    :math:`g_1^p \mapsto f(g)`. If a simplification is actually possible,
    :expr:`reduce_power()` must return :expr:`true` and write the linear
    function :math:`f(g)` into its second argument.
    Otherwise return :expr:`false`.

    The default implementation always returns :expr:`false`.

    .. note:: Simplifications for :math:`p = 2` must be carried out by
              :expr:`simplify_prod()`.

  .. rubric:: Other methods

  .. function:: virtual void conj(linear_function_t & f)

    Return the Hermitian conjugate of generator as a linear function of other
    generators (write the result into :expr:`f`). The default implementation
    returns the generator itself.

  .. function:: friend std::ostream & operator<<\
                (std::ostream & os, generator const& g)

     Output stream insertion operator. Calls :expr:`g.print(os)`.

  .. function:: protected virtual std::ostream & print(std::ostream & os) const

    Virtual stream output function to be overridden by the derived classes.


.. struct:: template<typename T> linear_function

  *Defined in <libcommute/utility.hpp>*

  A linear function of objects of type :expr:`T` with :expr:`double`
  coefficients,

  .. math::

    f(x_1, \ldots, x_n) = c + c_1 x_1 + \ldots + c_n x_n.

  .. type:: basis_type = T

  .. member:: double const_term = 0;

    Constant term :math:`c`.

  .. member:: std::vector<std::pair<T, double>> terms

    List of pairs :math:`(x_1, c_1), \ldots, (x_n, c_n)`.

  .. function:: linear_function() = default

    Construct an identically vanishing function :math:`f(x_1, \ldots, x_n) = 0`.

  .. function:: linear_function(double const_term)

    Construct a constant function :math:`f(x_1, \ldots, x_n) = c`.

  .. function:: linear_function(double const_term, Args&&... args)

    Construct a linear function from a sequence of arguments
    :math:`c, x_1, c_1, x_2, c_2, \ldots, x_n, c_n`.

  .. function:: linear_function(double const_term, \
                                std::vector<std::pair<T, double>> terms)

    Construct a linear function from a constant term and a list of pairs
    :math:`(x_1, c_1), \ldots, (x_n, c_n)`.

  .. function:: void set(double const_term, Args&&... args)

    Clear all terms and replace them with a sequence of arguments
    :math:`c, x_1, c_1, x_2, c_2, \ldots, x_n, c_n`.

  .. function:: bool vanishing() const

    Is this linear function identically zero?

.. _generator_fermion:

Fermions
--------

Fermionic algebra is generated by creation and annihilation operators
:math:`c_i^\dagger`/:math:`c_i` with canonical anti-commutation relations

.. math::

  \{c_i, c^\dagger_j\} &= \delta_{ij}, \\
  \{c_i, c_j\} &= 0, \\
  \{c^\dagger_i, c^\dagger_j\} &= 0.

The canonical order is defined according to

.. math::

  c^\dagger_{i_1} < c^\dagger_{i_2} < c^\dagger_{i_3}
  < c_{i_3} < c_{i_2} < c_{i_1},

where index sequences :math:`i_k` satisfy :math:`i_1 < i_2 < i_3`.
In other words,

- Creation operators precede annihilation operators;
- Creation operator with the smallest index sequence comes first;
- Annihilation operator with the smallest index sequence comes last.

.. class:: template<typename... IndexTypes> \
           generator_fermion : public generator<IndexTypes...>

  *Defined in <libcommute/expression/generator_fermion.hpp>*

  .. rubric:: Part of the interface not inherited from / identical to
              :type:`libcommute::generator`.

  .. function:: bool dagger() const

    Returns :expr:`true` for :math:`c^\dagger` and :expr:`false` for :math:`c`.

.. function:: template<typename... IndexTypes> \
              generator_fermion<IndexTypes...> \
              static_indices::make_fermion(bool dagger, \
              IndexTypes&&... indices)

  *Defined in <libcommute/expression/generator_fermion.hpp>*

  Make a fermionic creation (:expr:`dagger = true`) or annihilation
  (:expr:`dagger = false`) operator with given indices.

.. function:: template<typename... IndexTypes> \
              generator_fermion<IndexTypes...> \
              dynamic_indices::make_fermion(bool dagger, \
              IndexTypes&&... indices)

  *Defined in <libcommute/expression/generator_fermion.hpp>*

  Make a fermionic creation (:expr:`dagger = true`) or annihilation
  (:expr:`dagger = false`) operator with a given
  :ref:`dynamic index sequence <dyn_indices>`.

.. function:: template<typename... IndexTypes> \
              bool is_fermion(generator<IndexTypes...> const& gen)

  *Defined in <libcommute/expression/generator_fermion.hpp>*

  Detect if :expr:`gen` points to a generator of the fermionic algebra.

.. code-block:: cpp

  using namespace libcommute::static_indices;

  // Make c^\dagger_{1,up}
  auto g = make_fermion(true, 1, "up");

  // ...

  // If 'g' is a fermionic generator, print whether it is a creation
  // or annihilation operator.
  if(is_fermion(g)) {
    auto const& f = dynamic_cast<generator_fermion<int, std::string> const&>(g);
    std::cout << (f.dagger() ? "creation" : "annihilation") << std::endl;
  }

.. _generator_boson:

Bosons
------

Bosonic algebra is generated by creation and annihilation operators
:math:`a_i^\dagger`/:math:`a_i` with canonical commutation relations

.. math::

  [a_i, a^\dagger_j] &= \delta_{ij}, \\
  [a_i, a_j] &= 0, \\
  [a^\dagger_i, a^\dagger_j] &= 0.

The canonical order is defined according to

.. math::

  a^\dagger_{i_1} < a^\dagger_{i_2} < a^\dagger_{i_3}
  < a_{i_3} < a_{i_2} < a_{i_1},

where index sequences :math:`i_k` satisfy :math:`i_1 < i_2 < i_3`.
In other words,

- Creation operators precede annihilation operators;
- Creation operator with the smallest index sequence comes first;
- Annihilation operator with the smallest index sequence comes last.

.. class:: template<typename... IndexTypes> \
           generator_boson : public generator<IndexTypes...>

  *Defined in <libcommute/expression/generator_boson.hpp>*

  .. rubric:: Part of the interface not inherited from / identical to
              :type:`libcommute::generator`.

  .. function:: bool dagger() const

    Returns :expr:`true` for :math:`a^\dagger` and :expr:`false` for :math:`a`.

.. function:: template<typename... IndexTypes> \
              generator_boson<IndexTypes...> \
              static_indices::make_boson(bool dagger, \
              IndexTypes&&... indices)

  *Defined in <libcommute/expression/generator_boson.hpp>*

  Make a bosonic creation (:expr:`dagger = true`) or annihilation
  (:expr:`dagger = false`) operator with given indices.

.. function:: template<typename... IndexTypes> \
              generator_fermion<IndexTypes...> \
              dynamic_indices::make_boson(bool dagger, \
              IndexTypes&&... indices)

  *Defined in <libcommute/expression/generator_boson.hpp>*

  Make a bosonic creation (:expr:`dagger = true`) or annihilation
  (:expr:`dagger = false`) operator with a given
  :ref:`dynamic index sequence <dyn_indices>`.

.. function:: template<typename... IndexTypes> \
              bool is_boson(generator<IndexTypes...> const& gen)

  *Defined in <libcommute/expression/generator_boson.hpp>*

  Detect if :expr:`gen` points to a generator of the bosonic algebra.

.. code-block:: cpp

  using namespace libcommute::static_indices;

  // Make a^\dagger_1
  auto g = make_boson(true, 1);

  // ...

  // If 'g' is a bosonic generator, print whether it is a creation or
  // annihilation operator.
  if(is_boson(g)) {
    auto const& b = dynamic_cast<generator_boson<int> const&>(g);
    std::cout << (b.dagger() ? "creation" : "annihilation") << std::endl;
  }

.. _generator_spin:

Spins
-----

*libcommute* supports algebra of spin operators for :math:`S = 1/2` as well as
for higher integer and half-integer spins. Generators of spin algebras with
different :math:`S` share the same algebra ID and are distinguished by an extra
integer data member :expr:`multiplicity` equal to :math:`2S+1`.
For a fixed :math:`S` and a set of indices, the spin algebra is generated
by the triplet of operators :math:`S_+`, :math:`S_-` and :math:`S_z` subject to
the following commutation relations.

.. math::

  [S_+, S_-] &= 2 S_z, \\
  [S_z, S_+] &= S_+, \\
  [S_z, S_-] &= -S_-.

.. note::

  Using :math:`S_\pm` instead of :math:`S_x`, :math:`S_y` as algebra generators
  is beneficial because all coefficients in the commutation relations above are
  real. :math:`S_x`/:math:`S_y` would necessitate the complex scalar types in
  all :class:`libcommute::expression` objects.

The canonical order is defined according to

.. math::

  S_{1,+}^{S=1/2} < S_{1,-}^{S=1/2} < S_{1,z}^{S=1/2} <
  S_{2,+}^{S=1/2} < S_{2,-}^{S=1/2} < S_{2,z}^{S=1/2} < \\ <
  S_{2,+}^{S=3/2} < S_{2,-}^{S=3/2} < S_{2,z}^{S=3/2} <
  S_{2,+}^{S=3/2} < S_{2,-}^{S=3/2} < S_{2,z}^{S=3/2}.

In other words,

- Operators with lower :math:`S` precede operators with higher :math:`S`.
- Among operators with the same :math:`S`, the operator with the smallest
  index sequence comes first.
- Among operators with the same :math:`S` and index sequence, :math:`S_+` comes
  first followed by :math:`S_-` and eventually by :math:`S_z`.

.. enum:: spin_component : int

  Component of spin operator.

  .. enumerator:: plus = 0

    :math:`S_+`.

  .. enumerator:: minus = 1

    :math:`S_-`.

  .. enumerator:: z = 2

    :math:`S_z`.

.. class:: template<typename... IndexTypes> \
           generator_spin : public generator<IndexTypes...>

  *Defined in <libcommute/expression/generator_spin.hpp>*

  .. rubric:: Part of the interface not inherited from / identical to
              :type:`libcommute::generator`.

  .. function:: template<typename... Args> \
                generator_spin(spin_component c, Args&&... indices)

    Construct generator :math:`S_+`, :math:`S_-` or :math:`S_z` for spin
    :math:`S=1/2` with given indices.

  .. function:: template<typename... Args> \
                generator_spin(double spin, spin_component c, \
                Args&&... indices)

    Construct generator :math:`S_+`, :math:`S_-` or :math:`S_z` for spin
    :expr:`spin` with given indices.

  .. function:: double spin() const

    Read-only access to generator's spin :math:`S`.

  .. function:: int multiplicity() const

    Read-only access to generator's multiplicity :math:`2S+1`.

  .. function:: libcommute::spin_component component() const

    Is this generator :math:`S_+`, :math:`S_-` or :math:`S_z`?

.. function:: template<typename... IndexTypes> \
              generator_spin<IndexTypes...> \
              static_indices::make_spin( \
              spin_component c, IndexTypes&&... indices)

  *Defined in <libcommute/expression/generator_spin.hpp>*

  Make generator :math:`S_+`, :math:`S_-` or :math:`S_z` for spin
  :math:`S=1/2` with given indices.

.. function:: template<typename... IndexTypes> \
              generator_spin<IndexTypes...> \
              static_indices::make_spin(double spin, \
              spin_component c, IndexTypes&&... indices)

  *Defined in <libcommute/expression/generator_spin.hpp>*

  Make generator :math:`S_+`, :math:`S_-` or :math:`S_z` for spin
  :expr:`spin` with given indices.

.. function:: template<typename... IndexTypes> \
              generator_spin<dyn_indices> \
              dynamic_indices::make_spin( \
              spin_component c, IndexTypes&&... indices)

  *Defined in <libcommute/expression/generator_spin.hpp>*

  Make generator :math:`S_+`, :math:`S_-` or :math:`S_z` for spin
  :math:`S=1/2` with a given :ref:`dynamic index sequence <dyn_indices>`.

.. function:: template<typename... IndexTypes> \
              generator_spin<dyn_indices> \
              libcommute::dynamic_indices::make_spin( \
              double spin, libcommute::spin_component c, \
              IndexTypes&&... indices)

  *Defined in <libcommute/expression/generator_spin.hpp>*

  Make generator :math:`S_+`, :math:`S_-` or :math:`S_z` for spin
  :expr:`spin` with a given :ref:`dynamic index sequence <dyn_indices>`.

.. function:: template<typename... IndexTypes> \
              bool libcommute::is_spin( \
              libcommute::generator<IndexTypes...> const& gen)

  *Defined in <libcommute/expression/generator_spin.hpp>*

  Detect if :expr:`gen` points to a generator of the spin algebra.

.. code-block:: cpp

  using namespace libcommute::static_indices;

  // Make S^{J=1}_{1,+}
  auto g = make_spin(1.0, libcommute::plus, 1);

  // ...

  // If 'g' is a spin algebra generator, print its properties.
  if(is_spin(g)) {
    auto const& s = dynamic_cast<generator_spin<int> const&>(g);

    std::cout << "J = " << s.spin() << std::endl;
    std::cout << "2J+1 = " << s.multiplicity() << std::endl;
    switch(s.component()) {
      case libcommute::plus:
        std::cout << "+" << std::endl;
        break;
      case libcommute::minus:
        std::cout << "-" << std::endl;
        break;
      case libcommute::z:
        std::cout << "z" << std::endl;
        break;
    }
  }
