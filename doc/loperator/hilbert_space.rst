.. _hilbert_space:

Finite-dimensional Hilbert spaces
=================================

.. default-domain:: cpp
.. namespace:: libcommute

Class :class:`hilbert_space` is an abstraction of a finite-dimensional state
space of a quantum system. It contains information needed to construct a
:ref:`linear operator object <loperator>` out of a
:ref:`polynomial expression <expression>`.

A :class:`hilbert_space` is defined as an (ordered) direct product of
:ref:`elementary spaces <elementary_spaces>` :math:`\mathcal{H}_i`,

.. math::

  \mathcal{H} = \mathcal{H}_1 \otimes \mathcal{H}_2 \otimes \ldots \otimes
                \mathcal{H}_N.

An elementary space associated with an :ref:`algebra generator <generator>`
:math:`g` is a vector space of dimension :math:`2^b`, where :math:`b` is
the smallest integer sufficient to construct a matrix representation of
:math:`g`.
For instance, an elementary space associated with a fermionic
creation/annihilation operator is two-dimensional (:math:`b = 1`); It is
spanned by the occupation number states :math:`|0\rangle` and :math:`|1\rangle`.
For a spin-1 operator, we have a 4-dimensional elementary space.
There are three basis vectors :math:`|m=0\rangle`, :math:`|m=\pm 1\rangle` of
the irreducible representation, but one has to round the dimension up to the
nearest power of 2 (:math:`b = 2`). Finally, in the case of a bosonic generator
one has to truncate the infinite-dimensional state space and manually set
:math:`b\geq 1`.

The ordering of the elementary spaces in the product is established by the
algebra IDs of the generators :math:`g` these spaces are associated with
(the smaller algebra ID comes first). Two elementary spaces sharing the same ID
are ordered according to the :ref:`rules specific to the corresponding
algebra <elementary_spaces>`.

The requirement of dimensions of the elementary spaces being powers of two is
not arbitrary. It is justified by the way basis states of the full space
:math:`\mathcal{H}` are enumerated. Each basis state of :math:`\mathcal{H}` is
a direct product of basis vectors of the elementary spaces,

.. math::

  |n\rangle_\mathcal{H} = |n_1\rangle_{\mathcal{H}_1} \otimes
                          |n_2\rangle_{\mathcal{H}_2} \otimes\ldots\otimes
                          |n_N\rangle_{\mathcal{H}_N}.

In the code, the basis vectors are represented by 64-bit unsigned integers
(:type:`sv_index_type`).
The binary form of :math:`|n\rangle_\mathcal{H}` is
then a concatenation of binary forms of :math:`|n_i\rangle_{\mathcal{H}_i}`.
For example, the following picture shows memory representation of basis
state :math:`|90\rangle_\mathcal{H} = |0\rangle_{\mathcal{H}_1} \otimes
|5\rangle_{\mathcal{H}_2} \otimes |1\rangle_{\mathcal{H}_3} \otimes
|1\rangle_{\mathcal{H}_4}`. Cells (bits) of the same color belong to the same
elementary space and the higher blank cells 7-63 are unused -- set to zero.

.. image:: ../images/basis_state.svg
  :width: 800

The simplest way to construct a :class:`hilbert_space` object is by calling
:func:`make_hilbert_space()` on an expression.

.. code:: cpp

  using namespace libcommute;
  using namespace static_indices; // For n()

  auto H = 2.0 * n("up", 0) * n("dn", 0);

  // Construct a 4-dimensional Hilbert space, which is a product of two
  // fermionic elementary spaces (for {"up", 0} and {"dn", 0}).
  auto hs = make_hilbert_space(H);

:func:`make_hilbert_space()` is a convenience function that forwards its
arguments to one of :class:`hilbert_space`'s constructors. That constructor
iterates over all generators found in the expression and adds their associated
elementary spaces into the product :math:`\mathcal{H}`. Sometimes, it may
require extra bits of information to translate a generator into an elementary
space. Perhaps the most prominent example here is a bosonic elementary space,
whose truncated dimension must be set by the user. It is possible to
:ref:`customize the Hilbert space construction procedure <es_constructor>`
by passing an extra argument to :func:`make_hilbert_space()`. The following
code snippet shows how to set the truncated space dimension for all bosonic
generators at once.

.. code:: cpp

    using namespace libcommute;
    using namespace static_indices; // For a_dag() and a()

    auto H = 2.0 * (a_dag(0) * a(0) - 0.5) + 3.0 * (a_dag(1) * a(1) - 0.5);

    // hs is a direct product of two bosonic elementary spaces, each with
    // dimension 2^4 = 16.
    auto hs = make_hilbert_space(H, boson_es_constructor(4));

Other, more refined ways to create a Hilbert space are (a) to explicitly provide
a list of elementary spaces or (b) to start from an empty product and add
elementary spaces one by one. One might have to resort to these approaches when
some elementary spaces must be added into the product, but their corresponding
generators are not necessarily found in ``H``.

.. code:: cpp

    using namespace libcommute;
    using namespace static_indices; // For make_space_*()

    // A product of three elementary spaces
    // std::string and int are index types of generators
    hilbert_space<std::string, int> hs1(
      make_space_fermion("dn", 0),   // Fermion
      make_space_boson(4, "x", 0),   // Boson truncated to dim = 2^4
      make_space_spin(0.5, "i", 0)   // Spin-1/2
    );

    // Empty space, to be filled later
    hilbert_space<int> hs2;
    // Fill the space
    hs2.add(make_space_fermion(0));  // Add a fermion
    hs2.add(make_space_boson(4, 0)); // Add a boson


The order in which the elementary spaces are passed to the constructor or added
does not matter -- they will be reordered automatically.

.. class:: template<typename... IndexTypes> hilbert_space

  *Defined in <libcommute/loperator/hilbert_space.hpp>*

  State space of a quantum system as a direct product of
  :class:`elementary spaces <elementary_space>`.

  Parameter pack :type:`IndexTypes` must agree with that of
  :class:`elementary_space` and/or :class:`libcommute::expression` objects,
  which are used to construct this Hilbert space.

  .. rubric:: Constructors

  .. function:: expression() = default

    Construct an empty space.

  .. function:: template<typename... Args> \
                explicit hilbert_space(Args&&... args)

    Construct from a list of elementary spaces. The elementary spaces need not
    be given in any particular order.
    Throws :struct:`hilbert_space_too_big` if all elementary spaces together
    would overflow the 64-bit integer type of the basis state index.

  .. function:: template<typename ScalarType, \
                         typename ESConstructor = default_es_constructor> \
            hilbert_space(libcommute::expression<ScalarType, IndexTypes...> \
                            const& expr, \
                          ESConstructor&& es_constr = {})

    Inspect an expression :expr:`expr` and collect all elementary spaces
    associated with algebra generators found in :expr:`expr`.
    Construction of the elementary spaces is performed by the functor
    :expr:`es_constr`.
    Throws :struct:`hilbert_space_too_big` if all collected elementary spaces
    together would overflow the 64-bit integer type of the basis state index.


  .. rubric:: Copy/move-constructors and assignments

  .. function:: hilbert_space(hilbert_space const&)
  .. function:: hilbert_space(hilbert_space&&) noexcept = default
  .. function:: hilbert_space& operator=(hilbert_space const&)
  .. function:: hilbert_space& operator=(hilbert_space&&) noexcept = default

  .. rubric:: Other methods and friend functions

  .. function:: friend bool operator==(hilbert_space const& hs1, \
                                       hilbert_space const& hs2)
                friend bool operator!=(hilbert_space const& hs1, \
                                       hilbert_space const& hs2)

    Check that two Hilbert spaces have an identical/different structure.

  .. function:: void add(elementary_space<IndexTypes...> const& es)

    Insert a new elementary space into the product. Throws
    :struct:`elementary_space_exists` if an elementary space equivalent to
    :expr:`es` is already part of the product.
    Throws :struct:`hilbert_space_too_big` if adding :expr:`es` into the product
    would overflow the 64-bit integer type of the basis state index.

  .. function:: bool has(elementary_space<IndexTypes...> const& es) const

    Is elementary space :expr:`es` part of the product?

  .. function:: int index(elementary_space<IndexTypes...> const& es) const

    Position of a given elementary space in the product.

  .. function:: std::pair<int, int> bit_range( \
                elementary_space<IndexTypes...> const& es) const

    Returns the range of bits in the binary representation of a
    basis state index that is occupied by the elementary space :expr:`es`.
    The range is returned as a pair (first_bit, last_bit).
    Throws :struct:`elementary_space_not_found` if :expr:`es` is not part of
    the product.

  .. function:: std::pair<int, int> const& algebra_bit_range(int algebra_id) \
                const

    Return the range of bits in the binary representation of a basis state
    index that is occupied by all elementary spaces with a given
    :ref:`algebra ID <generator>`. The range is always contiguous because
    elementary spaces with the same algebra ID are grouped together in the
    product. Throws :class:`std::runtime_error` if there is no elementary spaces
    with the given ID in the product.

  .. function:: std::size_t size() const

    The number of elementary spaces in the product.

  .. function:: int total_n_bits() const

    The total number of used bits in the binary representation of a basis state
    index.

  .. function:: std::size_t dim() const
                friend std::size_t get_dim(hilbert_space const& hs)

    The dimension of this Hilbert space computed as a product of dimensions
    of the elementary spaces.

  .. function:: template<typename Functor> \
                friend void foreach(hilbert_space const& hs, Functor&& f)

    Apply functor :expr:`f` to all basis state indices in :expr:`hs`.
    :expr:`f` must accept one argument of type :type:`sv_index_type`.

  .. function:: sv_index_type \
                basis_state_index(elementary_space<IndexTypes...> const& es, \
                                  sv_index_type n)

    Given an elementary space :expr:`es` and an index :expr:`n` of a basis state
    within it, return the corresponding basis state index within the full
    Hilbert space.

  .. rubric:: Exception types

  .. struct:: elementary_space_exists : public std::runtime_error

    Thrown when one tries to add an elementary space that is already part of
    the product.

  .. struct:: elementary_space_not_found : public std::runtime_error

    Given elementary space is not part of the product.

  .. struct:: hilbert_space_too_big : public std::runtime_error

    The total basis state index size exceeds 64 bits.

.. function:: template<typename ScalarType, \
                       typename... IndexTypes, \
                       typename ESConstructor = default_es_constructor> \
              hilbert_space<IndexTypes...> \
              make_hilbert_space(\
                expression<ScalarType, IndexTypes...> const& expr, \
                ESConstructor&& es_constr = {})

  *Defined in <libcommute/loperator/hilbert_space.hpp>*

  A helper factory function that constructs an :class:`hilbert_space` instance
  from an expression :expr:`expr` using an
  :ref:`elementary space constructor <es_constructor>`. This function is a
  more convenient equivalent of one of :class:`hilbert_space`'s constructors.

.. _elementary_spaces:

Elementary spaces
-----------------

An elementary space has an :ref:`algebra ID <generator>` assigned to it and
carries a tuple of indices. Together, these two pieces of information link
the elementary space to algebra generators acting in it.

.. class:: template<typename... IndexTypes> elementary_space

  *Defined in <libcommute/loperator/elementary_space.hpp>*

  Abstract base class for elementary spaces. :type:`IndexTypes` are index types
  of the associated algebra generators.

  .. type:: index_types = std::tuple<IndexTypes...>

    Index tuple type.

  .. rubric:: Constructors

  .. function:: template<typename... Args> elementary_space(Args&&... indices)
                elementary_space(index_types const& indices)
                elementary_space(index_types && indices)

    Construct from a list/tuple of indices.

  .. rubric:: Copy/move-constructors, assignments and destructor

  .. function:: elementary_space(elementary_space const&) = default
  .. function:: elementary_space(elementary_space&&) noexcept = default
  .. function:: elementary_space& operator=(elementary_space const&) = default
  .. function:: elementary_space& operator=(elementary_space&&) noexcept \
                = default
  .. function:: virtual ~elementary_space()
  .. function:: virtual std::unique_ptr<elementary_space> clone() const = 0

    Virtual copy-constructor. Makes a copy of this elementary space managed by a
    unique pointer.

  .. rubric:: Algebra ID

  .. function:: virtual int algebra_id() const = 0

    :ref:`Algebra ID <generator>` of the generators associated with this
    elementary space.

  .. rubric:: Index sequence

  .. function:: index_types const& indices() const

    Read-only access to the index tuple carried by this elementary space.

  .. rubric:: Ordering within a direct product

  .. function:: protected virtual bool equal(elementary_space const& es) const
                protected virtual bool less(elementary_space const& es) const
                protected virtual bool greater(elementary_space const& es) const

    These methods can be overridden by the derived classes to establish
    the order of :expr:`es` w.r.t. :expr:`*this` assuming both elementary spaces
    are associated with the same algebra. The default implementation compares
    index tuples of :expr:`*this` and :expr:`es`.

  .. function:: friend bool operator==(generator const& es1, \
                                       generator const& es2)
                friend bool operator!=(generator const& es1, \
                                       generator const& es2)
                friend bool operator<(generator const& es1, \
                                      generator const& es2)
                friend bool operator>(generator const& es1, \
                                      generator const& es2)

    Comparison operators for a pair of elementary spaces. First, they compare
    algebra IDs of :expr:`es1` and :expr:`es2`. If those are equal,
    :expr:`es1.equal(es2)`, :expr:`es1.less(es2)` or :expr:`es1.greater(es2)`
    is called.

  .. rubric:: Binary representation of the basis state index

  .. function:: virtual int n_bits() const = 0

  The number :math:`b` of bits occupied by this elementary space (dimension of
  the space is :math:`2^b`).

.. rubric:: Predefined concrete elementary space types

.. class:: template<typename... IndexTypes> \
           elementary_space_fermion : public elementary_space<IndexTypes...>

  *Defined in <libcommute/loperator/elementary_space_fermion.hpp>*

  An elementary space associated with fermionic algebra generators. This
  elementary space is two-dimensional (:math:`b = 1`).

.. function:: template<typename... IndexTypes> \
              elementary_space_fermion<IndexTypes...> \
              libcommute::static_indices::make_space_fermion( \
              IndexTypes&&... indices)

  *Defined in <libcommute/loperator/elementary_space_fermion.hpp>*

  Make an elementary space associated with fermionic algebra generators
  with given indices.

.. function:: template<typename... IndexTypes> \
              elementary_space_fermion<dyn_indices> \
              libcommute::dynamic_indices::make_space_fermion( \
              IndexTypes&&... indices)

  *Defined in <libcommute/loperator/elementary_space_fermion.hpp>*

  Make an elementary space associated with fermionic algebra generators
  with a given dynamic index sequence.

.. class:: template<typename... IndexTypes> \
           elementary_space_boson : public elementary_space<IndexTypes...>

  *Defined in <libcommute/loperator/elementary_space_boson.hpp>*

  An elementary space associated with bosonic algebra generators. This
  elementary space is truncated and can have an arbitrary dimension of
  form :math:`2^b`.

  .. rubric:: Part of the interface not inherited from / identical to
              :class:`elementary_space`.

  .. function:: template<typename... Args> \
                elementary_space_boson(int n_bits, Args&&... indices)

  Construct a bosonic elementary space of dimension :math:`2^\text{n_bits}`.

.. function:: template<typename... IndexTypes> \
              elementary_space_boson<IndexTypes...> \
              libcommute::static_indices::make_space_boson(int n_bits, \
              IndexTypes&&... indices)

  *Defined in <libcommute/loperator/elementary_space_boson.hpp>*

  Make an elementary space of dimension :math:`2^\text{n_bits}` associated with
  bosonic algebra generators with given indices.

.. function:: template<typename... IndexTypes> \
              elementary_space_boson<dyn_indices> \
              libcommute::dynamic_indices::make_space_boson(int n_bits, \
              IndexTypes&&... indices)

  *Defined in <libcommute/loperator/elementary_space_boson.hpp>*

  Make an elementary space of dimension :math:`2^\text{n_bits}` associated with
  bosonic algebra generators with a given dynamic index sequence.

.. class:: template<typename... IndexTypes> \
           elementary_space_spin : public elementary_space<IndexTypes...>

  *Defined in <libcommute/loperator/elementary_space_spin.hpp>*

  An elementary space associated with spin algebra generators. Dimension of
  this elementary space depends on spin :math:`S`, and is computed as
  :math:`2S+1` rounded up to the nearest power of 2.

  .. rubric:: Part of the interface not inherited from / identical to
              :class:`elementary_space`.

  .. function:: template<typename... Args> \
                elementary_space_spin(double spin, Args&&... indices)

    Construct a spin elementary space with a given spin :math:`S` = :expr:`spin`.

.. function:: template<typename... IndexTypes> \
              elementary_space_spin<IndexTypes...> \
              libcommute::static_indices::make_space_spin(double spin, \
              IndexTypes&&... indices)

  *Defined in <libcommute/loperator/elementary_space_spin.hpp>*

  Make a spin elementary space with :math:`S` = :expr:`spin` and given indices.

.. function:: template<typename... IndexTypes> \
              elementary_space_spin<dyn_indices> \
              libcommute::dynamic_indices::make_space_spin(double spin, \
              IndexTypes&&... indices)

  *Defined in <libcommute/loperator/elementary_space_spin.hpp>*

  Make a spin elementary space with :math:`S` = :expr:`spin` and a given
  dynamic index sequence.

.. _es_constructor:

Advanced: Customization of automatic Hilbert space construction
---------------------------------------------------------------

:func:`make_hilbert_space()` delegates the task of translating algebra
generators into elementary spaces to the functor passed as its second (optional)
argument. It is possible to customize the translation process by giving
:func:`make_hilbert_space()` a callable object similar to the following one

.. code-block:: cpp

  // A custom elementary space constructor object
  struct my_es_constructor {

    template<typename... IndexTypes>
    std::unique_ptr<elementary_space<IndexTypes...>>
    operator()(generator<IndexTypes...> const& g) const {
      //
      // Create an elementary space associated with 'g' and return it
      // wrapped in a unique pointer.
      //
    }

    // Other members if needed ...
  };

This approach gives total control over elementary space creation. It works best
when expressions to be translated do not mix too many algebras and the body
of :expr:`my_es_constructor::operator()` can be kept relatively simple.

Now imagine a different, more common situation, when expressions mix generators
of various predefined algebras as well as generators of a new user-defined
algebra ``my_algebra``. It would be desirable to instruct
:func:`make_hilbert_space()` how to translate instances of
:type:`generator_my_algebra` into :type:`elementary_space_my_algebra` without
rewriting all the code needed to processed the predefined generators. This goal
can be achieved in a few steps by means of a special utility class
:class:`es_constructor`.

.. class:: template<int... AlgebraIDs> es_constructor

  *Defined in <libcommute/loperator/es_constructor.hpp>*

* Define a new algebra ID, e.g. :expr:`my_algebra_id`.

  .. code-block:: cpp

    // A unique integer >=min_user_defined_algebra_id
    static constexpr int my_algebra_id = 7;

* Specialize class :class:`libcommute::es_constructor` as follows

  .. code-block:: cpp

    template<> class es_constructor<my_algebra_id> {
    public:

      es_constructor() = default;

      template<typename... IndexTypes>
      std::unique_ptr<elementary_space<IndexTypes...>>
      operator()(generator<IndexTypes...> const& g) const {
        //
        // Create an elementary space associated with 'g' and return it
        // wrapped in a unique pointer. This method will be called only for
        // the generators of the new algebra, i.e. only when
        // g.algebra_id() == my_algebra_id
        //
      }
    };

  :type:`es_constructor<my_algebra_id>` is obviously a valid elementary space
  constructor for ``my_algebra``

* Instantiate :type:`es_constructor` with multiple template parameters
  (algebra IDs).

  .. code-block:: cpp

    auto es_constr = es_constructor<fermion, spin, my_algebra_id>();

  Now, :expr:`es_constr` knows how to process
  :var:`fermionic <libcommute::fermion>`,
  :var:`spin <libcommute::spin>` and ``my_algebra`` generators.

  .. warning:: The algebra IDs must come in the ascending order when used as
               template parameters of :class:`es_constructor`.

* Finally, call :func:`make_hilbert_space()` with two arguments.

  .. code-block:: cpp

    auto hs = make_hilbert_space(expr, es_constr);

It is worth noting that by default :func:`make_hilbert_space()` uses the
following constructor type as its second argument.

.. type:: default_es_constructor = es_constructor<fermion, spin>

In other words, it recognizes only fermionic and spin generators, and throws
:struct:`es_construction_failure` for all other algebra IDs. If there are
bosonic creation/annihilation operators found in the expression, one may
use another elementary space constructor,

.. type:: boson_es_constructor = es_constructor<fermion, boson, spin>

.. note:: Calling
          ``es_constructor<ID1, ID2, ..., IDN>(arg1, arg2, ..., argK)``
          will internally construct a series of objects
          :expr:`es_constructor<ID1>`, :expr:`es_constructor<ID2>`, ...,
          :expr:`es_constructor<IDN>`. The arguments will be 'fed' into
          constructors of the single-ID objects in order, *at most one
          argument per constructor*. For example,
          :expr:`es_constructor<fermion, boson, spin>(4)` will call
          :expr:`es_constructor<fermion>()`,
          :expr:`es_constructor<boson>(4)` and
          :expr:`es_constructor<spin>()` because the bosonic constructor is the
          first in the sequence accepting one argument.
