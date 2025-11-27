.. _monomial:

Monomial
========

.. default-domain:: cpp

Monomials are ordered products of :ref:`algebra generators <generator>`. Even
though the monomials stored as parts of a
:ref:`polynomial expression <expression>` are always canonically ordered,
that need not be true in general.

.. namespace:: libcommute

.. class:: template<typename... IndexTypes> monomial

  Monomial, an ordered product of zero or more algebra generators
  :math:`g_{i_1}, g_{i_2}, \ldots, g_{i_n}`,

  .. math::

    M^{(n)}_{i_1, i_2, \ldots, i_n} = g_{i_1} g_{i_2} \ldots g_{i_n}.

  :type:`IndexTypes` - types of indices carried by the algebra generators.

  .. rubric:: Member type aliases

  .. type:: index_types = std::tuple<IndexTypes...>
  .. type:: generator_type = generator<IndexTypes...>
  .. type:: gen_ptr_type = std::shared_ptr<const generator_type>

  .. rubric:: Constructors

  .. function:: monomial() = default

    Construct an empty (zero-length) monomial :math:`M^{(0)} \equiv 1`.

  .. function:: template<typename... GenTypes> \
                explicit monomial(GenTypes&&... generators)

    Construct a non-empty monomial from a list of generators. This constructor
    calls ``std::make_shared()`` for each of the provided generators, and the
    constructed monomial object takes ownership of the generators.

  .. function:: monomial(std::vector<gen_ptr_type> generators)

    Construct a monomial from a vector of shared pointers to generators. The
    constructed monomial object shares ownership of the generators.

  .. function:: monomial(std::initializer_list<gen_ptr_type> generators)

    Construct a monomial from a list of shared pointers to generators using the
    list initialization syntax. The constructed monomial object shares
    ownership of the generators.

  .. rubric:: Copy/move-constructors and assignments

  .. function:: monomial(monomial const& m)
  .. function:: monomial(monomial&&) noexcept = default
  .. function:: monomial& operator=(monomial const& m)
  .. function:: monomial& operator=(monomial&&) noexcept = default

  .. rubric:: :ref:`Iteration interface <expr_iteration>`

  .. class:: const_iterator

    Random access constant iterator to the list of algebra generators.
    The dereference type is :type:`generator_type`.

  .. function:: const_iterator begin() const noexcept
                const_iterator cbegin() const noexcept

    Constant iterator to the first algebra generator.

  .. function:: const_iterator end() const noexcept
                const_iterator cend() const noexcept

    Constant past-the-end iterator.

  .. type:: const_reverse_iterator = std::reverse_iterator<const_iterator>

    Reverse iteration version of :type:`const_iterator`.

  .. function:: const_reverse_iterator rbegin() const noexcept
                const_reverse_iterator crbegin() const noexcept

    Reverse constant iterator to last algebra generator.

  .. function:: const_reverse_iterator rend() const noexcept
                const_reverse_iterator crend() const noexcept

    Reverse constant past-the-end iterator.

  .. function:: generator_type const& operator[](std::size_t n) const

    Access an algebra generator by its position :var:`n` in the monomial.

  .. rubric:: Other methods and friend functions

  .. function:: std::size_t size() const

    Number of algebra generators in the monomial.

  .. function:: bool empty() const

    Is this a zero-length monomial?

  .. function:: friend bool operator==(monomial const& m1, monomial const& m2)
                friend bool operator!=(monomial const& m1, monomial const& m2)
                friend bool operator<(monomial const& m1, monomial const& m2)
                friend bool operator>(monomial const& m1, monomial const& m2)

    Compare two monomials. The lesser and greater comparisons check monomials'
    lengths first, and in the case of equal lengths perform the
    `lexicographical comparison
    <https://en.cppreference.com/w/cpp/algorithm/lexicographical_compare>`_ of
    the :func:`generator <generator::operator\<>` lists.

  .. function:: bool is_ordered() const

    Is this monomial canonically ordered? In other words, does its algebra
    generator list satisfy :math:`g_{i_1} < g_{i_2} < \ldots < g_{i_n}`?

  .. function:: void swap_generators(std::size_t n1, std::size_t n2)

    Swap algebra generators at positions :var:`n1` and :var:`n2` within the
    list.

  .. function:: void append(generator_type const& g)
                void append(monomial const& m)
                void append(std::pair<const_iterator, const_iterator> const& r)

    Append a generator, a monomial or a range of generators specified by a
    begin-end iterator pair :var:`r` to this monomial. All appended generators
    are expected to be managed by ``std::shared_ptr`` objects, and the monomial
    will take co-ownership of them by `calling \
    <https://en.cppreference.com/w/cpp/memory/enable_shared_from_this/
    shared_from_this.html>`_ ``shared_from_this()``.

  .. function:: template<typename... PartTypes> \
                friend monomial concatenate(PartTypes&&... parts)

    Concatenate a number of parts to build a new monomial. The parts can be
    monomials, algebra generators and ranges within monomials. The ranges
    should be passed as begin-end iterator pairs
    :expr:`std::pair\<const_iterator, const_iterator>`.

  .. function:: friend std::ostream& \
                operator<<(std::ostream& os, monomial const& m)

    Output stream insertion operator.
