TODO list
=========

* generator: add `virtual int nilpotent_power() const = 0;`
* generator: add `virtual int idempotent_power() const = 0;` (???)
* expression: add aliases

    template<typename... IndexTypes> using expr_real = expression<double, IndexTypes...>;

    template<typename... IndexTypes> using expr_complex = expression<std::complex<double>, IndexTypes...>;

    // C++17 only

    template<typename... IndexTypes> using dyn_expr_real = expression<double, dyn_indices>;

    template<typename... IndexTypes> using dyn_expr_complex = expression<std::complex<double>, dyn_indices>;
