TODO list
=========

* generator: add `virtual int nilpotent_power() const = 0;`
* generator: add `virtual int idempotent_power() const = 0;` (???)
* expression: add aliases

    // C++17 only

    template<typename... IndexTypes> using dyn_expr_real = expression<double, dyn_indices>;

    template<typename... IndexTypes> using dyn_expr_complex = expression<std::complex<double>, dyn_indices>;
