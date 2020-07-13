![libcommute logo](doc/logo.svg)

[![Build Status](https://travis-ci.org/krivenko/libcommute.svg?branch=master)](
https://travis-ci.org/krivenko/libcommute)
[![Documentation](https://img.shields.io/badge/docs-GitHub%20Pages-red)](
https://krivenko.github.io/libcommute)

*libcommute* is a C++11/14/17 template library including two major parts.

* A Domain-Specific Language (DSL) designed to easily construct and manipulate
  polynomial expressions with quantum-mechanical operators,
  especially those used in the quantum many-body theory.

* A fast intermediate representation of the quantum-mechanical operators
  that enables their action on finite-dimensional state vectors.
  This feature provides a basis for writing highly performant Exact
  Diagonalization (ED) codes without loss of flexibility.

Dependencies
------------

A C++11 conformant compiler. C++17 support is required for
the dynamic index sequence feature.

Installation
------------

TODO

Usage
-----

Once *libcommute* is installed, you can use it in your CMake project. Here is
a minimal example of an application `CMakeLists.txt` file.

```cmake
# TODO
```

Here is how `test.cpp` could look like.
```c++
// TODO
```

Citing
------

If you find this library useful for your research, you can help me by citing it
using the following BibTeX entry.

```
TODO: Zenodo BibTeX entry here
```

License
-------

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
