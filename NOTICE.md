This repository contains two code boundaries:

- Project-local code under `src/abag_rbfe` and `src/abag_pmx`, licensed under
  the top-level MIT license declared in `pyproject.toml`.
- Vendored upstream `pmx` source under `vendor/pmx`, preserved under its
  original upstream license and notices.

The vendored pmx tree is included as a standalone scientific core and is not
rewritten into the local package namespace.

