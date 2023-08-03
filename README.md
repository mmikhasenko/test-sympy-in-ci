# SymPy in Pluto

The repo test the sympy installation with Pluto notebooks in CI.
It used to be sufficient to add empty `PYTHON` env. variables to CI yaml file. 
```
        runs-on: ubuntu-latest
        env:
            PYTHON: ''
```

Let's see if it works.

For testing:

- [`html link`](crosscheck/N-001-cartesian-pions.html)
