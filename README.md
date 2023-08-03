# SymPy in Pluto

The repo test the `Sympy` installation with Pluto notebooks in CI.
<img width="862" alt="image" src="https://github.com/mmikhasenko/test-sympy-in-ci/assets/22725744/98e859d6-f291-45a6-887f-9196db8e21a9">

It used to be sufficient to add empty the `PYTHON` environment variable to [CI yaml file](.github/pages.yaml). 
```
        runs-on: ubuntu-latest
        env:
            PYTHON: ''
```

However, the notebooks appear to fail.
E.g. https://mmikhasenko.github.io/rescattering-3pi-2023-001/crosscheck/N-001-cartesian-pions.html

Using this repository, I found out that the probkem is caused by another `Project.toml` file in a subfolder.
Naively, it should not affer the top level file processing, but it does.

Here a comparison of CI cout for running [before](https://github.com/mmikhasenko/test-sympy-in-ci/actions/runs/5749591291/job/15584695559) and [after](https://github.com/mmikhasenko/test-sympy-in-ci/actions/runs/5749764081/job/15585199420) delecting the `KTMC.jl/Project.toml`
<img width="1883" alt="image" src="https://github.com/mmikhasenko/test-sympy-in-ci/assets/22725744/d6c7b7d2-7b44-4cdc-96a2-a67fc40d45c8">


For testing, once you are on the [Pages](https://mmikhasenko.github.io/test-sympy-in-ci/):

- [`html link`](N-001-cartesian-pions.html)
