# ![logo](logo/logo.svg) FairSquare

A tool to verify algorithmic fairness properties, as seen in the OOPSLA 2017
paper **Probabilistic Verification of Program Fairness** by
[Aws Albarghouthi](http://pages.cs.wisc.edu/~aws/),
[Loris D'Antoni](http://pages.cs.wisc.edu/~loris/),
[Samuel Drews](http://pages.cs.wisc.edu/~sdrews/), and
[Aditya Nori](https://www.microsoft.com/en-us/research/people/adityan/).

## Installation

FairSquare runs on Linux and uses Python 3.
Additionally, it uses the following dependencies:
[z3](http://github.com/Z3Prover/z3),
the [SciPy](http://scipy.org/) stack,
and the python packages `codegen` and `asteval`
(easily obtained using [pip](http://pypi.python.org/pypi/pip):
`pip install --user codegen asteval`)

By default, FairSqaure uses
[redlog](http://www.redlog.eu/get-redlog/)
for quantifier elimination, which must be downloaded separately:
FairSquare expects the downloaded files to be placed in `src/tools/`
(specifically, it runs `src/tools/reduce`,
which relies on the existence of `src/tools/reduce.img`).
Alternatively, FairSquare can use a DNF-based quantifier elimination
implemented with z3 by using the `-z` flag when running the tool.

## Running the tool

Once the dependencies are installed, navigate to the `src` directory.
Basic usage of the tool involves a command like

```python fairProve.py -f ex.fr -mi -mf```

`-f ex.fr` specifies that `ex.fr` contains the fairness verification problem to be analyzed
(in fact, `src/ex.fr` contains the illustrative example from the OOPSLA paper),
and `-mi -mf` are optimizations that should always be used.
A full list of command-line flags is available from
`python fairProve.py -h`.

### OOPSLA Benchmarks

You can reproduce our published benchmarks by navigating
to the `oopsla` directory (at the top level).
The benchmarks are divided, based on their fairness postconditions,
into `oopsla/qual` and `oopsla/noqual`.
The `oopsla/run.sh` script can be used to run batches of benchmarks;
see `oopsla/README` for more information.

## Contributors

* [Samuel Drews](http://pages.cs.wisc.edu/~sdrews/) (primary contact)
* [David Merrell](https://dpmerrell.github.io/about/)
* [Aws Albarghouthi](http://pages.cs.wisc.edu/~aws/)
* [Loris D'Antoni](http://pages.cs.wisc.edu/~loris/)
* [Aditya Nori](https://www.microsoft.com/en-us/research/people/adityan/)
