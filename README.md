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

## Running the tool

Once the dependencies are installed, navigate to the `src` directory.
Basic usage of the tool involves a command like

```python fairProve.py -f ex.fr -mi -mf```

A full list of command-line flags is available from
`python fairProve.py -h`.

### OOPSLA Benchmarks

You can reproduce our published benchmarks by navigating
to the `oopsla` directory (at the top level).
The benchmarks are divided, based on their fairness postconditions,
into `oopsla/qual` and `oopsla/noqual`.
The `oopsla/run.sh` script can be used to run batches of benchmarks;
see `oopsla/README` for more information.
