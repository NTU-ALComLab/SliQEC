# SliQEQ - A BDD-based Quantum Circuit Equivalence Checker

## Introduction
`SliQEC` is a BDD-based quantum circuit equivalence checker implemented in C/C++ on top of [CUDD](http://web.mit.edu/sage/export/tmp/y/usr/share/doc/polybori/cudd/cuddIntro.html) package. 
The circuit format being checked is `OpenQASM` used by IBM's [Qiskit](https://github.com/Qiskit/qiskit), and our gate set supported now contains Pauli-X, Pauli-Y, Pauli-Z, Hadamard, Phase and its inverse, π/8 and its inverse, Rotation-X with phase π/2, Rotation-Y with phase π/2, Controlled-NOT, Controlled-Z, Toffoli, SWAP, and Fredkin.

## Build
To build the checker, type 
```commandline
make
```
at the root directory.

## Execution
The help message states the details:

```commandline
$ ./SliQEC --help
Options:
  --help                produce help message.
  --print_info          print statistics such as runtime, memory, etc.
  --r arg (=32)         integer bit size.
  --reorder arg (=1)    allow variable reordering or not.
                        0: disable reordering.
                        1: enable reordering (default option).
  --p                   toggle conducting partial equivalence checking.
  --circuit1 arg        1st circuit for equivalence checking.
  --circuit2 arg        2nd circuit for equivalence checking.
  --nQin arg (=0)       the number of input qubits.
  --nQout arg (=0)      the number of output qubits.
  --s                   toggle using algorithm2 for partial equivalence 
                        checking.
```
