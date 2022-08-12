# SliQEQ - A BDD-based Quantum Circuit Equivalence Checker

## Introduction
`SliQEC` is a BDD-based quantum circuit equivalence checker implemented in C/C++ on top of [CUDD](http://web.mit.edu/sage/export/tmp/y/usr/share/doc/polybori/cudd/cuddIntro.html) package. In `SliQSim`, a bit-slicing technique based on BDDs is used to represent quantum state vectors. For more details of the simulator, please refer to the [paper](https://arxiv.org/abs/2007.09304).
The circuit format being checked is `OpenQASM` used by IBM's [Qiskit](https://github.com/Qiskit/qiskit), and our gate set supported now contains Pauli-X (x), Pauli-Y (y), Pauli-Z (z), Hadamard (h), Phase and its inverse (s and sdg), π/8 and its inverse (t and tdg), Rotation-X with phase π/2 (rx(pi/2)), Rotation-Y with phase π/2 (ry(pi/2)), Controlled-NOT (cx), Controlled-Z (cz), Toffoli (ccx and mcx), SWAP (swap), and Fredkin (cswap).

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
