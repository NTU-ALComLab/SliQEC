# SliQEC - A BDD-based Quantum Circuit Equivalence Checker

## Introduction
`SliQEC` is a BDD-based quantum circuit equivalence checker implemented in C/C++ on top of [CUDD](http://web.mit.edu/sage/export/tmp/y/usr/share/doc/polybori/cudd/cuddIntro.html) package. 
The concerned equivalence checking problem includes the full equivalence checking and partial equivalence checking.
For more details about problem formulations and theories, please refer to the [papers](##Citation).

## Build
To build the checker, type the command at the root directory.
```
$ make
```

## Execution
The circuit format being checked is `OpenQASM` used by IBM's [Qiskit](https://github.com/Qiskit/qiskit), and our supported gate set now contains Pauli-X, Pauli-Y, Pauli-Z, Hadamard, Phase and its inverse, π/8 and its inverse, Rotation-X with phase π/2, Rotation-Y with phase π/2, Controlled-NOT, Controlled-Z, Toffoli, SWAP, and Fredkin.

The help message concludes the details for execution:

``` 
$ ./SliQEC --help
Options:
  --help                produce help message.
  --reorder arg (=1)    allow variable reordering or not.
                        0: disable 1: enable
  --circuit1 arg        1st circuit for equivalence checking.
  --circuit2 arg        2nd circuit for equivalence checking.
  --p arg (=0)          conduct full or partial equivalence checking.
                        0: full 1: partial
  --nQd arg (=0)        (only for --p 1) the number of data qubits.
  --nQm arg (=0)        (only for --p 1) the number of measured qubits.
```

## Example
#### Full Equivalence Checking
For conducting full equivalence checking on [examples/FEC/bv_1.qasm](https://github.com/NTU-ALComLab/SliQEC/blob/main/examples/FEC/bv_1.qasm) and [examples/FEC/bv_2.qasm](https://github.com/NTU-ALComLab/SliQEC/blob/main/examples/FEC/bv_2.qasm), execute:
``` commandline
./SliQEC --circuit1 examples/FEC/bv_1.qasm --circuit2 examples/FEC/bv_2.qasm
```
Then the results will be shown:
``` 
{
	#Qubits (n): 10
	Gatecount of circuit1: 29
	Gatecount of circuit2: 62
	Is equivalent? Yes
}

Runtime: 0.019328 seconds
Peak memory usage: 12881920 bytes
```

#### Partial Equivalence Checking
For conducting partial equivalence checking on [examples/PEC/period_finding_1.qasm](https://github.com/NTU-ALComLab/SliQEC/blob/main/examples/PEC/period_finding_1.qasm) and [examples/PEC/period_finding_2.qasm](https://github.com/NTU-ALComLab/SliQEC/blob/main/examples/PEC/period_finding_2.qasm) with 3 data qubits and 3 measured qubits, execute:
``` commandline
./SliQEC --p 1 --circuit1 examples/PEC/period_finding_1.qasm --circuit2 examples/PEC/period_finding_2.qasm --nQd 3 --nQm 3
```
Then the results will be shown:
``` 
{
        #Qubits (n): 8
        #Data qubits (d): 3
        #Measured qubits (m): 3
        Gatecount of circuit1: 395
        Gatecount of circuit2: 437
        Is partially equivalent? Yes
}

Runtime: 1.86163 seconds
Peak memory usage: 14553088 bytes
```

## Citation

## Contact
If you have any questions or suggestions, feel free to [create an issue](https://github.com/NTU-ALComLab/SliQEC/issues), or contact us.
For full equivalence checking, please contact joey.cywei@gmail.com.
For partial equivalence checking, please contact  ghdftff542@gmail.com.
