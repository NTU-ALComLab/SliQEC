# SliQEC - A BDD-based Quantum Circuit Equivalence Checker

## Introduction
`SliQEC` is a BDD-based quantum circuit equivalence checker implemented in C/C++ on top of [CUDD](http://web.mit.edu/sage/export/tmp/y/usr/share/doc/polybori/cudd/cuddIntro.html) package. 
The concerned equivalence checking problem includes the full equivalence checking and partial equivalence checking.
For more details about problem formulations and theories, please refer to the [Citation](##Citation).

## Build
First configure CUDD:
```
cd cudd
./configure --enable-dddmp --enable-obj --enable-shared --enable-static
cd ..
```
Then build the checker, type the command at the root directory.
```
$ make
```

## Execution
The circuit format being checked is `OpenQASM` used by IBM's [Qiskit](https://github.com/Qiskit/qiskit), and our supported gate set now contains Pauli-X, Pauli-Y, Pauli-Z, Hadamard, Phase and its inverse, π/8 and its inverse, Rotation-X with phase π/2, Rotation-Y with phase π/2, Controlled-NOT, Controlled-Z, Toffoli, SWAP, and Fredkin. One can find some example benchmarks in [examples](https://github.com/NTU-ALComLab/SliQEC/tree/main/examples) folder.

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
                        0: full 
                        1: strong partial
                        2: unweighted zero-clean-ancilla strong partial
                        3: constant-probability weak partial
  --nQd arg (=0)        (only for --p 1/2/3) #data qubits.
  --nQkc arg (=0)       (only for --p 1/3)   #clean ancilla qubits.
  --nQkd arg (=0)       (only for --p 1/2/3) #dirty ancilla qubits.
  --nQm arg (=0)        (only for --p 1/2/3) #measured qubits.
  --nQw arg (=0)        (only for --p 1/3)   #weighted qubits.
  --nQp arg (=0)        (only for --p 1/2/3) #propagating qubits.
  --nQg arg (=0)        (only for --p 1/2/3) #garbage qubits.
  --nQkr arg (=0)       (only for --p 1/3)   #reverted clean ancilla qubits.
  --careSet arg         (only for --p 1/2/3) the care set of the circuits.
                        Can be omitted if the care set is the universal set.
  --weightFun1 arg      (only for --p 1/2/3 and nQw > 0)the weight function of
                        the 1st circuit
  --weightFun2 arg      (only for --p 1/2/3 and nQw > 0)the weight function of
                        the 2nd circuit
```

For partial equivalence checking, please note that we assume the input qubits follow the order of data qubits -> clean ancilla qubits -> dirty ancilla qubits, and the output qubits follow the order of measured qubits -> weighted qubits -> propagating qubits -> garbage qubits -> reverted clean ancilla qubits -> dirty ancilla qubits. 
This assumption is made without loss of generality, as swap gates can be used to adjust their positions if needed.

## Example
#### Full Equivalence Checking
For conducting full equivalence checking on [examples/FEC/bv_1.qasm](https://github.com/NTU-ALComLab/SliQEC/blob/main/examples/FEC/bv_1.qasm) and [examples/FEC/bv_2.qasm](https://github.com/NTU-ALComLab/SliQEC/blob/main/examples/FEC/bv_2.qasm), execute:
``` commandline
$ ./SliQEC --circuit1 examples/FEC/bv_1.qasm --circuit2 examples/FEC/bv_2.qasm
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
$ ./SliQEC --p 1 --circuit1 examples/PEC/period_finding_1.qasm --circuit2 examples/PEC/period_finding_2.qasm --nQd 3 --nQm 3 --nQkc 5 --nQg 5
```
Then the results will be shown:
``` 
{
        #Qubits (n): 8
        #Data qubits (n_d):                    3
        #Clean ancilla qubits (n_kc):          5
        #Dirty ancilla qubits (n_kd):          0
        #Measured qubits (n_m):                3
        #Weighted qubits (n_w):                0
        #Propagating qubits (n_p):             0
        #Garbage qubits (n_g):                 5
        #Reverted clean ancilla qubits (n_kr): 0
        Gatecount of circuit1: 395
        Gatecount of circuit2: 437
        Is strongly partially equivalent? Yes
}

Runtime: 0.638659 seconds
Peak memory usage: 14987264 bytes
```

## Citation
<summary>
    <a href="https://doi.org/10.1145/3489517.3530481">C.-Y. Wei, Y.-H. Tsai, C. -S. Jhang, and J.-H. R. Jiang, “Accurate BDD-based Unitary Operator Manipulation for Scalable and Robust Quantum Circuit Verification,” in Proceedings of the <em>Design Automation Conference (DAC)</em>, 2022, pp. 523-528. </a>
</summary>
<summary>
    <a href="https://doi.org/10.1109/QCE53715.2022.00082">T.-F. Chen, J.-H. R. Jiang, and M.-H. Hsieh, “Partial Equivalence Checking of Quantum Circuits,” in Proceedings of the <em>International Conference on Quantum Computing and Engineering (QCE)</em>, 2022, pp. 594-604. </a>
</summary>

## Contact
If you have any questions or suggestions, feel free to [create an issue](https://github.com/NTU-ALComLab/SliQEC/issues), or contact us.
For full equivalence checking, please contact joey.cywei@gmail.com.
For partial equivalence checking, please contact  ghdftff542@gmail.com.
