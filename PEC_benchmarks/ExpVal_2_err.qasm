// ===== Description =====
// Calaulating the expectation value of 
//   ZZ-basis measurement on q[0], q[1] after the circuit.
// Trying to use q[2] as a dirty ancilla qubit,
//   but q[2] is not always returned to its initial state 
//   due to wrong usage of the bridge gate.
// Paired with the weight function "ExpVal.weight".
// The final swap gate and Z gate should not affect the expectation value.
// =======================

OPENQASM 2.0;
include "qelib1.inc";
qreg q[3];
h q[0];
h q[1];
s q[0];
t q[1];
h q[0];
h q[1];
s q[0];
t q[0];
tdg q[1];
h q[0];
h q[1];
cx q[0] q[2];
cx q[2] q[1];
cx q[0] q[2];