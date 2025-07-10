// ===== Description =====
// Calaulating the expectation value of 
//   ZZ-basis measurement on q[0], q[1] after the circuit.
// Using q[2] as a dirty ancilla qubit.
// Paired with the weight function "ExpVal.weight".
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
cx q[2] q[1];
cx q[0] q[2];
cx q[2] q[1];
cx q[0] q[2];