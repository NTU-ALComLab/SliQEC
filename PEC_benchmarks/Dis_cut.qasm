// ===== Description =====
// The circuit after cut for distributed quantum computing.
// Using the method in the paper 
//   "Simulating large quantum circuits on a small quantum computer."
// Paired with the weight function "Dis_cut.weight".
// =======================

OPENQASM 2.0;
include "qelib1.inc";
qreg q[7];
swap q[2] q[3];
h q[4];
h q[5];
h q[6];
h q[0];		// U1
cx q[0], q[1];	// U1

x q[5];
mcsdg q[5] q[4] q[1];
x q[5];

rx(pi/2) q[1];
t q[1];
cx q[4] q[1];
cx q[5] q[1];
tdg q[1];
rx(-pi/2) q[1];

x q[6];
cx q[6] q[2];
x q[6];

rx(pi/2) q[2];
t q[2];
cx q[4] q[2];
cx q[5] q[2];
tdg q[2];
rx(-pi/2) q[2];

x q[5];
mcs q[5] q[4] q[2];
x q[5];

s q[2];		// U2
h q[3];		// U2
cx q[3], q[2];	// U2
swap q[0] q[1];
swap q[4] q[3];
swap q[3] q[2];
swap q[2] q[1];
swap q[5] q[4];
swap q[4] q[3];
swap q[3] q[2];
swap q[6] q[5];
swap q[5] q[4];
swap q[4] q[3];