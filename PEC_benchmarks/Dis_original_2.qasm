// ===== Description =====
// The original circuit for distributed quantum computing.
// Paired with the weight function "Dis_original.weight", 
//   which is empty because the original circuit is unweighted.
// =======================

OPENQASM 2.0;
include "qelib1.inc";
qreg q[8];
h q[0];
cx q[0], q[1];
// The circuit is cut here
s q[1];
h q[2];
cx q[2], q[1];
swap q[2] q[3];
swap q[3] q[4];
swap q[4] q[5];
swap q[5] q[6];
swap q[6] q[7];
swap q[1] q[2];
swap q[2] q[3];
swap q[3] q[4];
swap q[4] q[5];
swap q[5] q[6];
swap q[0] q[1];
swap q[1] q[2];
swap q[2] q[3];
swap q[3] q[4];
swap q[4] q[5];