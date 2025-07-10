// ===== Description =====
// A 3-controlled Toffoli gate with q[3] being the target bit.
// Using q[4] as a dirty ancilla qubit.
// =======================

OPENQASM 2.0;
include "qelib1.inc";
qreg q[5];
cx q[4] q[3];
mcx q[0] q[1] q[2] q[4];
cx q[4] q[3];
mcx q[0] q[1] q[2] q[4];