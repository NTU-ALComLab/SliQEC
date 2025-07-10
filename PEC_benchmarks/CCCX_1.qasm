// ===== Description =====
// A 3-controlled Toffoli gate with q[3] being the target bit.
// =======================

OPENQASM 2.0;
include "qelib1.inc";
qreg q[5];
mcx q[0] q[1] q[2] q[3];


