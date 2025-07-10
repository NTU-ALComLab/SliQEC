// ===== Description =====
// A wrongly implemented 3-controlled Toffoli gate with q[3] being the target bit.
// Trying to use q[4] as a reverted clean ancilla qubit,
//   but q[4] is not returned to |0>.
// =======================

OPENQASM 2.0;
include "qelib1.inc";
qreg q[5];
mcx q[0] q[1] q[2] q[4];
cx q[4] q[3];

