// ===== Description =====
// An repeat-until-success circuit implementing the T gate.
// Paired with the care set "RUS_T.careset".
// =======================

OPENQASM 2.0;
include "qelib1.inc";
qreg q[2];
cx q[0], q[1];
t q[1];
h q[1];
swap q[0], q[1];