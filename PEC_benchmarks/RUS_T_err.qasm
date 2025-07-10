// ===== Description =====
// An incorrect repeat-until-success circuit trying to implement the T gate.
// Paired with the care set "RUS_T.careset".
// =======================

OPENQASM 2.0;
include "qelib1.inc";
qreg q[2];
cx q[0], q[1];
tdg q[1];
h q[1];
swap q[0], q[1];
