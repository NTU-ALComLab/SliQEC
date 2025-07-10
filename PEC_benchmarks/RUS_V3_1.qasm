// ===== Description =====
// An repeat-until-success circuit implementing the V3 gate.
// Paired with the care set "RUS_V3.careset".
// =======================

OPENQASM 2.0;
include "qelib1.inc";
qreg q[3];
h q[1];
h q[2];
ccx q[1] q[2] q[0];
s q[0];
ccx q[1] q[2] q[0];
z q[0];
h q[1];
h q[2];
swap q[0] q[1];
swap q[1] q[2];
