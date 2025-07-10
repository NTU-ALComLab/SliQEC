// ===== Description =====
// An repeat-until-success circuit implementing the T gate.
// Paired with the care set "RUS_T.careset".
// =======================

OPENQASM 2.0;
include "qelib1.inc";
qreg q[2];
t q[0]; 
swap q[0], q[1]; 