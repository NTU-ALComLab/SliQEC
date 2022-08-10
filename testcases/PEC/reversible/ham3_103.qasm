OPENQASM 2.0;
include "qelib1.inc";
qreg q[3];
mcx q[1], q[2], q[0];
cx q[0], q[1];
swap q[1], q[2];
cx q[0], q[1];
