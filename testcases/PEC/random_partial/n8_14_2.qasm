OPENQASM 2.0;
include "qelib1.inc";
qreg q[9];
h q[0];
h q[1];
h q[2];
h q[3];
h q[4];
h q[5];
h q[6];
h q[7];
cx q[7], q[0];
ccx q[7], q[6], q[5];
h q[2];
s q[6];
ccx q[5], q[1], q[7];
t q[2];
s q[7];
cx q[2], q[0];
t q[2];
cx q[1], q[5];
s q[5];
h q[1];
t q[4];
t q[2];
ccx q[3], q[7], q[6];
ccx q[0], q[2], q[5];
cx q[1], q[2];
t q[6];
t q[4];
h q[7];
s q[4];
cx q[0], q[1];
h q[1];
h q[5];
tdg q[0];
x q[2];
cx q[6], q[7];
s q[7];
t q[6];
t q[7];
cx q[8], q[2];