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
h q[8];
ccx q[5], q[3], q[4];
ccx q[6], q[5], q[1];
ccx q[8], q[7], q[4];
h q[6];
cx q[4], q[5];
cx q[5], q[6];
cx q[1], q[2];
cx q[7], q[8];
s q[3];
s q[1];
t q[2];
ccx q[8], q[2], q[4];
ccx q[7], q[3], q[1];
ccx q[6], q[7], q[5];
s q[0];
ccx q[8], q[4], q[2];
ccx q[1], q[0], q[3];
t q[6];
h q[4];
cx q[0], q[7];
cx q[8], q[5];
cx q[6], q[4];
ccx q[0], q[1], q[8];
s q[3];
s q[4];
t q[8];
s q[8];
cx q[1], q[0];
sdg q[0];
cx q[1], q[0];
sdg q[0];
tdg q[0];
z q[2];
z q[3];
cx q[7], q[4];
s q[6];
ccx q[4], q[6], q[8];
cx q[7], q[8];
cx q[6], q[5];