OPENQASM 2.0;
include "qelib1.inc";
qreg q[10];
h q[0];
h q[1];
h q[2];
h q[3];
h q[4];
h q[5];
h q[6];
h q[7];
h q[8];
h q[9];
ccx q[1], q[3], q[5];
t q[0];
s q[7];
cx q[6], q[5];
ccx q[2], q[7], q[9];
h q[8];
ccx q[4], q[3], q[5];
t q[4];
h q[6];
ccx q[1], q[8], q[6];
s q[4];
s q[2];
s q[2];
h q[3];
t q[4];
h q[8];
h q[6];
t q[9];
ccx q[5], q[1], q[8];
cx q[1], q[8];
s q[0];
cx q[9], q[3];
s q[6];
cx q[4], q[6];
cx q[7], q[1];
s q[4];
h q[1];
cx q[7], q[4];
cx q[9], q[5];
ccx q[3], q[4], q[2];
sdg q[0];
x q[0];
cx q[1], q[0];
x q[0];
y q[2];
x q[3];
ccx q[7], q[9], q[8];
ccx q[9], q[8], q[6];
t q[6];
h q[9];
ccx q[7], q[5], q[8];