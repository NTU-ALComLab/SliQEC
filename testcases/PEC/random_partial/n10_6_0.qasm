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
s q[6];
s q[7];
h q[6];
s q[3];
ccx q[8], q[1], q[6];
h q[8];
ccx q[0], q[2], q[7];
t q[4];
ccx q[8], q[7], q[1];
ccx q[0], q[9], q[3];
ccx q[5], q[0], q[2];
ccx q[9], q[7], q[1];
ccx q[5], q[9], q[8];
t q[4];
ccx q[5], q[7], q[3];
cx q[0], q[5];
t q[6];
h q[3];
s q[7];
ccx q[7], q[4], q[5];
ccx q[2], q[0], q[3];
t q[5];
h q[7];
s q[6];
ccx q[1], q[3], q[9];
ccx q[2], q[1], q[4];
cx q[3], q[8];
ccx q[1], q[5], q[9];
h q[9];
s q[2];
y q[0];
y q[1];
x q[2];
tdg q[3];
cx q[3], q[4];
tdg q[3];
cx q[9], q[6];
ccx q[7], q[9], q[5];
cx q[8], q[5];
s q[7];
s q[7];