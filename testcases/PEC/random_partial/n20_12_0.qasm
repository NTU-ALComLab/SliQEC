OPENQASM 2.0;
include "qelib1.inc";
qreg q[20];
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
h q[10];
h q[11];
h q[12];
h q[13];
h q[14];
h q[15];
h q[16];
h q[17];
h q[18];
h q[19];
cx q[15], q[3];
h q[8];
t q[17];
ccx q[9], q[16], q[18];
t q[0];
s q[6];
s q[7];
h q[10];
h q[8];
cx q[10], q[17];
s q[5];
t q[5];
s q[18];
s q[6];
h q[7];
s q[13];
s q[1];
cx q[16], q[7];
ccx q[5], q[17], q[12];
s q[6];
ccx q[14], q[17], q[6];
ccx q[14], q[10], q[6];
s q[5];
ccx q[7], q[2], q[8];
t q[14];
t q[8];
h q[1];
h q[7];
h q[2];
t q[19];
ccx q[15], q[11], q[10];
s q[2];
cx q[11], q[5];
ccx q[7], q[14], q[6];
h q[5];
ccx q[6], q[18], q[0];
ccx q[6], q[16], q[9];
s q[9];
s q[16];
cx q[6], q[17];
ccx q[7], q[11], q[1];
t q[16];
ccx q[8], q[5], q[13];
s q[10];
t q[2];
h q[9];
ccx q[18], q[3], q[10];
h q[0];
cx q[1], q[15];
h q[9];
ccx q[8], q[4], q[10];
s q[1];
cx q[13], q[2];
s q[8];
ccx q[9], q[7], q[2];
ccx q[17], q[4], q[19];
s q[14];
s q[18];
h q[13];
h q[0];
z q[0];
y q[1];
cx q[3], q[2];
sdg q[2];
tdg q[5];
cx q[4], q[5];
sdg q[5];
sdg q[4];
sdg q[7];
ccx q[14], q[10], q[18];
cx q[12], q[10];
ccx q[19], q[17], q[15];
cx q[15], q[12];
cx q[16], q[10];
h q[12];
ccx q[15], q[12], q[13];
cx q[11], q[15];
h q[17];
cx q[19], q[17];