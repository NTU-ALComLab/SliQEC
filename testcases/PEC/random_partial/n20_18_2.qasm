OPENQASM 2.0;
include "qelib1.inc";
qreg q[23];
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
t q[3];
s q[1];
h q[6];
t q[3];
t q[2];
t q[5];
t q[3];
s q[1];
cx q[17], q[2];
s q[2];
s q[2];
s q[7];
t q[5];
t q[17];
s q[8];
ccx q[6], q[4], q[1];
ccx q[7], q[2], q[17];
ccx q[9], q[17], q[8];
s q[18];
h q[12];
cx q[12], q[11];
ccx q[1], q[17], q[4];
t q[19];
t q[17];
t q[9];
ccx q[18], q[9], q[7];
ccx q[5], q[4], q[7];
cx q[12], q[3];
ccx q[5], q[2], q[10];
cx q[1], q[10];
h q[9];
cx q[18], q[13];
h q[9];
h q[6];
ccx q[15], q[18], q[14];
ccx q[7], q[3], q[16];
s q[15];
s q[17];
s q[5];
ccx q[14], q[19], q[0];
ccx q[3], q[13], q[0];
h q[1];
h q[11];
cx q[8], q[18];
ccx q[6], q[11], q[17];
ccx q[9], q[13], q[6];
t q[0];
s q[5];
s q[8];
cx q[6], q[7];
ccx q[4], q[2], q[8];
h q[2];
s q[15];
ccx q[5], q[4], q[12];
h q[7];
ccx q[2], q[18], q[11];
cx q[4], q[10];
cx q[16], q[0];
h q[8];
h q[19];
y q[0];
cx q[2], q[1];
tdg q[1];
x q[1];
cx q[2], q[1];
x q[1];
x q[3];
cx q[4], q[3];
x q[3];
sdg q[5];
cx q[8], q[7];
tdg q[7];
z q[9];
s q[11];
ccx q[18], q[14], q[13];
ccx q[11], q[19], q[16];
t q[16];
t q[15];
h q[16];
h q[14];
cx q[12], q[14];
cx q[10], q[14];
ccx q[17], q[12], q[15];
cx q[20], q[6];
cx q[21], q[18];
cx q[22], q[9];