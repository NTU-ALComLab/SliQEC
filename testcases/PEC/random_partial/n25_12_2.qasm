OPENQASM 2.0;
include "qelib1.inc";
qreg q[30];
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
h q[20];
h q[21];
h q[22];
h q[23];
h q[24];
t q[12];
h q[17];
h q[7];
t q[6];
t q[16];
ccx q[21], q[10], q[3];
cx q[9], q[20];
h q[10];
h q[3];
ccx q[12], q[1], q[8];
cx q[12], q[24];
s q[15];
ccx q[18], q[0], q[20];
s q[4];
cx q[5], q[3];
cx q[22], q[6];
h q[10];
cx q[9], q[10];
t q[22];
ccx q[14], q[0], q[24];
h q[11];
ccx q[2], q[18], q[3];
h q[4];
s q[9];
s q[8];
h q[9];
t q[24];
s q[14];
h q[3];
s q[3];
t q[5];
t q[12];
h q[2];
h q[4];
t q[1];
cx q[24], q[1];
h q[20];
h q[23];
t q[9];
h q[11];
s q[2];
cx q[13], q[10];
s q[3];
cx q[19], q[7];
h q[18];
h q[8];
h q[3];
t q[11];
s q[20];
t q[4];
ccx q[12], q[10], q[14];
cx q[22], q[13];
t q[1];
cx q[11], q[1];
t q[14];
cx q[3], q[17];
h q[3];
cx q[11], q[20];
ccx q[11], q[9], q[23];
ccx q[12], q[16], q[9];
cx q[20], q[19];
s q[14];
cx q[23], q[2];
s q[18];
cx q[0], q[24];
cx q[11], q[22];
t q[10];
s q[20];
s q[8];
t q[23];
s q[0];
t q[19];
h q[12];
ccx q[13], q[11], q[1];
t q[11];
cx q[0], q[1];
sdg q[0];
sdg q[3];
sdg q[2];
x q[4];
tdg q[6];
cx q[6], q[5];
tdg q[5];
cx q[6], q[5];
sdg q[5];
tdg q[7];
cx q[7], q[8];
cx q[10], q[9];
tdg q[9];
t q[24];
t q[16];
ccx q[14], q[23], q[24];
s q[15];
s q[12];
ccx q[15], q[16], q[24];
t q[19];
h q[13];
t q[22];
ccx q[18], q[21], q[13];
ccx q[21], q[23], q[22];
ccx q[16], q[12], q[22];
ccx q[16], q[21], q[13];
cx q[25], q[16];
cx q[26], q[13];
cx q[27], q[14];
cx q[28], q[4];
cx q[29], q[22];