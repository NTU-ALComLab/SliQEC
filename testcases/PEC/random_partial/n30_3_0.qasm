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
h q[25];
h q[26];
h q[27];
h q[28];
h q[29];
s q[28];
ccx q[27], q[2], q[11];
cx q[27], q[23];
h q[26];
t q[4];
s q[5];
s q[18];
ccx q[23], q[0], q[8];
h q[1];
ccx q[5], q[14], q[25];
t q[23];
t q[12];
t q[28];
h q[21];
s q[6];
h q[5];
s q[26];
s q[1];
s q[19];
ccx q[11], q[18], q[24];
ccx q[10], q[6], q[25];
h q[3];
h q[4];
h q[27];
ccx q[27], q[21], q[3];
h q[13];
cx q[6], q[17];
ccx q[3], q[4], q[29];
ccx q[27], q[9], q[10];
ccx q[18], q[7], q[27];
s q[13];
h q[15];
s q[3];
cx q[23], q[8];
ccx q[8], q[11], q[13];
t q[4];
ccx q[26], q[28], q[19];
cx q[5], q[25];
h q[17];
h q[17];
cx q[28], q[11];
cx q[25], q[9];
h q[0];
s q[5];
h q[13];
cx q[0], q[2];
ccx q[21], q[28], q[25];
h q[2];
ccx q[4], q[21], q[10];
ccx q[8], q[9], q[12];
t q[17];
ccx q[25], q[10], q[11];
s q[10];
t q[26];
h q[24];
t q[9];
h q[25];
cx q[28], q[18];
ccx q[10], q[8], q[15];
s q[22];
s q[27];
ccx q[0], q[19], q[5];
h q[29];
t q[17];
t q[0];
s q[10];
h q[28];
t q[29];
cx q[24], q[19];
ccx q[29], q[20], q[11];
ccx q[19], q[15], q[10];
t q[6];
h q[18];
cx q[8], q[2];
t q[16];
h q[0];
cx q[14], q[6];
t q[23];
cx q[23], q[4];
t q[14];
h q[7];
h q[16];
h q[19];
s q[19];
ccx q[25], q[8], q[4];
cx q[22], q[27];
t q[0];
s q[27];
cx q[19], q[27];
cx q[14], q[2];
sdg q[0];
z q[2];
sdg q[3];
tdg q[4];
tdg q[3];
cx q[3], q[4];
tdg q[6];
sdg q[5];
cx q[5], q[6];
sdg q[5];
cx q[7], q[8];
sdg q[8];
cx q[7], q[8];
sdg q[7];
x q[9];
tdg q[10];
y q[12];
z q[13];
h q[28];
ccx q[19], q[23], q[17];
ccx q[23], q[26], q[24];
s q[23];
s q[23];
ccx q[25], q[29], q[27];
h q[18];
h q[24];
h q[19];
ccx q[17], q[25], q[27];
h q[26];
h q[16];
s q[18];
t q[27];
cx q[26], q[21];