OPENQASM 2.0;
include "qelib1.inc";
qreg q[35];
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
h q[30];
h q[31];
h q[32];
h q[33];
h q[34];
s q[20];
ccx q[11], q[18], q[10];
ccx q[26], q[34], q[2];
t q[7];
cx q[8], q[6];
ccx q[2], q[0], q[4];
h q[22];
h q[34];
t q[9];
cx q[13], q[6];
h q[32];
cx q[29], q[2];
t q[18];
cx q[12], q[33];
s q[23];
cx q[33], q[4];
h q[15];
cx q[33], q[23];
t q[3];
t q[25];
h q[24];
s q[5];
s q[17];
ccx q[16], q[15], q[25];
t q[2];
h q[18];
cx q[29], q[3];
cx q[22], q[31];
h q[31];
ccx q[12], q[16], q[28];
s q[9];
h q[25];
h q[30];
ccx q[30], q[29], q[7];
t q[5];
cx q[17], q[12];
t q[5];
cx q[26], q[32];
h q[27];
s q[26];
t q[5];
s q[1];
cx q[6], q[15];
t q[10];
cx q[17], q[9];
h q[10];
ccx q[4], q[24], q[30];
s q[10];
ccx q[17], q[32], q[22];
cx q[11], q[13];
cx q[8], q[15];
cx q[21], q[33];
cx q[11], q[12];
cx q[19], q[4];
ccx q[5], q[11], q[18];
s q[28];
t q[0];
t q[23];
h q[2];
cx q[6], q[0];
h q[1];
h q[10];
h q[3];
cx q[28], q[25];
s q[23];
cx q[3], q[4];
cx q[11], q[10];
cx q[2], q[27];
h q[22];
s q[34];
t q[0];
cx q[22], q[27];
cx q[18], q[9];
t q[28];
ccx q[4], q[20], q[25];
s q[10];
t q[0];
t q[5];
h q[30];
ccx q[14], q[4], q[10];
s q[6];
s q[1];
s q[12];
cx q[24], q[10];
ccx q[15], q[11], q[17];
ccx q[33], q[34], q[15];
h q[8];
s q[15];
t q[3];
h q[9];
ccx q[29], q[18], q[13];
ccx q[13], q[33], q[8];
cx q[5], q[31];
ccx q[22], q[32], q[10];
h q[26];
t q[30];
h q[33];
cx q[25], q[21];
h q[10];
cx q[6], q[3];
ccx q[27], q[12], q[32];
s q[3];
h q[17];
t q[23];
t q[14];
sdg q[1];
cx q[0], q[1];
tdg q[1];
cx q[0], q[1];
sdg q[0];
sdg q[2];
cx q[3], q[2];
tdg q[3];
sdg q[2];
x q[6];
sdg q[7];
cx q[8], q[7];
sdg q[7];
cx q[8], q[7];
tdg q[7];
tdg q[9];
cx q[11], q[12];
tdg q[11];
sdg q[12];
cx q[11], q[12];
sdg q[11];
sdg q[13];
cx q[14], q[13];
x q[13];
y q[15];
s q[20];
s q[33];
ccx q[34], q[23], q[18];
cx q[23], q[32];
s q[34];
ccx q[20], q[29], q[30];
ccx q[31], q[18], q[24];
s q[22];
t q[19];
s q[27];
t q[19];
cx q[22], q[31];
s q[29];
h q[29];
s q[22];
s q[23];
t q[27];
ccx q[33], q[32], q[24];