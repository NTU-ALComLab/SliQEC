OPENQASM 2.0;
include "qelib1.inc";
qreg q[47];
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
h q[35];
h q[36];
h q[37];
h q[38];
h q[39];
t q[2];
h q[16];
s q[15];
ccx q[30], q[5], q[8];
cx q[3], q[29];
t q[30];
cx q[0], q[17];
s q[30];
t q[18];
t q[3];
cx q[12], q[26];
cx q[20], q[19];
t q[35];
t q[14];
t q[39];
cx q[0], q[19];
h q[1];
s q[34];
cx q[3], q[7];
s q[13];
t q[35];
ccx q[1], q[21], q[17];
h q[31];
ccx q[28], q[12], q[25];
ccx q[20], q[13], q[10];
ccx q[11], q[15], q[5];
h q[35];
ccx q[6], q[3], q[34];
h q[22];
h q[39];
ccx q[32], q[18], q[27];
h q[29];
ccx q[19], q[21], q[17];
ccx q[38], q[19], q[6];
s q[38];
ccx q[24], q[2], q[30];
ccx q[8], q[12], q[15];
s q[12];
cx q[15], q[3];
ccx q[21], q[22], q[1];
h q[1];
s q[30];
ccx q[26], q[21], q[18];
h q[1];
s q[36];
h q[31];
ccx q[0], q[5], q[38];
ccx q[31], q[16], q[11];
t q[13];
s q[2];
s q[37];
s q[22];
cx q[10], q[37];
s q[38];
h q[19];
h q[32];
cx q[1], q[37];
s q[32];
t q[11];
t q[21];
ccx q[15], q[13], q[28];
ccx q[32], q[15], q[38];
ccx q[38], q[22], q[23];
cx q[30], q[31];
ccx q[2], q[35], q[14];
ccx q[20], q[24], q[19];
cx q[22], q[36];
h q[21];
cx q[9], q[36];
ccx q[29], q[1], q[23];
h q[0];
ccx q[11], q[12], q[39];
h q[24];
s q[4];
s q[11];
ccx q[19], q[2], q[5];
cx q[14], q[10];
h q[23];
cx q[28], q[30];
h q[29];
h q[2];
cx q[2], q[15];
ccx q[26], q[10], q[31];
t q[27];
cx q[2], q[13];
cx q[25], q[26];
ccx q[31], q[15], q[8];
cx q[9], q[24];
h q[12];
ccx q[12], q[18], q[7];
h q[36];
s q[26];
ccx q[29], q[22], q[26];
t q[8];
t q[35];
h q[38];
h q[8];
cx q[14], q[12];
h q[11];
h q[7];
t q[0];
ccx q[28], q[32], q[30];
h q[24];
s q[9];
cx q[5], q[3];
s q[35];
t q[30];
s q[5];
h q[21];
t q[19];
cx q[10], q[23];
t q[22];
cx q[5], q[1];
t q[33];
ccx q[0], q[11], q[26];
s q[30];
cx q[23], q[22];
h q[23];
ccx q[13], q[19], q[36];
t q[0];
sdg q[1];
cx q[0], q[1];
sdg q[0];
x q[2];
cx q[3], q[4];
tdg q[3];
sdg q[6];
sdg q[5];
z q[7];
y q[8];
sdg q[9];
tdg q[9];
y q[11];
cx q[12], q[13];
y q[14];
tdg q[15];
sdg q[16];
cx q[16], q[15];
tdg q[15];
sdg q[17];
s q[39];
t q[27];
t q[26];
h q[24];
t q[25];
ccx q[21], q[37], q[30];
h q[37];
s q[28];
s q[35];
t q[37];
t q[27];
ccx q[26], q[33], q[22];
t q[29];
ccx q[38], q[33], q[20];
t q[28];
cx q[38], q[32];
ccx q[32], q[30], q[23];
t q[22];
t q[36];
h q[32];
cx q[40], q[1];
cx q[41], q[7];
cx q[42], q[11];
cx q[43], q[13];
cx q[44], q[6];
cx q[45], q[35];
cx q[46], q[15];