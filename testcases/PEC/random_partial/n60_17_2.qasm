OPENQASM 2.0;
include "qelib1.inc";
qreg q[70];
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
h q[40];
h q[41];
h q[42];
h q[43];
h q[44];
h q[45];
h q[46];
h q[47];
h q[48];
h q[49];
h q[50];
h q[51];
h q[52];
h q[53];
h q[54];
h q[55];
h q[56];
h q[57];
h q[58];
h q[59];
s q[26];
ccx q[37], q[28], q[34];
h q[6];
ccx q[22], q[18], q[39];
t q[41];
t q[57];
ccx q[53], q[45], q[28];
cx q[52], q[20];
t q[16];
ccx q[20], q[31], q[38];
cx q[32], q[55];
t q[53];
t q[51];
t q[45];
h q[41];
t q[17];
h q[49];
t q[5];
cx q[1], q[43];
cx q[57], q[32];
cx q[54], q[29];
cx q[14], q[8];
h q[12];
ccx q[31], q[0], q[30];
s q[21];
ccx q[9], q[46], q[35];
cx q[3], q[28];
h q[6];
cx q[0], q[50];
h q[35];
h q[35];
h q[19];
t q[28];
ccx q[0], q[13], q[27];
ccx q[42], q[7], q[25];
cx q[26], q[52];
cx q[47], q[40];
t q[21];
t q[4];
h q[51];
t q[56];
ccx q[29], q[53], q[17];
cx q[8], q[24];
h q[17];
t q[39];
cx q[17], q[15];
cx q[13], q[9];
h q[12];
cx q[0], q[22];
t q[9];
s q[56];
ccx q[47], q[5], q[21];
cx q[26], q[21];
h q[24];
h q[49];
h q[28];
cx q[37], q[47];
h q[11];
h q[54];
s q[49];
s q[11];
cx q[42], q[53];
cx q[24], q[49];
h q[21];
t q[52];
s q[0];
ccx q[36], q[25], q[28];
h q[45];
ccx q[15], q[49], q[17];
cx q[21], q[53];
ccx q[34], q[29], q[53];
cx q[15], q[48];
t q[3];
t q[14];
ccx q[42], q[45], q[26];
ccx q[9], q[11], q[33];
cx q[53], q[36];
s q[27];
ccx q[16], q[27], q[54];
cx q[35], q[30];
h q[35];
t q[11];
cx q[33], q[51];
h q[35];
ccx q[30], q[27], q[40];
s q[52];
s q[8];
t q[34];
ccx q[53], q[51], q[23];
ccx q[54], q[9], q[37];
cx q[27], q[56];
ccx q[9], q[52], q[28];
cx q[59], q[12];
cx q[37], q[9];
t q[9];
t q[4];
s q[53];
t q[39];
t q[50];
cx q[35], q[33];
s q[39];
ccx q[59], q[55], q[34];
t q[29];
cx q[21], q[27];
s q[4];
ccx q[49], q[54], q[47];
ccx q[43], q[53], q[7];
t q[9];
ccx q[29], q[17], q[51];
ccx q[23], q[27], q[13];
cx q[29], q[49];
ccx q[57], q[53], q[0];
s q[54];
t q[15];
ccx q[39], q[55], q[30];
t q[7];
h q[45];
h q[58];
t q[37];
cx q[32], q[10];
h q[19];
ccx q[2], q[51], q[3];
h q[28];
t q[51];
t q[16];
t q[18];
t q[42];
h q[47];
ccx q[50], q[47], q[29];
s q[58];
cx q[22], q[45];
t q[51];
cx q[13], q[8];
s q[15];
ccx q[14], q[6], q[50];
s q[45];
cx q[38], q[36];
s q[51];
t q[19];
cx q[17], q[57];
cx q[57], q[38];
t q[9];
cx q[59], q[9];
cx q[8], q[23];
ccx q[38], q[32], q[27];
h q[15];
t q[28];
cx q[11], q[13];
s q[58];
s q[20];
ccx q[3], q[11], q[27];
t q[31];
t q[29];
t q[30];
ccx q[42], q[8], q[1];
s q[14];
ccx q[34], q[51], q[42];
ccx q[40], q[4], q[37];
s q[39];
s q[14];
ccx q[8], q[52], q[38];
s q[45];
h q[17];
h q[48];
t q[29];
ccx q[49], q[18], q[13];
h q[3];
ccx q[8], q[2], q[7];
cx q[37], q[10];
s q[11];
s q[31];
ccx q[50], q[52], q[55];
cx q[51], q[32];
h q[6];
ccx q[55], q[9], q[42];
t q[37];
ccx q[48], q[38], q[15];
ccx q[1], q[9], q[49];
h q[52];
h q[56];
tdg q[0];
cx q[1], q[0];
sdg q[0];
sdg q[2];
cx q[3], q[2];
sdg q[2];
x q[4];
sdg q[6];
tdg q[5];
x q[7];
cx q[8], q[9];
sdg q[8];
tdg q[10];
y q[12];
z q[13];
x q[14];
z q[15];
cx q[17], q[16];
tdg q[16];
sdg q[18];
cx q[18], q[19];
y q[20];
x q[22];
cx q[24], q[25];
sdg q[24];
tdg q[26];
y q[28];
cx q[51], q[52];
t q[43];
ccx q[44], q[48], q[52];
s q[35];
cx q[54], q[48];
t q[42];
cx q[39], q[37];
t q[35];
t q[46];
h q[43];
t q[47];
s q[45];
h q[47];
h q[38];
h q[40];
s q[59];
cx q[34], q[40];
ccx q[31], q[51], q[47];
ccx q[55], q[49], q[43];
h q[54];
h q[48];
h q[33];
s q[40];
ccx q[33], q[49], q[40];
h q[47];
h q[30];
s q[57];
t q[47];
ccx q[31], q[57], q[51];
ccx q[55], q[32], q[49];
cx q[60], q[38];
cx q[61], q[26];
cx q[62], q[45];
cx q[63], q[41];
cx q[64], q[1];
cx q[65], q[12];
cx q[66], q[31];
cx q[67], q[24];
cx q[68], q[28];
cx q[69], q[11];