OPENQASM 2.0;
include "qelib1.inc";
qreg q[60];
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
t q[52];
ccx q[3], q[2], q[27];
cx q[6], q[28];
t q[26];
ccx q[18], q[26], q[36];
cx q[41], q[49];
t q[2];
h q[58];
s q[54];
cx q[22], q[36];
s q[53];
t q[16];
h q[17];
cx q[26], q[40];
s q[58];
cx q[9], q[8];
s q[53];
h q[34];
cx q[21], q[37];
s q[24];
h q[26];
t q[20];
ccx q[29], q[41], q[18];
h q[53];
ccx q[33], q[52], q[16];
cx q[31], q[0];
ccx q[4], q[59], q[33];
s q[25];
t q[11];
h q[36];
h q[35];
s q[37];
cx q[4], q[38];
ccx q[9], q[25], q[47];
cx q[47], q[7];
ccx q[50], q[54], q[11];
ccx q[1], q[43], q[41];
h q[12];
s q[47];
ccx q[58], q[23], q[8];
cx q[35], q[18];
t q[17];
h q[17];
ccx q[6], q[58], q[1];
h q[57];
h q[51];
s q[24];
ccx q[57], q[30], q[35];
cx q[25], q[52];
ccx q[58], q[11], q[25];
t q[13];
ccx q[53], q[0], q[39];
s q[28];
ccx q[6], q[15], q[51];
ccx q[38], q[41], q[28];
ccx q[5], q[33], q[14];
h q[17];
t q[49];
ccx q[1], q[37], q[54];
s q[2];
ccx q[58], q[24], q[19];
ccx q[22], q[4], q[0];
s q[10];
cx q[45], q[9];
ccx q[29], q[23], q[17];
cx q[15], q[34];
s q[44];
s q[59];
cx q[25], q[14];
ccx q[39], q[42], q[58];
h q[1];
ccx q[15], q[37], q[56];
ccx q[24], q[7], q[51];
s q[9];
h q[14];
s q[6];
t q[52];
h q[57];
cx q[45], q[38];
h q[48];
h q[28];
t q[28];
ccx q[35], q[55], q[54];
s q[1];
h q[23];
h q[24];
s q[31];
ccx q[54], q[22], q[17];
s q[23];
t q[9];
t q[29];
cx q[22], q[10];
cx q[47], q[17];
cx q[43], q[15];
cx q[54], q[49];
ccx q[34], q[51], q[50];
cx q[4], q[6];
t q[52];
h q[15];
h q[57];
ccx q[14], q[21], q[33];
cx q[36], q[28];
cx q[36], q[23];
ccx q[30], q[38], q[16];
cx q[52], q[6];
h q[44];
t q[31];
s q[38];
t q[42];
s q[9];
ccx q[39], q[24], q[55];
s q[57];
s q[56];
h q[30];
s q[55];
t q[45];
s q[20];
s q[38];
cx q[0], q[56];
ccx q[48], q[58], q[6];
cx q[7], q[49];
t q[3];
ccx q[59], q[1], q[54];
h q[10];
s q[51];
s q[11];
h q[14];
h q[36];
ccx q[5], q[27], q[41];
s q[11];
cx q[10], q[53];
h q[18];
cx q[27], q[26];
s q[56];
h q[33];
ccx q[43], q[0], q[50];
h q[29];
t q[7];
h q[35];
h q[30];
s q[10];
cx q[38], q[5];
cx q[32], q[44];
ccx q[39], q[10], q[20];
cx q[50], q[1];
s q[40];
s q[7];
ccx q[48], q[45], q[39];
t q[11];
s q[54];
t q[2];
s q[16];
ccx q[11], q[49], q[16];
cx q[12], q[14];
s q[4];
s q[43];
ccx q[29], q[52], q[21];
ccx q[6], q[51], q[58];
s q[38];
ccx q[16], q[53], q[55];
cx q[0], q[32];
h q[38];
t q[33];
t q[52];
s q[39];
s q[38];
t q[17];
ccx q[9], q[12], q[53];
h q[7];
ccx q[1], q[44], q[24];
s q[19];
h q[6];
ccx q[49], q[23], q[59];
ccx q[53], q[58], q[56];
h q[16];
t q[41];
t q[1];
ccx q[11], q[4], q[34];
cx q[10], q[7];
t q[17];