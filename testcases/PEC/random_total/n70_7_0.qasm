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
h q[60];
h q[61];
h q[62];
h q[63];
h q[64];
h q[65];
h q[66];
h q[67];
h q[68];
h q[69];
cx q[7], q[11];
h q[17];
h q[36];
ccx q[23], q[32], q[3];
h q[39];
ccx q[0], q[2], q[43];
s q[2];
ccx q[16], q[43], q[19];
s q[2];
h q[39];
cx q[49], q[54];
t q[22];
h q[1];
h q[28];
cx q[60], q[58];
s q[54];
h q[12];
ccx q[68], q[44], q[59];
h q[35];
cx q[38], q[9];
ccx q[18], q[39], q[43];
t q[37];
cx q[6], q[18];
cx q[18], q[32];
cx q[23], q[46];
h q[8];
h q[45];
ccx q[69], q[23], q[0];
h q[46];
ccx q[54], q[63], q[19];
t q[66];
cx q[33], q[36];
h q[31];
h q[17];
t q[29];
s q[47];
t q[14];
t q[14];
ccx q[58], q[42], q[65];
t q[48];
ccx q[17], q[20], q[31];
ccx q[59], q[42], q[45];
t q[62];
ccx q[20], q[50], q[28];
s q[29];
t q[50];
ccx q[41], q[11], q[16];
h q[40];
s q[23];
t q[53];
ccx q[10], q[41], q[68];
h q[68];
h q[19];
h q[10];
t q[3];
cx q[44], q[1];
ccx q[28], q[20], q[23];
ccx q[20], q[35], q[55];
cx q[6], q[36];
h q[16];
t q[45];
ccx q[48], q[60], q[58];
ccx q[17], q[51], q[13];
t q[62];
ccx q[62], q[63], q[11];
t q[30];
cx q[22], q[48];
cx q[12], q[45];
cx q[3], q[15];
ccx q[48], q[42], q[31];
cx q[60], q[12];
h q[55];
h q[50];
s q[65];
cx q[17], q[3];
s q[61];
s q[65];
t q[32];
cx q[57], q[1];
t q[37];
t q[63];
h q[14];
h q[15];
cx q[41], q[62];
cx q[39], q[60];
ccx q[2], q[50], q[54];
cx q[62], q[0];
ccx q[49], q[68], q[51];
s q[62];
h q[35];
cx q[21], q[43];
h q[38];
cx q[4], q[60];
ccx q[30], q[44], q[2];
cx q[59], q[63];
cx q[30], q[26];
s q[29];
cx q[66], q[36];
s q[12];
cx q[2], q[5];
s q[33];
s q[67];
t q[39];
cx q[55], q[21];
h q[37];
cx q[16], q[4];
cx q[33], q[36];
t q[49];
ccx q[6], q[58], q[47];
t q[22];
cx q[8], q[42];
t q[51];
h q[27];
h q[52];
cx q[28], q[60];
ccx q[20], q[41], q[10];
cx q[44], q[17];
cx q[32], q[6];
ccx q[14], q[68], q[63];
h q[43];
s q[37];
cx q[53], q[33];
t q[50];
cx q[69], q[47];
h q[69];
s q[28];
s q[48];
cx q[67], q[50];
ccx q[32], q[4], q[3];
h q[51];
ccx q[60], q[40], q[46];
ccx q[16], q[26], q[3];
cx q[66], q[68];
ccx q[0], q[44], q[64];
ccx q[28], q[30], q[39];
s q[64];
s q[66];
s q[15];
cx q[62], q[61];
h q[21];
ccx q[37], q[33], q[38];
ccx q[58], q[23], q[54];
t q[9];
t q[39];
s q[32];
cx q[15], q[19];
s q[64];
s q[41];
ccx q[60], q[58], q[37];
cx q[68], q[50];
s q[4];
s q[67];
ccx q[41], q[35], q[28];
h q[63];
cx q[25], q[42];
cx q[50], q[67];
cx q[67], q[48];
s q[17];
ccx q[45], q[8], q[7];
ccx q[68], q[25], q[64];
s q[57];
h q[15];
h q[22];
h q[10];
cx q[60], q[47];
t q[9];
cx q[32], q[26];
cx q[24], q[36];
cx q[7], q[8];
t q[1];
s q[49];
t q[22];
cx q[27], q[55];
cx q[43], q[54];
ccx q[4], q[14], q[16];
cx q[19], q[40];
h q[51];
h q[35];
h q[52];
t q[6];
s q[30];
s q[13];
s q[24];
s q[47];
ccx q[69], q[37], q[59];
ccx q[43], q[62], q[49];
s q[30];
ccx q[33], q[35], q[19];
h q[22];
ccx q[39], q[18], q[0];
ccx q[11], q[29], q[69];
t q[27];
ccx q[25], q[29], q[46];
t q[45];
h q[22];
s q[1];
t q[32];
s q[8];
t q[35];
h q[68];
t q[14];
cx q[34], q[6];
h q[59];
s q[47];
cx q[67], q[37];
h q[49];
ccx q[28], q[41], q[38];
t q[55];
t q[53];
cx q[63], q[38];