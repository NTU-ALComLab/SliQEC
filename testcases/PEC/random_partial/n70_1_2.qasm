OPENQASM 2.0;
include "qelib1.inc";
qreg q[83];
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
t q[37];
cx q[9], q[41];
s q[55];
s q[7];
s q[17];
ccx q[6], q[53], q[51];
h q[17];
s q[7];
ccx q[48], q[67], q[16];
h q[52];
cx q[5], q[7];
ccx q[59], q[22], q[15];
h q[50];
s q[16];
h q[38];
t q[68];
t q[18];
cx q[13], q[51];
t q[26];
cx q[56], q[0];
ccx q[10], q[52], q[23];
h q[23];
s q[6];
t q[55];
t q[68];
s q[61];
cx q[50], q[21];
ccx q[50], q[26], q[41];
cx q[24], q[51];
h q[49];
h q[14];
cx q[8], q[24];
ccx q[0], q[68], q[52];
cx q[38], q[20];
cx q[44], q[34];
t q[24];
h q[24];
s q[35];
t q[54];
ccx q[54], q[6], q[53];
t q[34];
s q[17];
cx q[31], q[41];
cx q[18], q[50];
h q[68];
s q[58];
h q[40];
ccx q[24], q[55], q[37];
s q[18];
h q[46];
cx q[15], q[54];
cx q[46], q[28];
h q[11];
ccx q[60], q[43], q[14];
t q[54];
t q[29];
t q[52];
cx q[59], q[27];
s q[41];
t q[49];
h q[29];
ccx q[47], q[44], q[35];
h q[44];
t q[25];
t q[11];
s q[26];
h q[47];
t q[45];
h q[25];
ccx q[48], q[26], q[5];
h q[21];
s q[27];
cx q[48], q[46];
h q[67];
t q[35];
cx q[38], q[15];
ccx q[66], q[47], q[48];
cx q[0], q[49];
t q[49];
s q[43];
ccx q[61], q[30], q[41];
h q[54];
s q[16];
h q[56];
cx q[65], q[24];
t q[50];
s q[43];
t q[7];
ccx q[7], q[23], q[14];
s q[40];
ccx q[57], q[36], q[34];
s q[32];
s q[60];
s q[39];
ccx q[66], q[39], q[23];
cx q[48], q[38];
ccx q[17], q[27], q[36];
h q[5];
s q[32];
cx q[51], q[61];
t q[37];
cx q[25], q[16];
ccx q[26], q[1], q[40];
cx q[60], q[14];
h q[17];
ccx q[9], q[24], q[4];
t q[53];
h q[29];
cx q[5], q[36];
h q[54];
t q[32];
h q[40];
h q[60];
t q[13];
s q[31];
s q[50];
h q[7];
t q[32];
ccx q[51], q[32], q[34];
s q[60];
cx q[58], q[14];
cx q[2], q[50];
ccx q[47], q[18], q[37];
cx q[13], q[25];
h q[65];
t q[51];
ccx q[67], q[35], q[2];
t q[11];
cx q[38], q[41];
s q[34];
t q[54];
cx q[21], q[13];
s q[29];
s q[37];
s q[8];
h q[23];
t q[43];
h q[4];
h q[16];
cx q[41], q[58];
ccx q[7], q[3], q[27];
t q[41];
ccx q[0], q[39], q[66];
t q[29];
t q[65];
t q[2];
s q[51];
cx q[65], q[53];
ccx q[39], q[36], q[3];
t q[5];
h q[59];
h q[60];
t q[35];
ccx q[43], q[37], q[21];
h q[12];
s q[53];
s q[52];
t q[47];
cx q[20], q[2];
t q[10];
h q[38];
cx q[0], q[31];
t q[31];
s q[15];
t q[23];
cx q[41], q[56];
h q[7];
s q[53];
h q[63];
h q[14];
cx q[45], q[13];
s q[13];
t q[4];
s q[28];
t q[23];
ccx q[5], q[47], q[53];
ccx q[16], q[39], q[49];
s q[56];
h q[51];
h q[48];
t q[7];
h q[22];
ccx q[66], q[17], q[23];
cx q[9], q[64];
ccx q[36], q[26], q[39];
cx q[60], q[61];
cx q[68], q[26];
s q[41];
h q[26];
t q[18];
s q[51];
t q[49];
ccx q[13], q[55], q[1];
h q[63];
ccx q[11], q[4], q[38];
t q[14];
t q[27];
t q[54];
ccx q[61], q[1], q[17];
cx q[8], q[42];
t q[37];
cx q[24], q[52];
s q[56];
h q[59];
h q[30];
h q[53];
s q[36];
h q[0];
s q[44];
s q[26];
z q[1];
tdg q[3];
sdg q[2];
cx q[2], q[3];
cx q[4], q[5];
sdg q[5];
sdg q[4];
cx q[4], q[5];
tdg q[4];
cx q[6], q[7];
sdg q[8];
cx q[10], q[11];
tdg q[12];
cx q[13], q[12];
sdg q[12];
cx q[13], q[12];
tdg q[12];
x q[14];
cx q[15], q[16];
tdg q[16];
sdg q[15];
cx q[15], q[16];
x q[17];
cx q[18], q[19];
cx q[20], q[21];
sdg q[21];
cx q[20], q[21];
sdg q[21];
sdg q[20];
y q[22];
x q[23];
cx q[25], q[26];
cx q[28], q[27];
x q[27];
sdg q[27];
cx q[28], q[27];
x q[27];
cx q[29], q[30];
tdg q[30];
cx q[29], q[30];
tdg q[29];
tdg q[31];
cx q[31], q[32];
cx q[33], q[34];
tdg q[33];
s q[50];
h q[43];
t q[48];
h q[40];
s q[55];
ccx q[55], q[41], q[43];
cx q[43], q[49];
t q[67];
ccx q[61], q[53], q[38];
s q[57];
s q[40];
h q[46];
ccx q[61], q[63], q[42];
s q[52];
cx q[55], q[61];
cx q[52], q[44];
h q[36];
cx q[47], q[52];
t q[36];
h q[53];
h q[57];
ccx q[64], q[37], q[59];
ccx q[40], q[49], q[66];
t q[61];
ccx q[49], q[42], q[65];
cx q[59], q[41];
ccx q[52], q[68], q[65];
h q[53];
t q[52];
h q[46];
h q[62];
h q[57];
h q[50];
ccx q[54], q[56], q[39];
cx q[52], q[69];
cx q[70], q[33];
cx q[71], q[24];
cx q[72], q[17];
cx q[73], q[9];
cx q[74], q[54];
cx q[75], q[3];
cx q[76], q[39];
cx q[77], q[46];
cx q[78], q[31];
cx q[79], q[48];
cx q[80], q[0];
cx q[81], q[50];
cx q[82], q[38];