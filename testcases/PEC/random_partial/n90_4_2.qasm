OPENQASM 2.0;
include "qelib1.inc";
qreg q[106];
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
h q[70];
h q[71];
h q[72];
h q[73];
h q[74];
h q[75];
h q[76];
h q[77];
h q[78];
h q[79];
h q[80];
h q[81];
h q[82];
h q[83];
h q[84];
h q[85];
h q[86];
h q[87];
h q[88];
h q[89];
ccx q[80], q[37], q[1];
ccx q[5], q[76], q[72];
s q[79];
t q[43];
cx q[33], q[28];
ccx q[5], q[30], q[85];
ccx q[81], q[86], q[21];
t q[29];
s q[57];
s q[76];
s q[0];
s q[23];
s q[25];
ccx q[34], q[50], q[89];
cx q[48], q[65];
cx q[29], q[19];
h q[13];
ccx q[63], q[61], q[80];
cx q[13], q[42];
cx q[23], q[88];
h q[3];
ccx q[48], q[49], q[52];
t q[65];
h q[20];
s q[62];
t q[71];
h q[46];
t q[9];
s q[65];
t q[1];
h q[63];
ccx q[2], q[16], q[34];
ccx q[40], q[47], q[46];
t q[80];
ccx q[80], q[25], q[33];
h q[3];
h q[10];
s q[12];
s q[84];
h q[31];
s q[57];
h q[21];
cx q[73], q[56];
t q[45];
cx q[4], q[41];
cx q[28], q[54];
t q[35];
s q[35];
h q[51];
s q[7];
t q[45];
s q[83];
cx q[70], q[68];
s q[40];
t q[82];
cx q[15], q[78];
s q[76];
s q[40];
s q[4];
t q[28];
h q[55];
cx q[82], q[46];
s q[89];
cx q[38], q[6];
ccx q[23], q[8], q[4];
s q[32];
ccx q[87], q[8], q[6];
h q[24];
h q[0];
h q[32];
t q[2];
ccx q[32], q[19], q[29];
ccx q[39], q[25], q[1];
cx q[80], q[67];
t q[6];
t q[29];
t q[86];
t q[61];
h q[12];
ccx q[70], q[89], q[74];
s q[11];
t q[63];
cx q[33], q[61];
h q[48];
h q[51];
h q[7];
s q[42];
cx q[38], q[63];
h q[22];
h q[1];
cx q[58], q[71];
t q[37];
cx q[7], q[11];
cx q[79], q[72];
ccx q[76], q[69], q[23];
cx q[10], q[50];
ccx q[84], q[10], q[38];
ccx q[14], q[62], q[32];
s q[12];
s q[2];
s q[58];
t q[36];
h q[34];
cx q[83], q[56];
h q[13];
cx q[14], q[86];
s q[16];
cx q[67], q[56];
t q[52];
cx q[20], q[73];
h q[46];
s q[21];
h q[57];
s q[62];
t q[20];
ccx q[44], q[51], q[77];
s q[42];
ccx q[13], q[23], q[29];
cx q[0], q[56];
h q[10];
s q[6];
s q[61];
ccx q[77], q[82], q[10];
t q[34];
ccx q[33], q[15], q[87];
h q[37];
cx q[60], q[0];
ccx q[55], q[81], q[30];
t q[13];
h q[25];
t q[44];
t q[26];
t q[78];
cx q[50], q[82];
t q[27];
s q[1];
t q[24];
h q[14];
t q[68];
ccx q[42], q[60], q[4];
h q[14];
h q[53];
cx q[72], q[74];
s q[83];
cx q[4], q[64];
cx q[1], q[65];
ccx q[10], q[72], q[51];
ccx q[46], q[6], q[88];
ccx q[41], q[26], q[17];
cx q[13], q[65];
t q[52];
s q[3];
cx q[84], q[55];
ccx q[83], q[6], q[38];
s q[16];
s q[39];
cx q[53], q[45];
cx q[58], q[50];
s q[20];
ccx q[65], q[63], q[66];
h q[30];
ccx q[61], q[23], q[70];
t q[5];
ccx q[29], q[21], q[24];
t q[83];
h q[46];
t q[83];
h q[29];
ccx q[83], q[49], q[18];
cx q[21], q[53];
s q[49];
ccx q[80], q[55], q[34];
cx q[82], q[89];
t q[25];
ccx q[38], q[59], q[80];
cx q[86], q[76];
cx q[59], q[73];
cx q[32], q[59];
h q[19];
t q[4];
ccx q[62], q[88], q[0];
cx q[70], q[18];
h q[65];
ccx q[11], q[51], q[24];
h q[49];
t q[27];
t q[56];
cx q[3], q[44];
h q[14];
s q[75];
cx q[28], q[38];
t q[65];
t q[83];
s q[44];
h q[27];
t q[85];
cx q[34], q[50];
h q[87];
cx q[27], q[46];
t q[14];
t q[37];
h q[64];
t q[16];
t q[40];
cx q[89], q[17];
h q[14];
s q[56];
ccx q[32], q[58], q[72];
ccx q[13], q[89], q[54];
ccx q[24], q[53], q[19];
ccx q[36], q[29], q[37];
cx q[17], q[10];
cx q[76], q[45];
ccx q[75], q[85], q[58];
t q[53];
cx q[62], q[10];
t q[29];
t q[79];
t q[67];
s q[76];
ccx q[85], q[19], q[74];
s q[58];
cx q[70], q[74];
t q[4];
cx q[23], q[75];
ccx q[21], q[16], q[3];
ccx q[17], q[71], q[14];
cx q[89], q[30];
h q[45];
ccx q[49], q[75], q[34];
cx q[72], q[20];
s q[19];
h q[89];
t q[16];
cx q[65], q[63];
cx q[84], q[41];
ccx q[80], q[21], q[29];
ccx q[56], q[66], q[25];
s q[3];
t q[86];
s q[39];
s q[20];
s q[48];
s q[45];
s q[0];
t q[51];
s q[66];
cx q[42], q[81];
t q[28];
cx q[15], q[37];
h q[38];
s q[50];
h q[56];
s q[20];
ccx q[88], q[73], q[8];
t q[30];
h q[75];
ccx q[71], q[27], q[78];
ccx q[29], q[55], q[3];
t q[67];
s q[13];
ccx q[82], q[54], q[22];
t q[40];
t q[72];
cx q[41], q[89];
cx q[88], q[39];
ccx q[72], q[39], q[12];
cx q[7], q[42];
ccx q[19], q[79], q[27];
ccx q[54], q[61], q[37];
x q[0];
x q[1];
cx q[3], q[2];
sdg q[2];
cx q[3], q[2];
sdg q[2];
sdg q[5];
sdg q[4];
cx q[5], q[4];
sdg q[4];
tdg q[6];
cx q[6], q[7];
sdg q[8];
cx q[11], q[12];
cx q[13], q[14];
sdg q[13];
tdg q[15];
y q[18];
cx q[20], q[19];
tdg q[19];
sdg q[19];
cx q[20], q[19];
tdg q[19];
x q[21];
tdg q[21];
cx q[22], q[21];
x q[21];
cx q[24], q[23];
x q[23];
sdg q[25];
tdg q[28];
cx q[27], q[28];
tdg q[28];
cx q[27], q[28];
tdg q[27];
sdg q[29];
tdg q[30];
tdg q[29];
cx q[29], q[30];
tdg q[32];
cx q[32], q[31];
sdg q[31];
sdg q[33];
cx q[35], q[36];
cx q[37], q[38];
tdg q[39];
cx q[40], q[39];
x q[39];
z q[41];
tdg q[43];
sdg q[42];
cx q[42], q[43];
h q[53];
cx q[50], q[49];
ccx q[75], q[64], q[82];
s q[67];
t q[79];
cx q[52], q[56];
cx q[62], q[59];
cx q[54], q[65];
cx q[83], q[73];
s q[78];
cx q[54], q[79];
h q[62];
h q[58];
ccx q[69], q[46], q[72];
s q[80];
h q[71];
ccx q[64], q[80], q[67];
cx q[66], q[78];
cx q[72], q[86];
t q[85];
ccx q[86], q[68], q[88];
h q[73];
h q[51];
h q[62];
t q[80];
h q[55];
t q[87];
ccx q[71], q[46], q[87];
ccx q[64], q[75], q[68];
t q[60];
ccx q[68], q[45], q[74];
s q[55];
cx q[74], q[75];
s q[60];
s q[53];
s q[81];
h q[72];
t q[84];
s q[54];
h q[47];
h q[67];
t q[85];
s q[78];
h q[73];
ccx q[73], q[60], q[45];
cx q[90], q[34];
cx q[91], q[3];
cx q[92], q[49];
cx q[93], q[84];
cx q[94], q[40];
cx q[95], q[73];
cx q[96], q[45];
cx q[97], q[1];
cx q[98], q[63];
cx q[99], q[5];
cx q[100], q[22];
cx q[101], q[56];
cx q[102], q[70];
cx q[103], q[16];
cx q[104], q[80];
cx q[105], q[20];