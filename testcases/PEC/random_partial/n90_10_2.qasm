OPENQASM 2.0;
include "qelib1.inc";
qreg q[99];
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
cx q[65], q[24];
h q[9];
cx q[80], q[76];
ccx q[57], q[82], q[61];
t q[66];
h q[27];
cx q[71], q[83];
t q[80];
s q[81];
t q[78];
h q[56];
s q[25];
ccx q[22], q[47], q[78];
s q[9];
h q[28];
t q[68];
cx q[42], q[32];
t q[71];
h q[87];
t q[48];
t q[54];
t q[70];
h q[15];
ccx q[72], q[84], q[36];
h q[64];
h q[31];
h q[69];
ccx q[48], q[53], q[1];
t q[48];
cx q[78], q[21];
t q[34];
cx q[35], q[76];
s q[15];
s q[53];
s q[60];
t q[47];
t q[11];
t q[45];
ccx q[2], q[79], q[49];
ccx q[51], q[54], q[29];
cx q[43], q[29];
cx q[79], q[43];
ccx q[56], q[9], q[59];
h q[84];
ccx q[48], q[27], q[69];
t q[10];
t q[85];
s q[2];
ccx q[25], q[63], q[57];
t q[31];
s q[69];
h q[59];
cx q[34], q[23];
cx q[75], q[85];
s q[48];
t q[20];
t q[20];
t q[13];
ccx q[87], q[6], q[14];
t q[69];
h q[57];
t q[18];
t q[67];
t q[89];
cx q[25], q[66];
s q[85];
h q[42];
h q[35];
s q[58];
ccx q[30], q[38], q[44];
s q[49];
ccx q[61], q[13], q[47];
cx q[30], q[28];
s q[84];
cx q[61], q[22];
h q[25];
h q[38];
ccx q[9], q[67], q[29];
s q[30];
t q[6];
cx q[87], q[83];
ccx q[45], q[21], q[33];
cx q[25], q[26];
cx q[38], q[66];
ccx q[34], q[68], q[17];
t q[67];
t q[32];
cx q[89], q[28];
cx q[38], q[4];
t q[45];
cx q[36], q[14];
s q[72];
t q[47];
cx q[9], q[31];
cx q[77], q[55];
ccx q[33], q[52], q[26];
ccx q[81], q[43], q[50];
cx q[14], q[1];
s q[3];
ccx q[87], q[39], q[13];
h q[32];
h q[16];
cx q[35], q[12];
t q[24];
h q[30];
s q[8];
t q[41];
t q[80];
s q[58];
s q[87];
cx q[41], q[74];
ccx q[8], q[20], q[63];
cx q[66], q[55];
t q[30];
ccx q[71], q[5], q[84];
ccx q[50], q[64], q[10];
ccx q[22], q[87], q[88];
h q[51];
s q[61];
h q[41];
ccx q[55], q[65], q[45];
ccx q[17], q[49], q[66];
t q[52];
h q[11];
cx q[16], q[31];
h q[64];
t q[83];
s q[15];
t q[13];
h q[36];
s q[58];
cx q[88], q[85];
cx q[65], q[44];
cx q[35], q[71];
cx q[1], q[83];
ccx q[18], q[72], q[82];
s q[9];
cx q[35], q[74];
ccx q[44], q[81], q[0];
ccx q[32], q[54], q[1];
ccx q[58], q[46], q[16];
cx q[8], q[58];
cx q[57], q[56];
h q[64];
t q[15];
s q[25];
h q[68];
s q[84];
h q[33];
s q[13];
cx q[78], q[52];
s q[46];
cx q[43], q[32];
h q[12];
t q[47];
ccx q[66], q[23], q[13];
cx q[53], q[44];
cx q[67], q[5];
ccx q[24], q[63], q[77];
ccx q[49], q[24], q[87];
t q[56];
h q[12];
t q[23];
h q[51];
ccx q[5], q[84], q[4];
ccx q[25], q[1], q[64];
t q[9];
h q[14];
t q[38];
s q[17];
s q[33];
cx q[52], q[50];
t q[83];
h q[1];
s q[58];
s q[76];
ccx q[81], q[21], q[31];
s q[82];
s q[30];
t q[7];
ccx q[77], q[17], q[89];
t q[33];
cx q[55], q[72];
cx q[21], q[50];
cx q[12], q[81];
ccx q[67], q[75], q[85];
ccx q[83], q[76], q[49];
t q[51];
h q[4];
ccx q[50], q[15], q[75];
s q[79];
cx q[2], q[6];
t q[11];
cx q[35], q[59];
h q[12];
h q[6];
s q[32];
h q[57];
ccx q[3], q[16], q[25];
h q[87];
ccx q[11], q[29], q[87];
h q[3];
h q[39];
cx q[31], q[5];
ccx q[66], q[56], q[64];
t q[53];
s q[25];
t q[38];
cx q[16], q[42];
h q[52];
s q[39];
cx q[29], q[81];
ccx q[8], q[12], q[2];
h q[7];
h q[85];
h q[29];
s q[37];
cx q[47], q[24];
t q[73];
ccx q[56], q[78], q[60];
t q[84];
cx q[74], q[58];
ccx q[52], q[17], q[81];
s q[23];
s q[36];
s q[56];
cx q[26], q[85];
ccx q[10], q[65], q[8];
s q[12];
cx q[74], q[25];
h q[72];
t q[1];
t q[87];
cx q[48], q[39];
s q[43];
s q[44];
h q[38];
cx q[78], q[72];
s q[87];
t q[59];
ccx q[70], q[82], q[31];
cx q[26], q[1];
h q[65];
t q[12];
h q[82];
cx q[21], q[17];
cx q[44], q[37];
h q[8];
cx q[55], q[78];
s q[35];
ccx q[34], q[36], q[2];
ccx q[21], q[66], q[58];
t q[22];
s q[47];
s q[8];
s q[27];
ccx q[64], q[61], q[0];
cx q[35], q[69];
ccx q[42], q[22], q[71];
h q[76];
s q[70];
cx q[78], q[1];
h q[86];
t q[45];
ccx q[22], q[30], q[80];
s q[27];
s q[77];
t q[7];
t q[18];
h q[45];
sdg q[1];
sdg q[0];
sdg q[2];
cx q[5], q[4];
sdg q[4];
cx q[5], q[4];
tdg q[4];
x q[6];
sdg q[8];
cx q[8], q[7];
tdg q[7];
y q[9];
cx q[11], q[10];
tdg q[10];
z q[13];
tdg q[15];
sdg q[14];
cx q[14], q[15];
sdg q[16];
y q[18];
sdg q[19];
tdg q[19];
cx q[21], q[22];
tdg q[21];
sdg q[22];
cx q[21], q[22];
tdg q[21];
cx q[24], q[23];
tdg q[23];
cx q[24], q[23];
sdg q[23];
cx q[25], q[26];
sdg q[26];
sdg q[25];
sdg q[27];
cx q[27], q[28];
sdg q[29];
sdg q[30];
cx q[30], q[29];
tdg q[29];
x q[31];
z q[32];
sdg q[33];
cx q[34], q[33];
tdg q[34];
tdg q[33];
cx q[35], q[36];
sdg q[36];
cx q[35], q[36];
sdg q[35];
cx q[38], q[37];
tdg q[37];
tdg q[39];
cx q[42], q[41];
tdg q[41];
sdg q[43];
cx q[44], q[43];
x q[43];
h q[46];
s q[65];
cx q[74], q[62];
s q[50];
s q[89];
t q[86];
h q[51];
h q[60];
ccx q[66], q[89], q[75];
h q[78];
ccx q[58], q[61], q[68];
ccx q[63], q[65], q[84];
t q[89];
h q[65];
s q[68];
ccx q[85], q[86], q[57];
t q[87];
s q[58];
s q[85];
t q[70];
s q[71];
t q[47];
h q[64];
s q[63];
s q[78];
h q[59];
ccx q[62], q[65], q[54];
ccx q[50], q[47], q[59];
cx q[46], q[57];
h q[78];
ccx q[69], q[58], q[81];
s q[61];
s q[67];
ccx q[89], q[79], q[86];
s q[87];
ccx q[87], q[50], q[75];
cx q[52], q[88];
h q[57];
s q[85];
s q[74];
cx q[88], q[58];
s q[64];
ccx q[56], q[59], q[74];
ccx q[58], q[52], q[70];
ccx q[81], q[74], q[80];
cx q[90], q[70];
cx q[91], q[23];
cx q[92], q[87];
cx q[93], q[52];
cx q[94], q[44];
cx q[95], q[1];
cx q[96], q[81];
cx q[97], q[18];
cx q[98], q[41];