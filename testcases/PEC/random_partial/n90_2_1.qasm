OPENQASM 2.0;
include "qelib1.inc";
qreg q[90];
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
h q[87];
h q[88];
h q[82];
h q[10];
h q[3];
cx q[3], q[46];
tdg q[46];
cx q[12], q[46];
t q[46];
cx q[3], q[46];
tdg q[46];
cx q[12], q[46];
t q[46];
cx q[12], q[3];
tdg q[3];
cx q[12], q[3];
t q[12];
t q[3];
h q[3];
t q[75];
s q[15];
t q[52];
s q[12];
s q[80];
t q[66];
t q[33];
cx q[36], q[72];
t q[45];
s q[28];
s q[2];
t q[85];
cx q[38], q[27];
h q[12];
s q[42];
h q[11];
cx q[11], q[76];
tdg q[76];
cx q[41], q[76];
t q[76];
cx q[11], q[76];
tdg q[76];
cx q[41], q[76];
t q[76];
cx q[41], q[11];
tdg q[11];
cx q[41], q[11];
t q[41];
t q[11];
h q[11];
h q[22];
cx q[22], q[38];
tdg q[38];
cx q[66], q[38];
t q[38];
cx q[22], q[38];
tdg q[38];
cx q[66], q[38];
t q[38];
cx q[66], q[22];
tdg q[22];
cx q[66], q[22];
t q[66];
t q[22];
h q[22];
s q[72];
h q[70];
h q[32];
cx q[32], q[37];
tdg q[37];
cx q[86], q[37];
t q[37];
cx q[32], q[37];
tdg q[37];
cx q[86], q[37];
t q[37];
cx q[86], q[32];
tdg q[32];
cx q[86], q[32];
t q[86];
t q[32];
h q[32];
h q[18];
cx q[18], q[29];
tdg q[29];
cx q[5], q[29];
t q[29];
cx q[18], q[29];
tdg q[29];
cx q[5], q[29];
t q[29];
cx q[5], q[18];
tdg q[18];
cx q[5], q[18];
t q[5];
t q[18];
h q[18];
t q[14];
cx q[87], q[51];
cx q[83], q[16];
h q[21];
cx q[21], q[83];
tdg q[83];
cx q[26], q[83];
t q[83];
cx q[21], q[83];
tdg q[83];
cx q[26], q[83];
t q[83];
cx q[26], q[21];
tdg q[21];
cx q[26], q[21];
t q[26];
t q[21];
h q[21];
cx q[83], q[6];
cx q[27], q[53];
s q[83];
s q[15];
cx q[83], q[67];
cx q[31], q[20];
h q[69];
cx q[72], q[44];
s q[60];
cx q[46], q[79];
t q[18];
h q[86];
cx q[86], q[74];
tdg q[74];
cx q[22], q[74];
t q[74];
cx q[86], q[74];
tdg q[74];
cx q[22], q[74];
t q[74];
cx q[22], q[86];
tdg q[86];
cx q[22], q[86];
t q[22];
t q[86];
h q[86];
h q[53];
cx q[53], q[65];
tdg q[65];
cx q[15], q[65];
t q[65];
cx q[53], q[65];
tdg q[65];
cx q[15], q[65];
t q[65];
cx q[15], q[53];
tdg q[53];
cx q[15], q[53];
t q[15];
t q[53];
h q[53];
cx q[9], q[14];
h q[39];
cx q[39], q[25];
tdg q[25];
cx q[79], q[25];
t q[25];
cx q[39], q[25];
tdg q[25];
cx q[79], q[25];
t q[25];
cx q[79], q[39];
tdg q[39];
cx q[79], q[39];
t q[79];
t q[39];
h q[39];
s q[16];
t q[28];
t q[36];
h q[73];
h q[12];
cx q[12], q[22];
tdg q[22];
cx q[23], q[22];
t q[22];
cx q[12], q[22];
tdg q[22];
cx q[23], q[22];
t q[22];
cx q[23], q[12];
tdg q[12];
cx q[23], q[12];
t q[23];
t q[12];
h q[12];
cx q[32], q[52];
cx q[10], q[48];
t q[1];
h q[56];
cx q[56], q[46];
tdg q[46];
cx q[12], q[46];
t q[46];
cx q[56], q[46];
tdg q[46];
cx q[12], q[46];
t q[46];
cx q[12], q[56];
tdg q[56];
cx q[12], q[56];
t q[12];
t q[56];
h q[56];
t q[34];
h q[61];
h q[77];
cx q[77], q[73];
tdg q[73];
cx q[49], q[73];
t q[73];
cx q[77], q[73];
tdg q[73];
cx q[49], q[73];
t q[73];
cx q[49], q[77];
tdg q[77];
cx q[49], q[77];
t q[49];
t q[77];
h q[77];
s q[6];
h q[64];
cx q[64], q[83];
tdg q[83];
cx q[18], q[83];
t q[83];
cx q[64], q[83];
tdg q[83];
cx q[18], q[83];
t q[83];
cx q[18], q[64];
tdg q[64];
cx q[18], q[64];
t q[18];
t q[64];
h q[64];
s q[73];
h q[1];
cx q[15], q[55];
h q[70];
s q[0];
s q[34];
t q[82];
s q[24];
s q[54];
t q[84];
cx q[8], q[76];
h q[54];
cx q[54], q[59];
tdg q[59];
cx q[76], q[59];
t q[59];
cx q[54], q[59];
tdg q[59];
cx q[76], q[59];
t q[59];
cx q[76], q[54];
tdg q[54];
cx q[76], q[54];
t q[76];
t q[54];
h q[54];
h q[73];
cx q[73], q[23];
tdg q[23];
cx q[82], q[23];
t q[23];
cx q[73], q[23];
tdg q[23];
cx q[82], q[23];
t q[23];
cx q[82], q[73];
tdg q[73];
cx q[82], q[73];
t q[82];
t q[73];
h q[73];
h q[71];
h q[76];
cx q[76], q[13];
tdg q[13];
cx q[34], q[13];
t q[13];
cx q[76], q[13];
tdg q[13];
cx q[34], q[13];
t q[13];
cx q[34], q[76];
tdg q[76];
cx q[34], q[76];
t q[34];
t q[76];
h q[76];
h q[9];
cx q[9], q[87];
tdg q[87];
cx q[16], q[87];
t q[87];
cx q[9], q[87];
tdg q[87];
cx q[16], q[87];
t q[87];
cx q[16], q[9];
tdg q[9];
cx q[16], q[9];
t q[16];
t q[9];
h q[9];
h q[29];
cx q[29], q[79];
tdg q[79];
cx q[50], q[79];
t q[79];
cx q[29], q[79];
tdg q[79];
cx q[50], q[79];
t q[79];
cx q[50], q[29];
tdg q[29];
cx q[50], q[29];
t q[50];
t q[29];
h q[29];
h q[38];
h q[52];
s q[58];
t q[18];
cx q[3], q[37];
h q[42];
cx q[36], q[69];
s q[21];
h q[11];
cx q[11], q[60];
tdg q[60];
cx q[73], q[60];
t q[60];
cx q[11], q[60];
tdg q[60];
cx q[73], q[60];
t q[60];
cx q[73], q[11];
tdg q[11];
cx q[73], q[11];
t q[73];
t q[11];
h q[11];
s q[31];
h q[87];
s q[22];
s q[77];
t q[81];
t q[11];
t q[84];
h q[72];
cx q[72], q[48];
s q[48];
cx q[84], q[72];
t q[49];
h q[39];
h q[75];
cx q[4], q[62];
h q[22];
h q[35];
cx q[40], q[16];
t q[73];
s q[81];
h q[9];
h q[11];
h q[59];
s q[81];
h q[68];
t q[66];
h q[52];
cx q[52], q[54];
tdg q[54];
cx q[74], q[54];
t q[54];
cx q[52], q[54];
tdg q[54];
cx q[74], q[54];
t q[54];
cx q[74], q[52];
tdg q[52];
cx q[74], q[52];
t q[74];
t q[52];
h q[52];
t q[35];
s q[36];
cx q[0], q[11];
h q[73];
cx q[73], q[27];
t q[69];
h q[10];
cx q[10], q[9];
tdg q[9];
cx q[77], q[9];
t q[9];
cx q[10], q[9];
tdg q[9];
cx q[77], q[9];
t q[9];
cx q[77], q[10];
tdg q[10];
cx q[77], q[10];
t q[77];
t q[10];
h q[10];
cx q[85], q[8];
h q[29];
s q[20];
t q[82];
s q[22];
h q[62];
cx q[62], q[0];
tdg q[0];
cx q[82], q[0];
t q[0];
cx q[62], q[0];
tdg q[0];
cx q[82], q[0];
t q[0];
cx q[82], q[62];
tdg q[62];
cx q[82], q[62];
t q[82];
t q[62];
h q[62];
cx q[54], q[87];
h q[1];
h q[87];
cx q[87], q[9];
tdg q[9];
cx q[32], q[9];
t q[9];
cx q[87], q[9];
tdg q[9];
cx q[32], q[9];
t q[9];
cx q[32], q[87];
tdg q[87];
cx q[32], q[87];
t q[32];
t q[87];
h q[87];
h q[71];
cx q[61], q[45];
t q[53];
h q[5];
s q[5];
cx q[52], q[1];
cx q[24], q[77];
h q[64];
h q[73];
cx q[73], q[64];
tdg q[64];
cx q[53], q[64];
t q[64];
cx q[73], q[64];
tdg q[64];
cx q[53], q[64];
t q[64];
cx q[53], q[73];
tdg q[73];
cx q[53], q[73];
t q[53];
t q[73];
h q[73];
t q[42];
h q[86];
h q[54];
cx q[54], q[43];
tdg q[43];
cx q[41], q[43];
t q[43];
cx q[54], q[43];
tdg q[43];
cx q[41], q[43];
t q[43];
cx q[41], q[54];
tdg q[54];
cx q[41], q[54];
t q[41];
t q[54];
h q[54];
h q[53];
cx q[53], q[76];
tdg q[76];
cx q[7], q[76];
t q[76];
cx q[53], q[76];
tdg q[76];
cx q[7], q[76];
t q[76];
cx q[7], q[53];
tdg q[53];
cx q[7], q[53];
t q[7];
t q[53];
h q[53];
s q[22];
cx q[24], q[2];
s q[11];
s q[29];
s q[62];
t q[4];
cx q[67], q[60];
s q[58];
h q[1];
h q[21];
cx q[21], q[66];
tdg q[66];
cx q[41], q[66];
t q[66];
cx q[21], q[66];
tdg q[66];
cx q[41], q[66];
t q[66];
cx q[41], q[21];
tdg q[21];
cx q[41], q[21];
t q[41];
t q[21];
h q[21];
t q[84];
cx q[11], q[47];
cx q[66], q[69];
h q[71];
t q[13];
s q[0];
s q[16];
s q[57];
s q[60];
t q[47];
cx q[80], q[66];
s q[40];
cx q[64], q[81];
s q[85];
s q[23];
t q[1];
cx q[0], q[42];
s q[5];
t q[5];
h q[62];
s q[88];
h q[52];
h q[5];
s q[0];
cx q[45], q[66];
s q[37];
s q[9];
cx q[84], q[69];
t q[17];
s q[53];
h q[78];
cx q[78], q[39];
tdg q[39];
cx q[10], q[39];
t q[39];
cx q[78], q[39];
tdg q[39];
cx q[10], q[39];
t q[39];
cx q[10], q[78];
tdg q[78];
cx q[10], q[78];
t q[10];
t q[78];
h q[78];
cx q[55], q[52];
t q[77];
t q[4];
h q[32];
cx q[32], q[86];
tdg q[86];
cx q[74], q[86];
t q[86];
cx q[32], q[86];
tdg q[86];
cx q[74], q[86];
t q[86];
cx q[74], q[32];
tdg q[32];
cx q[74], q[32];
t q[74];
t q[32];
h q[32];
t q[79];
t q[4];
t q[46];
h q[48];
h q[79];
cx q[87], q[1];
h q[56];
cx q[86], q[26];
h q[13];
h q[47];
h q[1];
h q[22];
h q[45];
s q[20];
t q[30];
h q[27];
cx q[27], q[1];
tdg q[1];
cx q[18], q[1];
t q[1];
cx q[27], q[1];
tdg q[1];
cx q[18], q[1];
t q[1];
cx q[18], q[27];
tdg q[27];
cx q[18], q[27];
t q[18];
t q[27];
h q[27];
h q[35];
cx q[35], q[53];
tdg q[53];
cx q[44], q[53];
t q[53];
cx q[35], q[53];
tdg q[53];
cx q[44], q[53];
t q[53];
cx q[44], q[35];
tdg q[35];
cx q[44], q[35];
t q[44];
t q[35];
h q[35];
h q[50];
t q[81];
s q[3];
cx q[43], q[27];
cx q[20], q[15];
cx q[13], q[33];
h q[69];
cx q[69], q[85];
tdg q[85];
cx q[22], q[85];
t q[85];
cx q[69], q[85];
tdg q[85];
cx q[22], q[85];
t q[85];
cx q[22], q[69];
tdg q[69];
cx q[22], q[69];
t q[22];
t q[69];
h q[69];
h q[38];
cx q[38], q[85];
tdg q[85];
cx q[20], q[85];
t q[85];
cx q[38], q[85];
tdg q[85];
cx q[20], q[85];
t q[85];
cx q[20], q[38];
tdg q[38];
cx q[20], q[38];
t q[20];
t q[38];
h q[38];
h q[89];
cx q[76], q[40];
s q[80];
h q[77];
h q[76];
cx q[76], q[17];
tdg q[17];
cx q[46], q[17];
t q[17];
cx q[76], q[17];
tdg q[17];
cx q[46], q[17];
t q[17];
cx q[46], q[76];
tdg q[76];
cx q[46], q[76];
t q[46];
t q[76];
h q[76];
s q[76];
cx q[43], q[25];
t q[6];
s q[88];
cx q[70], q[6];
s q[86];
t q[66];
h q[43];
cx q[43], q[70];
tdg q[70];
cx q[62], q[70];
t q[70];
cx q[43], q[70];
tdg q[70];
cx q[62], q[70];
t q[70];
cx q[62], q[43];
tdg q[43];
cx q[62], q[43];
t q[62];
t q[43];
h q[43];
s q[30];
h q[34];
cx q[34], q[26];
tdg q[26];
cx q[42], q[26];
t q[26];
cx q[34], q[26];
tdg q[26];
cx q[42], q[26];
t q[26];
cx q[42], q[34];
tdg q[34];
cx q[42], q[34];
t q[42];
t q[34];
h q[34];
s q[77];
cx q[18], q[9];
h q[82];
s q[74];
t q[50];
h q[43];
t q[64];
h q[25];
cx q[81], q[41];
h q[61];
h q[15];
h q[4];
h q[66];
cx q[66], q[53];
tdg q[53];
cx q[22], q[53];
t q[53];
cx q[66], q[53];
tdg q[53];
cx q[22], q[53];
t q[53];
cx q[22], q[66];
tdg q[66];
cx q[22], q[66];
t q[22];
t q[66];
h q[66];
t q[53];
t q[71];
t q[82];
h q[10];
h q[80];
s q[12];
h q[20];
h q[49];
cx q[49], q[14];
tdg q[14];
cx q[76], q[14];
t q[14];
cx q[49], q[14];
tdg q[14];
cx q[76], q[14];
t q[14];
cx q[76], q[49];
tdg q[49];
cx q[76], q[49];
t q[76];
t q[49];
h q[49];
t q[5];
cx q[41], q[45];
t q[34];
h q[13];
cx q[13], q[42];
tdg q[42];
cx q[88], q[42];
t q[42];
cx q[13], q[42];
tdg q[42];
cx q[88], q[42];
t q[42];
cx q[88], q[13];
tdg q[13];
cx q[88], q[13];
t q[88];
t q[13];
h q[13];
cx q[66], q[45];
cx q[5], q[1];
t q[0];
h q[46];
cx q[46], q[39];
tdg q[39];
cx q[8], q[39];
t q[39];
cx q[46], q[39];
tdg q[39];
cx q[8], q[39];
t q[39];
cx q[8], q[46];
tdg q[46];
cx q[8], q[46];
t q[8];
t q[46];
h q[46];
t q[71];
t q[42];
h q[63];
cx q[63], q[39];
tdg q[39];
cx q[9], q[39];
t q[39];
cx q[63], q[39];
tdg q[39];
cx q[9], q[39];
t q[39];
cx q[9], q[63];
tdg q[63];
cx q[9], q[63];
t q[9];
t q[63];
h q[63];
h q[46];
cx q[46], q[2];
tdg q[2];
cx q[31], q[2];
t q[2];
cx q[46], q[2];
tdg q[2];
cx q[31], q[2];
t q[2];
cx q[31], q[46];
tdg q[46];
cx q[31], q[46];
t q[31];
t q[46];
h q[46];
s q[68];
cx q[63], q[74];
t q[48];
h q[13];
t q[1];
cx q[33], q[60];
h q[1];
s q[34];
h q[53];
cx q[53], q[36];
tdg q[36];
cx q[34], q[36];
t q[36];
cx q[53], q[36];
tdg q[36];
cx q[34], q[36];
t q[36];
cx q[34], q[53];
tdg q[53];
cx q[34], q[53];
t q[34];
t q[53];
h q[53];
t q[86];
h q[69];
cx q[0], q[1];
t q[3];
z q[4];
z q[5];
y q[6];
cx q[8], q[9];
s q[9];
cx q[8], q[9];
t q[8];
s q[10];
cx q[15], q[14];
t q[14];
cx q[15], q[14];
z q[16];
z q[17];
cx q[19], q[18];
s q[18];
cx q[19], q[18];
s q[18];
x q[20];
cx q[22], q[21];
cx q[24], q[23];
z q[25];
t q[26];
s q[27];
cx q[27], q[26];
t q[29];
cx q[28], q[29];
s q[31];
s q[30];
cx q[30], q[31];
s q[33];
t q[34];
s q[34];
cx q[33], q[34];
cx q[35], q[36];
s q[36];
cx q[35], q[36];
t q[35];
cx q[40], q[39];
y q[41];
cx q[43], q[42];
t q[42];
s q[42];
cx q[43], q[42];
z q[44];
h q[88];
t q[60];
t q[72];
t q[56];
cx q[49], q[51];
cx q[80], q[55];
h q[85];
h q[66];
cx q[52], q[82];
s q[71];
h q[51];
cx q[51], q[70];
tdg q[70];
cx q[86], q[70];
t q[70];
cx q[51], q[70];
tdg q[70];
cx q[86], q[70];
t q[70];
cx q[86], q[51];
tdg q[51];
cx q[86], q[51];
t q[86];
t q[51];
h q[51];
h q[71];
cx q[71], q[86];
tdg q[86];
cx q[80], q[86];
t q[86];
cx q[71], q[86];
tdg q[86];
cx q[80], q[86];
t q[86];
cx q[80], q[71];
tdg q[71];
cx q[80], q[71];
t q[80];
t q[71];
h q[71];
h q[84];
t q[74];
h q[80];
s q[73];
s q[64];
h q[59];
t q[89];
h q[51];
s q[80];
h q[69];
cx q[69], q[54];
tdg q[54];
cx q[71], q[54];
t q[54];
cx q[69], q[54];
tdg q[54];
cx q[71], q[54];
t q[54];
cx q[71], q[69];
tdg q[69];
cx q[71], q[69];
t q[71];
t q[69];
h q[69];
h q[81];
h q[77];
s q[86];
cx q[76], q[74];
cx q[56], q[47];
h q[65];
cx q[65], q[86];
tdg q[86];
cx q[71], q[86];
t q[86];
cx q[65], q[86];
tdg q[86];
cx q[71], q[86];
t q[86];
cx q[71], q[65];
tdg q[65];
cx q[71], q[65];
t q[71];
t q[65];
h q[65];
s q[86];
s q[78];
s q[47];
h q[76];
t q[77];
s q[52];
cx q[50], q[52];
t q[60];
s q[65];
h q[86];
cx q[69], q[81];
s q[84];
h q[78];
cx q[78], q[67];
tdg q[67];
cx q[70], q[67];
t q[67];
cx q[78], q[67];
tdg q[67];
cx q[70], q[67];
t q[67];
cx q[70], q[78];
tdg q[78];
cx q[70], q[78];
t q[70];
t q[78];
h q[78];
cx q[73], q[46];
cx q[60], q[59];
h q[77];
t q[48];