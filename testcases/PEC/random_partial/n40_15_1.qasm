OPENQASM 2.0;
include "qelib1.inc";
qreg q[40];
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
t q[13];
cx q[11], q[24];
t q[31];
t q[11];
t q[18];
s q[28];
s q[9];
h q[15];
cx q[15], q[31];
tdg q[31];
cx q[11], q[31];
t q[31];
cx q[15], q[31];
tdg q[31];
cx q[11], q[31];
t q[31];
cx q[11], q[15];
tdg q[15];
cx q[11], q[15];
t q[11];
t q[15];
h q[15];
cx q[6], q[16];
h q[30];
cx q[30], q[11];
tdg q[11];
cx q[0], q[11];
t q[11];
cx q[30], q[11];
tdg q[11];
cx q[0], q[11];
t q[11];
cx q[0], q[30];
tdg q[30];
cx q[0], q[30];
t q[0];
t q[30];
h q[30];
s q[30];
t q[24];
cx q[19], q[37];
cx q[30], q[31];
h q[34];
t q[37];
h q[5];
cx q[5], q[4];
tdg q[4];
cx q[2], q[4];
t q[4];
cx q[5], q[4];
tdg q[4];
cx q[2], q[4];
t q[4];
cx q[2], q[5];
tdg q[5];
cx q[2], q[5];
t q[2];
t q[5];
h q[5];
h q[20];
h q[24];
s q[38];
h q[37];
s q[3];
cx q[21], q[38];
h q[32];
h q[16];
cx q[16], q[13];
tdg q[13];
cx q[14], q[13];
t q[13];
cx q[16], q[13];
tdg q[13];
cx q[14], q[13];
t q[13];
cx q[14], q[16];
tdg q[16];
cx q[14], q[16];
t q[14];
t q[16];
h q[16];
t q[30];
t q[31];
h q[8];
t q[36];
h q[9];
cx q[9], q[36];
tdg q[36];
cx q[11], q[36];
t q[36];
cx q[9], q[36];
tdg q[36];
cx q[11], q[36];
t q[36];
cx q[11], q[9];
tdg q[9];
cx q[11], q[9];
t q[11];
t q[9];
h q[9];
t q[30];
s q[8];
cx q[24], q[1];
t q[28];
t q[21];
h q[31];
h q[4];
s q[3];
h q[5];
t q[1];
t q[8];
s q[8];
t q[38];
h q[8];
cx q[8], q[28];
tdg q[28];
cx q[35], q[28];
t q[28];
cx q[8], q[28];
tdg q[28];
cx q[35], q[28];
t q[28];
cx q[35], q[8];
tdg q[8];
cx q[35], q[8];
t q[35];
t q[8];
h q[8];
t q[30];
cx q[13], q[2];
h q[28];
cx q[28], q[11];
tdg q[11];
cx q[16], q[11];
t q[11];
cx q[28], q[11];
tdg q[11];
cx q[16], q[11];
t q[11];
cx q[16], q[28];
tdg q[28];
cx q[16], q[28];
t q[16];
t q[28];
h q[28];
h q[3];
cx q[3], q[20];
tdg q[20];
cx q[17], q[20];
t q[20];
cx q[3], q[20];
tdg q[20];
cx q[17], q[20];
t q[20];
cx q[17], q[3];
tdg q[3];
cx q[17], q[3];
t q[17];
t q[3];
h q[3];
t q[6];
s q[36];
cx q[1], q[18];
h q[39];
cx q[39], q[19];
tdg q[19];
cx q[17], q[19];
t q[19];
cx q[39], q[19];
tdg q[19];
cx q[17], q[19];
t q[19];
cx q[17], q[39];
tdg q[39];
cx q[17], q[39];
t q[17];
t q[39];
h q[39];
h q[29];
h q[5];
cx q[5], q[19];
tdg q[19];
cx q[26], q[19];
t q[19];
cx q[5], q[19];
tdg q[19];
cx q[26], q[19];
t q[19];
cx q[26], q[5];
tdg q[5];
cx q[26], q[5];
t q[26];
t q[5];
h q[5];
s q[30];
h q[8];
cx q[8], q[3];
tdg q[3];
cx q[4], q[3];
t q[3];
cx q[8], q[3];
tdg q[3];
cx q[4], q[3];
t q[3];
cx q[4], q[8];
tdg q[8];
cx q[4], q[8];
t q[4];
t q[8];
h q[8];
h q[38];
h q[31];
cx q[31], q[16];
tdg q[16];
cx q[3], q[16];
t q[16];
cx q[31], q[16];
tdg q[16];
cx q[3], q[16];
t q[16];
cx q[3], q[31];
tdg q[31];
cx q[3], q[31];
t q[3];
t q[31];
h q[31];
h q[4];
cx q[4], q[26];
tdg q[26];
cx q[34], q[26];
t q[26];
cx q[4], q[26];
tdg q[26];
cx q[34], q[26];
t q[26];
cx q[34], q[4];
tdg q[4];
cx q[34], q[4];
t q[34];
t q[4];
h q[4];
cx q[36], q[27];
cx q[35], q[24];
h q[2];
cx q[2], q[21];
tdg q[21];
cx q[23], q[21];
t q[21];
cx q[2], q[21];
tdg q[21];
cx q[23], q[21];
t q[21];
cx q[23], q[2];
tdg q[2];
cx q[23], q[2];
t q[23];
t q[2];
h q[2];
cx q[37], q[8];
s q[19];
cx q[4], q[7];
cx q[3], q[17];
t q[38];
s q[2];
h q[1];
cx q[1], q[31];
tdg q[31];
cx q[0], q[31];
t q[31];
cx q[1], q[31];
tdg q[31];
cx q[0], q[31];
t q[31];
cx q[0], q[1];
tdg q[1];
cx q[0], q[1];
t q[0];
t q[1];
h q[1];
t q[28];
h q[25];
cx q[25], q[35];
tdg q[35];
cx q[28], q[35];
t q[35];
cx q[25], q[35];
tdg q[35];
cx q[28], q[35];
t q[35];
cx q[28], q[25];
tdg q[25];
cx q[28], q[25];
t q[28];
t q[25];
h q[25];
h q[0];
h q[32];
h q[9];
cx q[34], q[23];
s q[5];
h q[22];
cx q[22], q[27];
tdg q[27];
cx q[17], q[27];
t q[27];
cx q[22], q[27];
tdg q[27];
cx q[17], q[27];
t q[27];
cx q[17], q[22];
tdg q[22];
cx q[17], q[22];
t q[17];
t q[22];
h q[22];
t q[9];
cx q[0], q[37];
h q[34];
t q[32];
h q[4];
cx q[4], q[39];
tdg q[39];
cx q[33], q[39];
t q[39];
cx q[4], q[39];
tdg q[39];
cx q[33], q[39];
t q[39];
cx q[33], q[4];
tdg q[4];
cx q[33], q[4];
t q[33];
t q[4];
h q[4];
s q[39];
h q[17];
s q[18];
h q[14];
cx q[14], q[16];
tdg q[16];
cx q[29], q[16];
t q[16];
cx q[14], q[16];
tdg q[16];
cx q[29], q[16];
t q[16];
cx q[29], q[14];
tdg q[14];
cx q[29], q[14];
t q[29];
t q[14];
h q[14];
h q[38];
cx q[38], q[8];
tdg q[8];
cx q[29], q[8];
t q[8];
cx q[38], q[8];
tdg q[8];
cx q[29], q[8];
t q[8];
cx q[29], q[38];
tdg q[38];
cx q[29], q[38];
t q[29];
t q[38];
h q[38];
h q[3];
cx q[3], q[10];
tdg q[10];
cx q[38], q[10];
t q[10];
cx q[3], q[10];
tdg q[10];
cx q[38], q[10];
t q[10];
cx q[38], q[3];
tdg q[3];
cx q[38], q[3];
t q[38];
t q[3];
h q[3];
h q[17];
cx q[17], q[12];
tdg q[12];
cx q[15], q[12];
t q[12];
cx q[17], q[12];
tdg q[12];
cx q[15], q[12];
t q[12];
cx q[15], q[17];
tdg q[17];
cx q[15], q[17];
t q[15];
t q[17];
h q[17];
t q[22];
cx q[23], q[12];
t q[36];
t q[19];
s q[30];
t q[19];
s q[29];
cx q[4], q[30];
cx q[26], q[10];
cx q[31], q[36];
cx q[38], q[34];
s q[23];
h q[3];
h q[14];
h q[26];
cx q[26], q[8];
tdg q[8];
cx q[25], q[8];
t q[8];
cx q[26], q[8];
tdg q[8];
cx q[25], q[8];
t q[8];
cx q[25], q[26];
tdg q[26];
cx q[25], q[26];
t q[25];
t q[26];
h q[26];
h q[19];
t q[35];
s q[22];
t q[16];
t q[17];
h q[11];
cx q[11], q[25];
tdg q[25];
cx q[18], q[25];
t q[25];
cx q[11], q[25];
tdg q[25];
cx q[18], q[25];
t q[25];
cx q[18], q[11];
tdg q[11];
cx q[18], q[11];
t q[18];
t q[11];
h q[11];
h q[32];
cx q[32], q[8];
tdg q[8];
cx q[2], q[8];
t q[8];
cx q[32], q[8];
tdg q[8];
cx q[2], q[8];
t q[8];
cx q[2], q[32];
tdg q[32];
cx q[2], q[32];
t q[2];
t q[32];
h q[32];
cx q[23], q[33];
t q[32];
s q[4];
t q[2];
cx q[27], q[32];
h q[13];
cx q[25], q[12];
t q[39];
h q[0];
cx q[0], q[6];
tdg q[6];
cx q[28], q[6];
t q[6];
cx q[0], q[6];
tdg q[6];
cx q[28], q[6];
t q[6];
cx q[28], q[0];
tdg q[0];
cx q[28], q[0];
t q[28];
t q[0];
h q[0];
x q[0];
y q[3];
z q[4];
t q[5];
cx q[6], q[5];
s q[6];
t q[8];
cx q[9], q[8];
y q[10];
cx q[11], q[12];
s q[11];
t q[12];
cx q[11], q[12];
t q[13];
s q[14];
cx q[13], q[14];
x q[15];
s q[17];
t q[18];
cx q[17], q[18];
t q[18];
s q[36];
h q[34];
h q[39];
h q[29];
cx q[29], q[35];
tdg q[35];
cx q[28], q[35];
t q[35];
cx q[29], q[35];
tdg q[35];
cx q[28], q[35];
t q[35];
cx q[28], q[29];
tdg q[29];
cx q[28], q[29];
t q[28];
t q[29];
h q[29];
h q[31];
cx q[31], q[34];
tdg q[34];
cx q[30], q[34];
t q[34];
cx q[31], q[34];
tdg q[34];
cx q[30], q[34];
t q[34];
cx q[30], q[31];
tdg q[31];
cx q[30], q[31];
t q[30];
t q[31];
h q[31];
s q[26];
t q[37];
t q[32];
h q[31];
cx q[31], q[33];
tdg q[33];
cx q[36], q[33];
t q[33];
cx q[31], q[33];
tdg q[33];
cx q[36], q[33];
t q[33];
cx q[36], q[31];
tdg q[31];
cx q[36], q[31];
t q[36];
t q[31];
h q[31];
h q[26];
h q[35];
h q[30];
cx q[30], q[35];
tdg q[35];
cx q[36], q[35];
t q[35];
cx q[30], q[35];
tdg q[35];
cx q[36], q[35];
t q[35];
cx q[36], q[30];
tdg q[30];
cx q[36], q[30];
t q[36];
t q[30];
h q[30];
t q[35];
h q[27];
t q[23];
h q[37];
cx q[30], q[38];
t q[27];
h q[23];
s q[33];