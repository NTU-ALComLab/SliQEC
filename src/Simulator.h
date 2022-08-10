#ifndef _SIMULATOR_H_
#define _SIMULATOR_H_

#include <iostream>
#include <stdio.h> // FILE
#include <unordered_map>
#include <sys/time.h> //estimate time
#include <fstream> //fstream
#include <sstream> // int to string
#include <cstdlib> //atoi
#include <string> //string
#include <sstream>
#include <random>
#include <cmath>
#include <vector>
#include <algorithm>
#include "../cudd/cudd/cudd.h"
#include "../cudd/cudd/cuddInt.h"
#include "../cudd/util/util.h"


#define PI 3.14159265358979323846264338327950288419716939937510582097494459230781640628620899

class BDDbased
{
    friend class PartialEquivalenceChecker;
    friend class PartialEquivalenceCheckerSpecial;
public:
    // constructor and destructor
    BDDbased(int circuit_number, int bitSize, bool reorder) : nCircuit(circuit_number), n(0), r(bitSize), w(4), inc(3), NodeCount(0), gatecount(0), isReorder(reorder){}
    ~BDDbased()  {
        clear();
    }

    /* gates */
    void Toffoli(int ithCircuit, int targ, std::vector<int> cont, std::vector<int> ncont);
    void Fredkin(int ithCircuit, int swapA , int swapB, std::vector<int> cont);
    //void Peres(int a, int b, int c);
    //void Peres_i(int a, int b, int c);
    void Hadamard(int ithCircuit, int iqubit);
    void rx_pi_2(int ithCircuit, int iqubit, bool dagger);
    void ry_pi_2(int ithCircuit, int iqubit, bool tanspose);
    void Phase_shift(int ithCircuit, int phase, int iqubit); // phase can only be 2 to the power of an integer
    void Phase_shift_dagger(int ithCircuit, int phase, int iqubit);
    void PauliX(int ithCircuit, int iqubit);
    void PauliY(int ithCircuit, int iqubit, bool transpose);
    void PauliZ(int ithCircuit, std::vector<int> iqubit); // Z or CZ

    /* misc */
    void reorder();

private:
    DdManager *manager;
    DdNode ****All_BDD;     // [circuit_index][w=4][r]
    int *k; // k in algebraic representation
    int nCircuit; // # of circuits
    int n; // # of qubits
    int r; // resolution of integers
    int w; // # of integers
    int inc; // add inc BDDs when overflow occurs, used in alloc_BDD
    bool isReorder;

    unsigned long gatecount;
    unsigned long NodeCount;

    /* misc */
    void init_identity();
    void alloc_BDD(DdNode ***Bdd, bool extend);
    int overflow3(DdNode *g, DdNode *h, DdNode *crin);
    int overflow2(DdNode *g, DdNode *crin);
    void nodecount();

    // Clean up BDDbased
    void clear() {
        for (int i = 0; i < nCircuit; i++)
            for (int j = 0; j < w; j++)
                for (int k = 0; k < r; k++)
                    Cudd_RecursiveDeref(manager, All_BDD[i][j][k]);

        for (int i = 0; i < nCircuit; i++){
            for (int j = 0; j < w; j++){
                delete[] All_BDD[i][j];
            }
            delete[] All_BDD[i];
        }
        delete[] All_BDD;

        Cudd_Quit(manager);
    };
};

enum gateType {X, Y, Z, H, S, SDG, T, TDG, RX_PI_2, RX_PI_2_DG, RY_PI_2, RY_PI_2_DG, CX, CZ, CCX, SWAP, CSWAP};

class PartialEquivalenceChecker : public BDDbased
{
public:
    // constructor and destructor
    PartialEquivalenceChecker(std::string G, std::string G_p, int inputNumber, int outputNumber, int bitSize, bool reorder);
    ~PartialEquivalenceChecker()  {
        clear();
    }

    void runPEC();
    void runPEC_special();
    void print_info(double runtime, size_t memPeak);

private:
    int nInput;
    int nOutput;
    DdNode *zero;
    std::vector<std::vector<gateType> > gates; // 2 * #gates
    std::vector<std::vector<std::vector<int> > > qubits; // 2 * #gates * #qubits
    int ratio; // gatecount ratio: |G_p|/|G|
    bool isPEC;

    void setupDD(bool special);
    void circuitParser(std::string G, std::vector<gateType> &gate, std::vector<std::vector<int> > &qubit);
    void invertCircuit(std::vector<gateType> &gate);
    void initBaseBDD();
    void sim_miter();
    void applyGate(int ithCircuit, gateType type, std::vector<int> qubit, bool right);
    bool checkPEC();
    bool checkPEC_special();
    void getResult(bool special);
    void extract(int ithCircuit);

    // Clean up EquivalenceChecker
    void clear() {
        gates[0].clear();
        gates[1].clear();
        gates.clear();
        qubits[0].clear();
        qubits[1].clear();
        qubits.clear();
    };
};

#endif
