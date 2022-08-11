#ifndef _PEQCHECKER_H_
#define _PEQCHECKER_H_

#include "bddSystem.h"

class PartialEquivalenceChecker : public BDDSystem
{
public:
    // constructor and destructor
    PartialEquivalenceChecker(std::string circuit1Str, std::string circuit2Str, int nQubitIn, int nQibitOut, int r, bool isReorder);
    ~PartialEquivalenceChecker()  {
        clear();
    }

    void runPEC();
    void runPECSpecial();
    void printInfo(double runtime, size_t memPeak);

private:
    int _nInput;
    int _nOutput;
    std::vector<std::vector<GateType> > _gates; // 2 * #gates
    std::vector<std::vector<std::vector<int> > > _qubits; // 2 * #gates * #qubits
    int _ratio; // gatecount ratio: |G_p|/|G|
    bool _isPEC;

    void setupDD(bool special);
    void invertCircuit(std::vector<GateType> &gate);
    void initBaseBDD();
    void calculateMiter();
    void applyGate(int ithCircuit, GateType type, std::vector<int> &qubit, bool right);
    bool checkPEC();
    bool checkPECSpecial();
    void getResult(bool special);
    void extract(int ithCircuit);

    // Clean up EquivalenceChecker
    void clear() {
        _gates[0].clear();
        _gates[1].clear();
        _gates.clear();
        _qubits[0].clear();
        _qubits[1].clear();
        _qubits.clear();
    };
};

#endif
