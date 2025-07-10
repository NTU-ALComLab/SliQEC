#ifndef _EQCHECKER_H_
#define _EQCHECKER_H_

#include "bddSystem.h"


enum class EqType 
{ 
    Feq,    // Full equivalence 
    SPeq,   // Strong Partial equivalence 
    SPeqS,  // Strong Partial equivalence Special case (Unweighted Zero-Clean-Ancilla)
    CWPeq   // Constant-probability Weak Partial equialence
};

class EquivalenceChecker : public BDDSystem
{
public:
    // Constructor and Destructor
    EquivalenceChecker
    (
        std::vector<std::vector<GateType>>& gates,
        std::vector<std::vector<std::vector<int>>>& qubits,
        int n,
        int nQd,
        int nQkc,
        int nQkd,
        int nQm, 
        int nQw,  
        int nQp,  
        int nQg, 
        int nQkr,
        std::vector<std::string>& careSet,
        std::vector<std::vector<std::pair<std::vector<std::string>, std::pair<int, int>>>>& weightFuns,
        bool isReorder,
        EqType eqType
    );

    ~EquivalenceChecker()  
    {
        clear();
    }

    void check();
    void printInfo(double runtime, size_t memPeak) const;

private:
    std::vector<std::vector<GateType>> _gates;              // gates in circuits. [nCircuit]*[#gate]
    std::vector<std::vector<std::vector<int>>> _qubits;     // ith qubits of gates in circuits. [nCircuit]*[#gates]*[#qubits]
    int _ratio;                                             // gate count ratio. |circuit2|/|circuit1|
    int _isEq;                                              // if the result is equivalent or not.
    bool _isGatesSwap;                                      // if circuit1 and circuit2 are swapped.
    EqType _eqType;                                         // the equivalence checking type (Feq/Peq/PeqS)  
    
    int _nQd;                                               // #data qubits. (n_d in the paper)
    int _nQkc;                                              // #clean ancilla qubits. (n_kc in the paper)
    int _nQkd;                                              // #dirty ancilla qubits. (n_kd in the paper)
    int _nQm;                                               // #measured qubits. (n_m in the paper)     
    int _nQw;                                               // #weighted qubits. (n_w in the paper)
    int _nQp;                                               // #propagating qubits. (n_p in the paper)
    int _nQg;                                               // #garbage qubits. (n_g in the paper)   
    int _nQkr;                                              // #reverted clean ancilla qubits. (n_kr in the paper) 
    int _nQextra;
    std::vector<std::string> _careSet;                                                                // care set of circuits
    std::vector<std::vector<std::pair<std::vector<std::string>, std::pair<int, int>>>> _weightFuns;   // weight functions of circuits  // <postfix Boolean function, <sign, power>>
    int _weightFunBias; 

    void invertCircuit(std::vector<GateType> &gate);
    void init();
    void applyGate(int ithCircuit, GateType type, std::vector<int> qubit, bool right);
    void calculateMiter();
    void checkFeq();
    DdNode* func2node(std::vector<std::string>& func, int bias);
    
    void buildMatrix(int ithCircuit);
    void clearCleanAnc(int ithCircuit, int start, int end);
    bool checkDAnc(int ithCircuit);
    bool checkRCAnc(int ithCircuit);
    void prepare_init(int ithCircuit);
    void extract(int ithCircuit);
    void checkPeq();
    void checkPeqS();
    void printResult() const;

    // Clean up EquivalenceChecker
    void clear() 
    {
        _gates[0].clear();
        _gates[1].clear();
        _gates.clear();
        _qubits[0].clear();
        _qubits[1].clear();
        _qubits.clear();
    };
};

#endif
