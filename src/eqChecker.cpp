#include "eqChecker.h"

EquivalenceChecker::EquivalenceChecker
(
    std::vector<std::vector<GateType>>& gates,
    std::vector<std::vector<std::vector<int>>>& qubits,
    int n,
    int nQin, 
    int nQout, 
    int r, 
    bool isReorder,
    EqType eqType
)
:   BDDSystem
    ( 
        (eqType == EqType::Peq)? 2 : 1, // nCircuit
        r, 
        isReorder
    )
{
    _gates = gates;
    _qubits = qubits;
    _n = n;
    _nQin = nQin;
    _nQout = nQout;
    _eqType = eqType;

    // the longer circuit (in gatecount) is stored in gates[1]
    if (_gates[0].size() > _gates[1].size())
    {
        _gates[0].swap(_gates[1]);
        _qubits[0].swap(_qubits[1]);
    }

    // compute ratio
    _ratio = round(((double) _gates[1].size()) / ((double) _gates[0].size()));
}

/**Function*************************************************************

  Synopsis    [Run the checking procedure.]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void EquivalenceChecker::check()
{
    init();
    if(_eqType == EqType::Peq)
    {
        extract(0);
        extract(1);
        checkPeq();
    }
    else if(_eqType == EqType::PeqS)
    {
        calculateMiter();
        checkPeqS();
    }
    else if(_eqType == EqType::Feq)
    {
        calculateMiter();
        checkFeq();
    }
    printResult();
}

/**Function*************************************************************

  Synopsis    [Invert the gates in the given circuit]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/

void EquivalenceChecker::invertCircuit(std::vector<GateType> &gate)
{
    for (int i = 0; i < gate.size(); i++)
    {
        if (gate[i] == GateType::S) gate[i] = GateType::SDG;
        else if (gate[i] == GateType::SDG) gate[i] = GateType::S;
        else if (gate[i] == GateType::T) gate[i] = GateType::TDG;
        else if (gate[i] == GateType::TDG) gate[i] = GateType::T;
        else if (gate[i] == GateType::RX_PI_2) gate[i] = GateType::RX_PI_2_DG;
        else if (gate[i] == GateType::RX_PI_2_DG) gate[i] = GateType::RX_PI_2;
        else if (gate[i] == GateType::RY_PI_2) gate[i] = GateType::RY_PI_2_DG;
        else if (gate[i] == GateType::RY_PI_2_DG) gate[i] = GateType::RY_PI_2;
    }
}

/**Function*************************************************************

  Synopsis    [Initialize]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/

void EquivalenceChecker::init()
{
    if (_eqType == EqType::PeqS)
    {
        int nAncilla = _n - _nQin;
        if (_nQout > nAncilla)
        {
            _n += _nQout - nAncilla;
        }
    }

    ddInitialize();
    initIdentity();
}

/**Function*************************************************************

  Synopsis    [Apply a gate]

  Description [Apply a gate to the left side of the matrix if right = 0, the right side otherwise.]

  SideEffects []

  SeeAlso     []

***********************************************************************/
void EquivalenceChecker::applyGate(int ithCircuit, GateType type, std::vector<int> qubit, bool right)
{
    if (right) for (int i = 0; i < qubit.size(); i++) qubit[i] += _n;

    if (type == GateType::X) PauliX(ithCircuit, qubit[0]);
    else if (type == GateType::Y) PauliY(ithCircuit, qubit[0], right);
    else if (type == GateType::Z) PauliZ(ithCircuit, qubit);
    else if (type == GateType::H) Hadamard(ithCircuit, qubit[0]);
    else if (type == GateType::S) Phase_shift(ithCircuit, 2, qubit[0]);
    else if (type == GateType::SDG) Phase_shift_dagger(ithCircuit, -2, qubit[0]);
    else if (type == GateType::T) Phase_shift(ithCircuit, 4, qubit[0]);
    else if (type == GateType::TDG) Phase_shift_dagger(ithCircuit, -4, qubit[0]);
    else if (type == GateType::RX_PI_2) rx_pi_2(ithCircuit, qubit[0], false);
    else if (type == GateType::RX_PI_2_DG) rx_pi_2(ithCircuit, qubit[0], true);
    else if (type == GateType::RY_PI_2) ry_pi_2(ithCircuit, qubit[0], right^false);
    else if (type == GateType::RY_PI_2_DG) ry_pi_2(ithCircuit, qubit[0], right^true);
    else if (type == GateType::CX)
    {
        std::vector<int> ncont(0);
        int targ = qubit[1];
        qubit.pop_back();
        Toffoli(ithCircuit, targ, qubit, ncont);
        ncont.clear();
    }
    else if (type == GateType::CZ) PauliZ(ithCircuit, qubit);
    else if (type == GateType::SWAP)
    {
        std::vector<int> cont(0);
        Fredkin(ithCircuit, qubit[0], qubit[1], cont);
        cont.clear();
    }
    else if (type == GateType::CSWAP)
    {
        int swapA = qubit[1], swapB = qubit[2];
        qubit.pop_back();
        qubit.pop_back();
        Fredkin(ithCircuit, swapA, swapB, qubit);
    }
    else if (type == GateType::CCX)
    {
        std::vector<int> ncont(0);
        int targ = qubit.back();
        qubit.pop_back();
        Toffoli(ithCircuit, targ, qubit, ncont);
        ncont.clear();
    }

    if (_ddManager != NULL)
        updateNodeCount();
}

/**Function*************************************************************

  Synopsis    [Calculate the miter]

  Description [
               Apply gates in G and G' to evolve matrx from identity interleavingly.
               The longer gateuit (gates[1]) is seen as G (applied to the left side of the matrix)
               to save computation overhead
              ]

  SideEffects []

  SeeAlso     []

***********************************************************************/
void EquivalenceChecker::calculateMiter()
{
    int cntCir0 = 0, cntCir1 = 0;

    if (_isReorder) Cudd_AutodynEnable(_ddManager, CUDD_REORDER_SYMM_SIFT);

    invertCircuit(_gates[1]);
    while (cntCir0 < _gates[0].size() || cntCir1 < _gates[1].size())
    {
        // apply 1 gate from gates[0]
        if (cntCir0 < _gates[0].size())
        {
            applyGate(0, _gates[0][cntCir0], _qubits[0][cntCir0], false);
            cntCir0++;
        }
        // apply ratio gate(s) from gates[1]
        while(  cntCir1 * _gates[0].size() < cntCir0 * _gates[1].size()  &&  cntCir1 < _gates[1].size()   )  
        {
            applyGate(0, _gates[1][cntCir1], _qubits[1][cntCir1], true);
            cntCir1++;
        }
    }
    invertCircuit(_gates[1]);

    if (_isReorder) Cudd_AutodynDisable(_ddManager);
}

/**Function*************************************************************

  Synopsis    [Extract information]

  Description [Extract the information needed to compare between two circuits.]

  SideEffects []

  SeeAlso     []

***********************************************************************/

void EquivalenceChecker::extract(int ithCircuit){
    if (_isReorder) Cudd_AutodynEnable(_ddManager, CUDD_REORDER_SYMM_SIFT);

    for(int cntCir = 0; cntCir < _gates[ithCircuit].size(); cntCir ++){
        applyGate(ithCircuit, _gates[ithCircuit][cntCir], _qubits[ithCircuit][cntCir], 0);
    }
    // setup G

    for(int i = 0; i < _w; i++){
        for(int j = 0; j < _r; j++){
            for(int variable_1 = _nQin + _n; variable_1 < 2*_n; variable_1++){      // index of 1-variable
                DdNode *temp1, *temp2;

                temp1 = Cudd_Not(Cudd_bddIthVar(_ddManager, variable_1));
                Cudd_Ref(temp1);

                temp2 = Cudd_Cofactor(_ddManager, _allBDD[ithCircuit][i][j], temp1);
                Cudd_Ref(temp2);
                Cudd_RecursiveDeref(_ddManager, temp1);

                Cudd_RecursiveDeref(_ddManager, _allBDD[ithCircuit][i][j]);
                _allBDD[ithCircuit][i][j] = temp2;
            }
        }
    }
    // cofactor

    for(int i = 0; i < _w; i++){
        for(int j = 0; j < _r; j++){
            for(int index = 0; index < _nQout; index++){
                int variable_0 = index;                                    // index of 0-variable
                int variable_1 = _n + (_n - _nQout + index);                // index of 1-variable
                DdNode *temp1, *temp2, *temp3;

                temp1 = Cudd_bddXor(_ddManager, Cudd_bddIthVar(_ddManager, variable_0), Cudd_bddIthVar(_ddManager, variable_1));
                Cudd_Ref(temp1);

                temp2 = Cudd_Not(temp1);
                Cudd_Ref(temp2);
                Cudd_RecursiveDeref(_ddManager, temp1);

                temp3 = Cudd_bddAnd(_ddManager, _allBDD[ithCircuit][i][j], temp2);
                Cudd_Ref(temp3);
                Cudd_RecursiveDeref(_ddManager, temp2);

                Cudd_RecursiveDeref(_ddManager, _allBDD[ithCircuit][i][j]);
                _allBDD[ithCircuit][i][j] = temp3;
            }
        }
    }
    //

    for(int i = 0; i < _w; i++){
        for(int j = 0; j < _r; j++){
            for(int variable_1 = _nQin + _n; variable_1 < (_n - _nQout) + _n; variable_1++){     // index of 1-variable
                DdNode *temp1, *temp2;

                temp1 = Cudd_Not(Cudd_bddIthVar(_ddManager, variable_1));
                Cudd_Ref(temp1);

                temp2 = Cudd_bddAnd(_ddManager, _allBDD[ithCircuit][i][j], temp1);
                Cudd_Ref(temp2);
                Cudd_RecursiveDeref(_ddManager, temp1);

                Cudd_RecursiveDeref(_ddManager, _allBDD[ithCircuit][i][j]);
                _allBDD[ithCircuit][i][j] = temp2;
            }
        }
    }
    //

    invertCircuit(_gates[ithCircuit]);
    for(int cntCir = _gates[ithCircuit].size() - 1; cntCir >= 0 ; cntCir --){
        applyGate(ithCircuit, _gates[ithCircuit][cntCir], _qubits[ithCircuit][cntCir], 0);
    }
    invertCircuit(_gates[ithCircuit]);
    // apply G^-1

    for(int i = 0; i < _w; i++){
        for(int j = 0; j < _r; j++){
            for(int variable_0 = _nQin; variable_0 < _n; variable_0++){
                DdNode *temp1, *temp2;

                temp1 = Cudd_Not( Cudd_bddIthVar(_ddManager, variable_0) );
                Cudd_Ref(temp1);

                temp2 = Cudd_bddAnd(_ddManager, _allBDD[ithCircuit][i][j], temp1);
                Cudd_Ref(temp2);
                Cudd_RecursiveDeref(_ddManager, temp1);

                Cudd_RecursiveDeref(_ddManager, _allBDD[ithCircuit][i][j]);
                _allBDD[ithCircuit][i][j] = temp2;
            }
        }
    }
    //

    if (_isReorder) Cudd_AutodynDisable(_ddManager);
}

/**Function*************************************************************

  Synopsis    [Check if two circuits are equivalent]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void EquivalenceChecker::checkFeq()
{
    for (int i = 0; i < _w; i++)
    {
        for (int j = 0; j < _r; j++)
        {
            if(!(_allBDD[0][i][j] == _identityNode || _allBDD[0][i][j] == _zeroNode))
            {
                _isEq = 0;
                return;
            }
        }
    }
    
    _isEq = 1;
}

/**Function*************************************************************

  Synopsis    [Check if two circuits are partially equivalent]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/

void EquivalenceChecker::checkPeq()
{
    if ( (_k[0]-_k[1])%2 != 0){   // _k[0] - _k[1] must be even for two matrices to be equivalent
        assert(false);
        // this condition should not appear, because each gate will be applied pairwisely in (U^-1)U
    }

    int small, large, dk;
    if (_k[0] >= _k[1]){
        small = 1;
        large = 0;
        dk = (_k[0] - _k[1]) / 2;
    }
    else{
        small = 0;
        large = 1;
        dk = (_k[1] - _k[0]) / 2;
    }                               // need to consider different k

    for(int i = 0; i < _w; i++){
        for (int j = 0; j < dk; j++){
            if (_allBDD[small][i][_r - dk + j] != _allBDD[small][i][_r - dk - 1]){     // for 1's complement, the higher bits should be filled with the same bit
                _isEq = false;
                return;
            }
            if (_allBDD[large][i][j] != _zeroNode){
                _isEq = false;
                return;
            }
        }
        for (int j = 0; j < _r - dk; j++){
            if (_allBDD[small][i][j] != _allBDD[large][i][j + dk]){
                _isEq = false;
                return;
            }
        }
    }

    _isEq = true;
    return;
}

/**Function*************************************************************

  Synopsis    [Check if two circuits are partially equivalent in the special case.]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/

void EquivalenceChecker::checkPeqS()
{
    DdNode *mask = _zeroNode;
    Cudd_Ref(mask);
    for (int index = 0; index < _nQout; index++){
        int variable_0 = index;
        int variable_1 = _n + index;

        DdNode *temp1, *temp2;

        temp1 = Cudd_bddXor(_ddManager, Cudd_bddIthVar(_ddManager, variable_0), Cudd_bddIthVar(_ddManager, variable_1));
        Cudd_Ref(temp1);

        temp2 = Cudd_bddOr(_ddManager, mask, temp1);
        Cudd_RecursiveDeref(_ddManager, temp1);
        Cudd_Ref(temp2);

        Cudd_RecursiveDeref(_ddManager, mask);
        mask = temp2;
    }

    for(int i = 0; i < _w; i++){
        for(int j = 0; j < _r; j++){
            DdNode *temp;
            temp = Cudd_bddAnd(_ddManager, _allBDD[0][i][j], mask);
            Cudd_Ref(temp);

            if (temp != _zeroNode){
                Cudd_RecursiveDeref(_ddManager, mask);
                Cudd_RecursiveDeref(_ddManager, temp);
                _isEq = false;
                return;
            }
            Cudd_RecursiveDeref(_ddManager, temp);
        }
    }
    Cudd_RecursiveDeref(_ddManager, mask);
    _isEq = true;
    return;
}

/**Function*************************************************************

  Synopsis    [Print the final result]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void EquivalenceChecker::printResult() const
{
    std::cout << "  #Qubits: " << _n << std::endl;
    if(_eqType != EqType::Feq) std::cout << "  #Input Qubits: " << _nQin << std::endl;
    if(_eqType != EqType::Feq) std::cout << "  #Output Qubits: " << _nQout << std::endl;
    std::cout << "  Gatecount of circuit1 : " << _gates[0].size() << std::endl;
    std::cout << "  Gatecount of circuit2: " << _gates[1].size() << std::endl;
    printf("  |circuit2|/|circuit1|: %.2f\n", ((double) _gates[1].size()) / ((double) _gates[0].size()));
    if(_eqType == EqType::Feq) std::cout << "  Is equivalent? ";
    else std::cout << "  Is partially equivalent? ";
    if (_isEq) std::cout << "Yes" << std::endl;
    else std::cout << "No" << std::endl;
}

/**Function*************************************************************

  Synopsis    [Print statistics.]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void EquivalenceChecker::printInfo(double runtime, size_t memPeak) const
{
    std::cout << "  Runtime: " << runtime << " seconds" << std::endl;
    std::cout << "  Peak memory usage: " << memPeak << " bytes" << std::endl; //unit in bytes
    std::cout << "  Max #nodes: " << _nodeCount << std::endl;
    std::cout << "  Integer bit size: " << _r << std::endl;
}