#include <stack>
#include "eqChecker.h"

// Constructor
EquivalenceChecker::EquivalenceChecker
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
)
:   BDDSystem
    ( 
        (eqType == EqType::Feq) ? 1 : 2, // nCircuit
        isReorder
    )
{
    _gates = gates;
    _qubits = qubits;
    _n = n;
    _nQd  = nQd;
    _nQkc = nQkc;
    _nQkd = nQkd;
    _nQm  = nQm;
    _nQw  = nQw;
    _nQp  = nQp;
    _nQg  = nQg;
    _nQkr = nQkr;
    _careSet = careSet;
    _weightFuns = weightFuns;
    _eqType = eqType;
    _isGatesSwap = false;
    
    // the circuit with larger gatecount is stored in gates[1]
    if (_gates[0].size() > _gates[1].size())
    {
        _isGatesSwap = true;
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
    if(_eqType == EqType::Feq)
    {
        calculateMiter();
        checkFeq();
    }
    else
    {
        std::vector<int> _io_backup = {_nQkc, _nQkd, _nQkr, _nQg};
        
        if(_eqType != EqType::SPeqS || _nQkd != 0) 
        {
            buildMatrix(0);
            buildMatrix(1);
            clearCleanAnc(0, _nQd,  _nQd + _nQkc);
            clearCleanAnc(1, _nQd,  _nQd + _nQkc);
            
            // check dirty ancilla qubits & view dirty ancilla qubits as reverted clean ancilla qubits
            if (!checkDAnc(0) || !checkDAnc(1))
            {
                _isEq = -1;
                printResult();
                return;
            }
            clearCleanAnc(0, _nQd + _nQkc,  _nQd + _nQkc + _nQkd);
            clearCleanAnc(1, _nQd + _nQkc,  _nQd + _nQkc + _nQkd);
            _nQkc = _nQkc + _nQkd;
            _nQkr = _nQkr + _nQkd;
            _nQkd = 0;
            
            // check reverted ancilla qubits & view reverted ancilla qubits as garbage qubits
            if (!checkRCAnc(0) || !checkRCAnc(1))
            {
                _isEq = -1;
                _nQkc = _io_backup[0];
                _nQkd = _io_backup[1];
                _nQkr = _io_backup[2];
                printResult(); 
                return;
            }
            _nQg  = _nQg + _nQkr;
            _nQkr = 0;
            
            // add extra column variables
            if (_nQkc > _nQm + _nQw + 2*_nQp)
            {
                _nQextra = 0;
            }
            else
            {
                _nQextra = _nQm + _nQw + 2*_nQp - _nQkc;
                for (int i = 0; i < _nQextra; ++i)
                {
                    Cudd_bddNewVar(_ddManager);
                }
            }
            
            // add extra copy variables
            if (_eqType == EqType::CWPeq)
            {
                for (int i = 0; i < _nQm; ++i)
                {
                    Cudd_bddNewVar(_ddManager);
                }
            }
        }
        
        // main algorithms for SPeq and CWPeq
        if(_eqType == EqType::SPeq || _eqType == EqType::CWPeq)
        {
            if (_eqType == EqType::CWPeq)    // exchange BDDs
            {
                for(int i = 0; i < _w; ++i)
                {
                    for(int j = 0; j < _r; ++j)
                    {
                        DdNode *temp = _allBDD[0][i][j];
                        _allBDD[0][i][j] = _allBDD[1][i][j];
                        _allBDD[1][i][j] = temp;
                    }
                }
                int temp = _k[0];
                _k[0] = _k[1];
                _k[1] = temp;
            }
            
            extract(0);
            extract(1);
            
            checkPeq();
        }
        else if(_eqType == EqType::SPeqS)
        {
            BDDSystem::clear();
            init(); // re-initialize
            calculateMiter();
            checkPeqS();
        }
        
        _nQkc = _io_backup[0];
        _nQkd = _io_backup[1];
        _nQkr = _io_backup[2];
        _nQg  = _io_backup[3];
    }
    printResult();
}

/**Function*************************************************************

  Synopsis    [Invert the gates in the given circuit.]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/

void EquivalenceChecker::invertCircuit(std::vector<GateType> &gate)
{
    for (int i = 0; i < gate.size(); i++)
    {
        if (gate[i] == GateType::MCS) gate[i] = GateType::MCSDG;
        else if (gate[i] == GateType::MCSDG) gate[i] = GateType::MCS;
        else if (gate[i] == GateType::MCT) gate[i] = GateType::MCTDG;
        else if (gate[i] == GateType::MCTDG) gate[i] = GateType::MCT;
        else if (gate[i] == GateType::RX_PI_2) gate[i] = GateType::RX_PI_2_DG;
        else if (gate[i] == GateType::RX_PI_2_DG) gate[i] = GateType::RX_PI_2;
        else if (gate[i] == GateType::RY_PI_2) gate[i] = GateType::RY_PI_2_DG;
        else if (gate[i] == GateType::RY_PI_2_DG) gate[i] = GateType::RY_PI_2;
    }
}

/**Function*************************************************************

  Synopsis    [Initialize checker.]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/

void EquivalenceChecker::init()
{
    ddInitialize();
    initIdentity();
}

/**Function*************************************************************

  Synopsis    [Apply a gate.]

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
    else if (type == GateType::MCS) 
    {
        std::vector<int> ncont(0);
        Controlled_Phase_shift(ithCircuit, 2, qubit, ncont);
    }
    else if (type == GateType::MCSDG)
    {
        std::vector<int> ncont(0);
        Controlled_Phase_shift_dagger(ithCircuit, -2, qubit, ncont);
    }
    else if (type == GateType::MCT) 
    {
        std::vector<int> ncont(0);
        Controlled_Phase_shift(ithCircuit, 4, qubit, ncont); 
    }
    else if (type == GateType::MCTDG)
    {
        std::vector<int> ncont(0);
        Controlled_Phase_shift_dagger(ithCircuit, -4, qubit, ncont);
    }
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
               Apply gates in circuit1 and circuit2 to evolve matrx from identity interleavingly.
               The longer gateuit (gates[1]) is seen as circuit1 (applied to the left side of the matrix)
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

  Synopsis    [Build the matrix of a circuit.]

  Description []

  SideEffects [Unnecessary columns are cleared. 
               The built matrix is Rd(U) in the paper. ]

  SeeAlso     []

***********************************************************************/

void EquivalenceChecker::buildMatrix(int ithCircuit){
    if (_isReorder) Cudd_AutodynEnable(_ddManager, CUDD_REORDER_SYMM_SIFT);

    for(int cntCir = 0; cntCir < _gates[ithCircuit].size(); cntCir ++)
        applyGate(ithCircuit, _gates[ithCircuit][cntCir], _qubits[ithCircuit][cntCir], 0);
        
    if (_isReorder) Cudd_AutodynDisable(_ddManager);
}

/**Function*************************************************************

  Synopsis    [Clear columns with clean ancilla != 0.]

  Description ["end" is exclusive]

  SideEffects []

  SeeAlso     []

***********************************************************************/

void EquivalenceChecker::clearCleanAnc(int ithCircuit, int start, int end){
    if (_isReorder) Cudd_AutodynEnable(_ddManager, CUDD_REORDER_SYMM_SIFT);
    
    DdNode *mask = Cudd_Not(_zeroNode);
    Cudd_Ref(mask);
    for(int ith_var = start; ith_var < end; ++ith_var)
    {        
        DdNode *column_variable = Cudd_bddIthVar(_ddManager, ith_var + _n);
        Cudd_Ref(column_variable);
        DdNode *temp1, *temp2;

        temp1 = Cudd_Not(column_variable);
        Cudd_Ref(temp1);
        Cudd_RecursiveDeref(_ddManager, column_variable);

        temp2 = Cudd_bddAnd(_ddManager, mask, temp1);
        Cudd_Ref(temp2);
        Cudd_RecursiveDeref(_ddManager, temp1);

        Cudd_RecursiveDeref(_ddManager, mask);
        mask = temp2;
    }
    
    for(int i = 0; i < _w; ++i)
    {
        for(int j = 0; j < _r; ++j)
        {
            DdNode *temp;
            
            temp = Cudd_bddAnd(_ddManager, _allBDD[ithCircuit][i][j], mask);
            Cudd_Ref(temp);
            Cudd_RecursiveDeref(_ddManager, _allBDD[ithCircuit][i][j]);
            
            _allBDD[ithCircuit][i][j] = temp;
        }
    }
    Cudd_RecursiveDeref(_ddManager, mask);
    
    if (_isReorder) Cudd_AutodynDisable(_ddManager);
}

/**Function*************************************************************

  Synopsis    [Check whether all dirty ancilla qubits are always returned.]

  Description [Dirty ancilla qubits should return to initial state.       ]

  SideEffects []

  SeeAlso     []

***********************************************************************/

bool EquivalenceChecker::checkDAnc(int ithCircuit){
    if (_isReorder) Cudd_AutodynEnable(_ddManager, CUDD_REORDER_SYMM_SIFT);
    
    for(int i = 0; i < _w; ++i)
    {
        for(int j = 0; j < _r; ++j)
        {
            // main
            for (int ith_var = _n - 1; ith_var >= _n - _nQkd; --ith_var)
            {
                DdNode *row_variable = Cudd_bddIthVar(_ddManager, ith_var);
                Cudd_Ref(row_variable);
                DdNode *column_variable = Cudd_bddIthVar(_ddManager, ith_var + _n);
                Cudd_Ref(column_variable);
                
                DdNode *temp_0, *temp_1;
                temp_0 = Cudd_Cofactor(_ddManager, _allBDD[ithCircuit][i][j], Cudd_Not(row_variable));
                Cudd_Ref(temp_0);
                temp_1 = Cudd_Cofactor(_ddManager, _allBDD[ithCircuit][i][j], row_variable);
                Cudd_Ref(temp_1);
                
                DdNode *temp_01, *temp_10;
                temp_01 = Cudd_Cofactor(_ddManager, temp_0, column_variable);
                Cudd_Ref(temp_01);
                temp_10 = Cudd_Cofactor(_ddManager, temp_1, Cudd_Not(column_variable));
                Cudd_Ref(temp_10);
                if (temp_01 != _zeroNode || temp_10 != _zeroNode)
                {
                    Cudd_RecursiveDeref(_ddManager, temp_01);
                    Cudd_RecursiveDeref(_ddManager, temp_10);
                    Cudd_RecursiveDeref(_ddManager, row_variable);
                    Cudd_RecursiveDeref(_ddManager, column_variable);
                    Cudd_RecursiveDeref(_ddManager, temp_0);
                    Cudd_RecursiveDeref(_ddManager, temp_1);
                    return false;
                }
                Cudd_RecursiveDeref(_ddManager, temp_01);
                Cudd_RecursiveDeref(_ddManager, temp_10);
                    
                DdNode *temp_00, *temp_11;
                temp_00 = Cudd_Cofactor(_ddManager, temp_0, Cudd_Not(column_variable));
                Cudd_Ref(temp_00);
                temp_11 = Cudd_Cofactor(_ddManager, temp_1, column_variable);
                Cudd_Ref(temp_11);
                if (temp_00 != temp_11)
                {
                    Cudd_RecursiveDeref(_ddManager, temp_00);
                    Cudd_RecursiveDeref(_ddManager, temp_11);
                    Cudd_RecursiveDeref(_ddManager, row_variable);
                    Cudd_RecursiveDeref(_ddManager, column_variable);
                    Cudd_RecursiveDeref(_ddManager, temp_0);
                    Cudd_RecursiveDeref(_ddManager, temp_1);
                    return false;
                }
                Cudd_RecursiveDeref(_ddManager, temp_00);
                Cudd_RecursiveDeref(_ddManager, temp_11);
                
                Cudd_RecursiveDeref(_ddManager, row_variable);
                Cudd_RecursiveDeref(_ddManager, column_variable);
                Cudd_RecursiveDeref(_ddManager, temp_0);
                    Cudd_RecursiveDeref(_ddManager, temp_1);
            }
        }
    }
        
    if (_isReorder) Cudd_AutodynDisable(_ddManager);
    return true;
}

/**Function*************************************************************

  Synopsis    [Check whether all reverted clean ancilla qubits are always returned.]

  Description [Reverted clean ancilla qubits should return to |0>.                 ]

  SideEffects []

  SeeAlso     []

***********************************************************************/

bool EquivalenceChecker::checkRCAnc(int ithCircuit){
    if (_isReorder) Cudd_AutodynEnable(_ddManager, CUDD_REORDER_SYMM_SIFT);
    
    // check dirty ancilla qubits
    for(int i = 0; i < _w; ++i)
    {
        for(int j = 0; j < _r; ++j)
        {
            // main
            DdNode *copy = _allBDD[ithCircuit][i][j];
            Cudd_Ref(copy);
            
            for (int ith_var = _n - 1; ith_var >= _n - _nQkr; --ith_var)
            {
                DdNode *row_variable = Cudd_bddIthVar(_ddManager, ith_var);
                Cudd_Ref(row_variable);
                DdNode *column_variable = Cudd_bddIthVar(_ddManager, ith_var + _n);
                Cudd_Ref(column_variable);
                    
                DdNode *temp_00;
                temp_00 = Cudd_bddAnd(_ddManager, Cudd_bddAnd(_ddManager, _allBDD[ithCircuit][i][j], Cudd_Not(row_variable)), Cudd_Not(column_variable));
                Cudd_Ref(temp_00);
                Cudd_RecursiveDeref(_ddManager, copy);
                copy = temp_00;
                
                Cudd_RecursiveDeref(_ddManager, row_variable);
                Cudd_RecursiveDeref(_ddManager, column_variable);
            }
            
            if (copy != _allBDD[ithCircuit][i][j])
            {
                Cudd_RecursiveDeref(_ddManager, copy);
                return false;
            }
                
            Cudd_RecursiveDeref(_ddManager, copy);
        }
    }
        
    if (_isReorder) Cudd_AutodynDisable(_ddManager);
    return true;
}


/**Function*************************************************************

  Synopsis    [Extra preparation steps for CWPeq.]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/

void EquivalenceChecker::prepare_init(int ithCircuit){
    if (_isReorder) Cudd_AutodynEnable(_ddManager, CUDD_REORDER_SYMM_SIFT);
        
    DdNode *mask;
    if (_careSet.empty())
    {
        mask = Cudd_Not(_zeroNode);
        Cudd_Ref(mask);
    }
    else
    {
        mask = func2node(_careSet, _n); // already referenced
    }
    
    for(int i = 0; i < _w; ++i)
    {
        for(int j = 0; j < _r; ++j)
        {
            // build C(1-i)|0>  // Note that the BDDs have been swapped
            for(int ith_var = 0; ith_var < _n; ++ith_var)
            {        
                DdNode *column_variable = Cudd_bddIthVar(_ddManager, ith_var + _n);
                Cudd_Ref(column_variable);
                
                DdNode *temp = Cudd_Cofactor(_ddManager, _allBDD[ithCircuit][i][j], Cudd_Not(column_variable));
                Cudd_Ref(temp);
                Cudd_RecursiveDeref(_ddManager, _allBDD[ithCircuit][i][j]);
                Cudd_RecursiveDeref(_ddManager, column_variable);
                
                _allBDD[ithCircuit][i][j] = temp;
            }
            
            // F <- F AND [X]=[Y]
            for(int ith_var = 0; ith_var < _nQm + _nQw; ++ith_var)
            {        
                DdNode *row_variable = Cudd_bddIthVar(_ddManager, ith_var);
                Cudd_Ref(row_variable);
                DdNode *column_variable = Cudd_bddIthVar(_ddManager, ith_var + _n);
                Cudd_Ref(column_variable);
                
                DdNode *temp1 = Cudd_bddXnor(_ddManager, row_variable, column_variable);
                Cudd_Ref(temp1);
                Cudd_RecursiveDeref(_ddManager, row_variable);
                Cudd_RecursiveDeref(_ddManager, column_variable);
                
                DdNode *temp2 = Cudd_bddAnd(_ddManager, _allBDD[ithCircuit][i][j], temp1);
                Cudd_Ref(temp2);
                Cudd_RecursiveDeref(_ddManager, temp1);
                Cudd_RecursiveDeref(_ddManager, _allBDD[ithCircuit][i][j]);
                
                _allBDD[ithCircuit][i][j] = temp2;
            }
            
            // F <- F AND [Y] \in S
            DdNode *temp_mask = Cudd_bddAnd(_ddManager, _allBDD[ithCircuit][i][j], mask);
            Cudd_Ref(temp_mask);
            Cudd_RecursiveDeref(_ddManager, _allBDD[ithCircuit][i][j]);
            _allBDD[ithCircuit][i][j] = temp_mask; 
        }
    }
        
    // Apply C(1-i)^-1
    invertCircuit(_gates[1 - ithCircuit]);
    for(int cntCir = _gates[1 - ithCircuit].size() - 1; cntCir >= 0 ; cntCir --)
        applyGate(ithCircuit, _gates[1 - ithCircuit][cntCir], _qubits[1 - ithCircuit][cntCir], 0);
    invertCircuit(_gates[1 - ithCircuit]);
        
    for(int i = 0; i < _w; ++i)
    {
        for(int j = 0; j < _r; ++j)
        {            
            // F <- F|xi
            for(int ith_var = 0; ith_var < _n; ++ith_var)
            {        
                DdNode *row_variable = Cudd_bddIthVar(_ddManager, ith_var);
                Cudd_Ref(row_variable);
                
                DdNode *temp = Cudd_Cofactor(_ddManager, _allBDD[ithCircuit][i][j], Cudd_Not(row_variable));
                Cudd_Ref(temp);
                Cudd_RecursiveDeref(_ddManager, row_variable);
                Cudd_RecursiveDeref(_ddManager, _allBDD[ithCircuit][i][j]);
                
                _allBDD[ithCircuit][i][j] = temp;
            }
        }
    }
  
    Cudd_RecursiveDeref(_ddManager, mask);
        
    // weighted 
    for (std::pair<std::vector<std::string>, std::pair<int, int>>& item : _weightFuns[1 - ithCircuit])
    {
        DdNode *func = func2node(item.first, _n); // already referenced
        int sign = item.second.first;    // 0 for 0, 1 for +, -1 for -
        int power = item.second.second;  // i for 2^i
        
        if (sign == 0)
        {  
            for(int i = 0; i < _w; ++i)
            {
                for(int j = 0; j < _r; ++j)
                {
                    DdNode* temp = Cudd_bddAnd(_ddManager, _allBDD[ithCircuit][i][j], Cudd_Not(func));
                    Cudd_Ref(temp);
                    Cudd_RecursiveDeref(_ddManager, _allBDD[ithCircuit][i][j]);
                    _allBDD[ithCircuit][i][j] = temp;
                } 
            }
        }
        else
        {
            if (power < 0)
            {
                power = -power;
                _k[ithCircuit] += 2 * power;
                DdNode* temp = Cudd_Not(func);
                Cudd_Ref(temp);
                Cudd_RecursiveDeref(_ddManager, func);
                func = temp;
            }
            
            // extend BDDs
            _r += power;
            for (int indexOfCircuit = 0; indexOfCircuit < _nCircuit; indexOfCircuit++)   // each BDD of each circuit should be extend together
                allocBDD(_allBDD[indexOfCircuit], power, true); // add new BDDs
            
            // shift
            for(int i = 0; i < _w; ++i)
            {
                for(int j = _r - 1; j >= 0; --j)
                {
                    if (j >= power)
                    {
                        DdNode* shifted = Cudd_bddAnd(_ddManager, _allBDD[ithCircuit][i][j - power], func);
                        Cudd_Ref(shifted);
                        DdNode* unaffected = Cudd_bddAnd(_ddManager, _allBDD[ithCircuit][i][j], Cudd_Not(func));
                        Cudd_Ref(unaffected);
                        
                        DdNode* temp = Cudd_bddOr(_ddManager, shifted, unaffected);
                        Cudd_Ref(temp);
                        Cudd_RecursiveDeref(_ddManager, shifted);
                        Cudd_RecursiveDeref(_ddManager, unaffected);
                        Cudd_RecursiveDeref(_ddManager, _allBDD[ithCircuit][i][j]);
                        _allBDD[ithCircuit][i][j] = temp;
                    }
                    else 
                    {
                        DdNode* unaffected = Cudd_bddAnd(_ddManager, _allBDD[ithCircuit][i][j], Cudd_Not(func));
                        Cudd_Ref(unaffected);
                        Cudd_RecursiveDeref(_ddManager, _allBDD[ithCircuit][i][j]);
                        _allBDD[ithCircuit][i][j] = unaffected;
                    }
                    
                } 
            }
            
            // negate    // please refer to Pauli Z
            if (sign == -1)
            {
                int overflow_done = 0;
                for(int i = 0; i < _w; ++i)
                {
                    DdNode* carry = func;
                    Cudd_Ref(carry);
                    for(int j = 0; j < _r; ++j)
                    {
                        DdNode* temp1 = Cudd_bddAnd(_ddManager, _allBDD[ithCircuit][i][j], Cudd_Not(func));
                        Cudd_Ref(temp1);
                        DdNode* temp2 = Cudd_bddAnd(_ddManager, Cudd_Not(_allBDD[ithCircuit][i][j]), func);
                        Cudd_Ref(temp2);
                        DdNode* inter = Cudd_bddOr(_ddManager, temp1, temp2);
                        Cudd_Ref(inter);
                        Cudd_RecursiveDeref(_ddManager, temp1);
                        Cudd_RecursiveDeref(_ddManager, temp2);
                        Cudd_RecursiveDeref(_ddManager, _allBDD[ithCircuit][i][j]);
                        
                  			// Detect overflow
                  			if ((j == _r - 1) && !overflow_done)
                        {
                    				if (overflow2(inter, carry))
                    				{
                      					_r += _inc;
                      					for (int indexOfCircuit = 0; indexOfCircuit < _nCircuit; indexOfCircuit++)   // each BDD of each circuit should be extend together
                      						  allocBDD(_allBDD[indexOfCircuit], _inc, true); // add new BDDs
                      					overflow_done = 1;
                    				}
                        }
              
                        // Plus 1
                        if (carry == _zeroNode)
                            _allBDD[ithCircuit][i][j] = inter;
                        else
                        {
                            // Sum
                            _allBDD[ithCircuit][i][j] = Cudd_bddXor(_ddManager, inter, carry);
                            Cudd_Ref(_allBDD[ithCircuit][i][j]);
                            
                            // Carry
                            if (j == _r - 1)
                    				{
                    					  Cudd_RecursiveDeref(_ddManager, carry);
                                Cudd_RecursiveDeref(_ddManager, inter);
                    				}
                            else
                            {
                                DdNode* tmp = Cudd_bddAnd(_ddManager, inter, carry);
                                Cudd_Ref(tmp);
                                Cudd_RecursiveDeref(_ddManager, carry);
                                Cudd_RecursiveDeref(_ddManager, inter);
                                carry = tmp;
                            }
                        }

                    }
                }
            } // end of "if (sign == -1)"
        } // end of "if (sign == 0) else"  
        Cudd_RecursiveDeref(_ddManager, func);
    }
    
    // sum // please refer to Hadamard
    for (int ith_var = _nQm; ith_var < _nQm + _nQw; ++ith_var) {
        DdNode *column_variable = Cudd_bddIthVar(_ddManager, ith_var + _n);
        Cudd_Ref(column_variable);
        
        int overflow_done = 0;
        DdNode *g, *d, *c, *tmp, *term1, *term2;
        
        for (int i = 0; i < _w; i++) 
        {
            c = _zeroNode; 
            Cudd_Ref(c);
            for (int j = 0; j < _r; j++)
            {
                // G
                g = Cudd_Cofactor(_ddManager, _allBDD[ithCircuit][i][j], column_variable);
                Cudd_Ref(g);
    
                // D
                d = Cudd_Cofactor(_ddManager, _allBDD[ithCircuit][i][j], Cudd_Not(column_variable));
                Cudd_Ref(d);
    
                // Detect overflow
                if ((j == _r - 1) && !overflow_done)
                    if (overflow3(g, d, c))
                    {
                        _r += _inc;
                        for (int indexOfCircuit = 0; indexOfCircuit < _nCircuit; indexOfCircuit++)   // each BDD of each circuit should be extend together
                            allocBDD(_allBDD[indexOfCircuit], _inc, true); // add new BDDs
                        overflow_done = 1;
                    }
    
                // Sum
                Cudd_RecursiveDeref(_ddManager, _allBDD[ithCircuit][i][j]);
                _allBDD[ithCircuit][i][j] = Cudd_bddXor(_ddManager, g, d);
                Cudd_Ref(_allBDD[ithCircuit][i][j]);
                tmp = Cudd_bddXor(_ddManager, _allBDD[ithCircuit][i][j], c);
                Cudd_Ref(tmp);
                Cudd_RecursiveDeref(_ddManager, _allBDD[ithCircuit][i][j]);
                _allBDD[ithCircuit][i][j] = tmp;
    
                // Carry
                if (j == _r - 1)
                {
                    Cudd_RecursiveDeref(_ddManager, c);
                    Cudd_RecursiveDeref(_ddManager, g);
                    Cudd_RecursiveDeref(_ddManager, d);
                }
                else
                {
                    term1 = Cudd_bddAnd(_ddManager, g, d);
                    Cudd_Ref(term1);
                    term2 = Cudd_bddOr(_ddManager, g, d);
                    Cudd_Ref(term2);
                    Cudd_RecursiveDeref(_ddManager, g);
                    Cudd_RecursiveDeref(_ddManager, d);
                    tmp = Cudd_bddAnd(_ddManager, term2, c);
                    Cudd_Ref(tmp);
                    Cudd_RecursiveDeref(_ddManager, term2);
                    Cudd_RecursiveDeref(_ddManager, c);
                    term2 = tmp;
                    c = Cudd_bddOr(_ddManager, term1, term2);
                    Cudd_Ref(c);
                    Cudd_RecursiveDeref(_ddManager, term1);
                    Cudd_RecursiveDeref(_ddManager, term2);
                }
            }
        }
        Cudd_RecursiveDeref(_ddManager, column_variable);
    }
           
    // Compose
    for (int ith_var = 0; ith_var < _nQm; ++ith_var) {
        int column_variable = ith_var + _n;
        DdNode *copy_variable = Cudd_bddIthVar(_ddManager, ith_var + _n + _n + _nQextra);
        Cudd_Ref(copy_variable);
        
        for (int i = 0; i < _w; i++) 
        {
            for (int j = 0; j < _r; j++)
            {
                DdNode* temp = Cudd_bddCompose(_ddManager, _allBDD[ithCircuit][i][j], copy_variable, column_variable);
                Cudd_Ref(temp);
                Cudd_RecursiveDeref(_ddManager, _allBDD[ithCircuit][i][j]);
                _allBDD[ithCircuit][i][j] = temp;
            }
        }
        
        Cudd_RecursiveDeref(_ddManager, copy_variable);
    }
    
    // F <- F AND [X] = [Y]
    for(int i = 0; i < _w; ++i)
    {
        for(int j = 0; j < _r; ++j)
        {
            for(int ith_var = 0; ith_var < _n; ++ith_var)
            {        
                DdNode *row_variable = Cudd_bddIthVar(_ddManager, ith_var);
                Cudd_Ref(row_variable);
                DdNode *column_variable = Cudd_bddIthVar(_ddManager, ith_var + _n);
                Cudd_Ref(column_variable);
                
                DdNode *temp1 = Cudd_bddXnor(_ddManager, row_variable, column_variable);
                Cudd_Ref(temp1);
                Cudd_RecursiveDeref(_ddManager, row_variable);
                Cudd_RecursiveDeref(_ddManager, column_variable);
                
                DdNode *temp2 = Cudd_bddAnd(_ddManager, _allBDD[ithCircuit][i][j], temp1);
                Cudd_Ref(temp2);
                Cudd_RecursiveDeref(_ddManager, temp1);
                Cudd_RecursiveDeref(_ddManager, _allBDD[ithCircuit][i][j]);
                
                _allBDD[ithCircuit][i][j] = temp2;
            }
        }
    }
    
    if (_isReorder) Cudd_AutodynDisable(_ddManager);
}

/**Function*************************************************************

  Synopsis    [Extract the information needed to be compared between two circuits.]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/

void EquivalenceChecker::extract(int ithCircuit){
    if (_isReorder) Cudd_AutodynEnable(_ddManager, CUDD_REORDER_SYMM_SIFT);
    
    if (_eqType == EqType::CWPeq)
    {
        prepare_init(ithCircuit);
        
        // Apply Ci
        for(int cntCir = 0; cntCir <  _gates[ithCircuit].size(); ++cntCir)
            applyGate(ithCircuit, _gates[ithCircuit][cntCir], _qubits[ithCircuit][cntCir], 0);            
    }
    
    DdNode* mask;
    if (_careSet.empty())
    {
        mask = Cudd_Not(_zeroNode);
        Cudd_Ref(mask);
    }
    else
    {
        mask = func2node(_careSet, _n + _nQd); // already referenced
    }
    
    for(int i = 0; i < _w; ++i)
    {
        for(int j = 0; j < _r; ++j)
        {
            // F <- F|yi
            for(int ith_var = _nQd; ith_var < _nQd + _nQkc; ++ith_var)
            {        
                DdNode *column_variable = Cudd_bddIthVar(_ddManager, ith_var + _n);
                Cudd_Ref(column_variable);
                
                DdNode *temp = Cudd_Cofactor(_ddManager, _allBDD[ithCircuit][i][j], Cudd_Not(column_variable));
                Cudd_Ref(temp);
                Cudd_RecursiveDeref(_ddManager, _allBDD[ithCircuit][i][j]);
                Cudd_RecursiveDeref(_ddManager, column_variable);
                
                _allBDD[ithCircuit][i][j] = temp;
            }
            
            // Compose
            for(int ith_var = 0; ith_var < _nQp; ++ith_var)
            {        
                int row_variable = _nQm + _nQw + ith_var;
                DdNode *column_variable = Cudd_bddIthVar(_ddManager, _nQd + _nQm + _nQw + _nQp + ith_var + _n);
                Cudd_Ref(column_variable);
                
                DdNode *temp = Cudd_bddCompose(_ddManager, _allBDD[ithCircuit][i][j], column_variable, row_variable);
                Cudd_Ref(temp);
                Cudd_RecursiveDeref(_ddManager, _allBDD[ithCircuit][i][j]);
                Cudd_RecursiveDeref(_ddManager, column_variable);
                _allBDD[ithCircuit][i][j] = temp;
            }
            
            // F <- F AND [X]=[Y]
            for(int ith_var = 0; ith_var < _nQm + _nQw + _nQp; ++ith_var)
            {        
                DdNode *row_variable = Cudd_bddIthVar(_ddManager, ith_var);
                Cudd_Ref(row_variable);
                DdNode *column_variable = Cudd_bddIthVar(_ddManager, ith_var + _nQd + _n);
                Cudd_Ref(column_variable);
                
                DdNode *temp1 = Cudd_bddXnor(_ddManager, row_variable, column_variable);
                Cudd_Ref(temp1);
                Cudd_RecursiveDeref(_ddManager, row_variable);
                Cudd_RecursiveDeref(_ddManager, column_variable);
                
                DdNode *temp2 = Cudd_bddAnd(_ddManager, _allBDD[ithCircuit][i][j], temp1);
                Cudd_Ref(temp2);
                Cudd_RecursiveDeref(_ddManager, temp1);
                Cudd_RecursiveDeref(_ddManager, _allBDD[ithCircuit][i][j]);
                
                _allBDD[ithCircuit][i][j] = temp2;
            }
                                    
            // F <- F AND [Y] \in S
            DdNode *temp_mask = Cudd_bddAnd(_ddManager, _allBDD[ithCircuit][i][j], mask);
            Cudd_Ref(temp_mask);
            Cudd_RecursiveDeref(_ddManager, _allBDD[ithCircuit][i][j]);
            _allBDD[ithCircuit][i][j] = temp_mask; 
        }
    }      
    
    // Apply C^-1
    invertCircuit(_gates[ithCircuit]);
    for(int cntCir = _gates[ithCircuit].size() - 1; cntCir >= 0 ; cntCir --)
        applyGate(ithCircuit, _gates[ithCircuit][cntCir], _qubits[ithCircuit][cntCir], 0);
    invertCircuit(_gates[ithCircuit]);
        
    // F <- F AND NOT xi
    for(int i = 0; i < _w; ++i)
    {
        for(int j = 0; j < _r; ++j)
        {    
            for(int ith_var = _nQd; ith_var < _nQd + _nQkc; ++ith_var)
            {        
                DdNode *row_variable = Cudd_bddIthVar(_ddManager, ith_var);
                Cudd_Ref(row_variable);
                
                DdNode *temp = Cudd_bddAnd(_ddManager, _allBDD[ithCircuit][i][j], Cudd_Not(row_variable));
                Cudd_Ref(temp);
                Cudd_RecursiveDeref(_ddManager, row_variable);
                Cudd_RecursiveDeref(_ddManager, _allBDD[ithCircuit][i][j]);
                
                _allBDD[ithCircuit][i][j] = temp;
            }
        }
    }
    Cudd_RecursiveDeref(_ddManager, mask);

    // weighted 
    for (std::pair<std::vector<std::string>, std::pair<int, int>>& item : _weightFuns[ithCircuit])
    {
        DdNode *func = func2node(item.first, _nQd + _n);   // already referenced
        int sign = item.second.first;    // 0 for 0, 1 for +, -1 for -
        int power = item.second.second;  // i for 2^is
        
        if (sign == 0)
        {  
            for(int i = 0; i < _w; ++i)
            {
                for(int j = 0; j < _r; ++j)
                {
                    DdNode* temp = Cudd_bddAnd(_ddManager, _allBDD[ithCircuit][i][j], Cudd_Not(func));
                    Cudd_Ref(temp);
                    Cudd_RecursiveDeref(_ddManager, _allBDD[ithCircuit][i][j]);
                    _allBDD[ithCircuit][i][j] = temp;
                } 
            }
        }
        else
        {
            if (power < 0)
            {
                power = -power;
                _k[ithCircuit] += 2 * power;
                DdNode* temp = Cudd_Not(func);
                Cudd_Ref(temp);
                Cudd_RecursiveDeref(_ddManager, func);
                func = temp;
            }
            
            // extend BDDs
            _r += power;
            for (int indexOfCircuit = 0; indexOfCircuit < _nCircuit; indexOfCircuit++)   // each BDD of each circuit should be extend together
                allocBDD(_allBDD[indexOfCircuit], power, true); // add new BDDs
            
            // shift
            for(int i = 0; i < _w; ++i)
            {
                for(int j = _r - 1; j >= 0; --j)
                {
                    if (j >= power)
                    {
                        DdNode* shifted = Cudd_bddAnd(_ddManager, _allBDD[ithCircuit][i][j - power], func);
                        Cudd_Ref(shifted);
                        DdNode* unaffected = Cudd_bddAnd(_ddManager, _allBDD[ithCircuit][i][j], Cudd_Not(func));
                        Cudd_Ref(unaffected);
                        
                        DdNode* temp = Cudd_bddOr(_ddManager, shifted, unaffected);
                        Cudd_Ref(temp);
                        Cudd_RecursiveDeref(_ddManager, shifted);
                        Cudd_RecursiveDeref(_ddManager, unaffected);
                        Cudd_RecursiveDeref(_ddManager, _allBDD[ithCircuit][i][j]);
                        _allBDD[ithCircuit][i][j] = temp;
                    }
                    else 
                    {
                        DdNode* unaffected = Cudd_bddAnd(_ddManager, _allBDD[ithCircuit][i][j], Cudd_Not(func));
                        Cudd_Ref(unaffected);
                        Cudd_RecursiveDeref(_ddManager, _allBDD[ithCircuit][i][j]);
                        _allBDD[ithCircuit][i][j] = unaffected;
                    }
                    
                } 
            }
            
            // negate    // please refer to Pauli Z
            if (sign == -1)
            {
                int overflow_done = 0;
                for(int i = 0; i < _w; ++i)
                {
                    DdNode* carry = func;
                    Cudd_Ref(carry);
                    for(int j = 0; j < _r; ++j)
                    {
                        DdNode* temp1 = Cudd_bddAnd(_ddManager, _allBDD[ithCircuit][i][j], Cudd_Not(func));
                        Cudd_Ref(temp1);
                        DdNode* temp2 = Cudd_bddAnd(_ddManager, Cudd_Not(_allBDD[ithCircuit][i][j]), func);
                        Cudd_Ref(temp2);
                        DdNode* inter = Cudd_bddOr(_ddManager, temp1, temp2);
                        Cudd_Ref(inter);
                        Cudd_RecursiveDeref(_ddManager, temp1);
                        Cudd_RecursiveDeref(_ddManager, temp2);
                        Cudd_RecursiveDeref(_ddManager, _allBDD[ithCircuit][i][j]);
                        
                  			// Detect overflow
                  			if ((j == _r - 1) && !overflow_done)
                        {
                    				if (overflow2(inter, carry))
                    				{
                      					_r += _inc;
                      					for (int indexOfCircuit = 0; indexOfCircuit < _nCircuit; indexOfCircuit++)   // each BDD of each circuit should be extend together
                      						  allocBDD(_allBDD[indexOfCircuit], _inc, true); // add new BDDs
                      					overflow_done = 1;
                    				}
                        }
              
                        // Plus 1
                        if (carry == _zeroNode)
                            _allBDD[ithCircuit][i][j] = inter;
                        else
                        {
                            // Sum
                            _allBDD[ithCircuit][i][j] = Cudd_bddXor(_ddManager, inter, carry);
                            Cudd_Ref(_allBDD[ithCircuit][i][j]);
                            
                            // Carry
                            if (j == _r - 1)
                    				{
                    					  Cudd_RecursiveDeref(_ddManager, carry);
                                Cudd_RecursiveDeref(_ddManager, inter);
                    				}
                            else
                            {
                                DdNode* tmp = Cudd_bddAnd(_ddManager, inter, carry);
                                Cudd_Ref(tmp);
                                Cudd_RecursiveDeref(_ddManager, carry);
                                Cudd_RecursiveDeref(_ddManager, inter);
                                carry = tmp;
                            }
                        }

                    }
                }
            } // end of "if (sign == -1)"
        } // end of "if (sign == 0) else"  
        Cudd_RecursiveDeref(_ddManager, func);
    }
    
    // sum // please refer to Hadamard
    for (int ith_var = _nQm; ith_var < _nQm + _nQw; ++ith_var) {
        DdNode *column_variable = Cudd_bddIthVar(_ddManager, ith_var + _nQd + _n);
        Cudd_Ref(column_variable);
        
        int overflow_done = 0;
        DdNode *g, *d, *c, *tmp, *term1, *term2;
        
        for (int i = 0; i < _w; i++) 
        {
            c = _zeroNode; 
            Cudd_Ref(c);
            for (int j = 0; j < _r; j++)
            {
                // G
                g = Cudd_Cofactor(_ddManager, _allBDD[ithCircuit][i][j], column_variable);
                Cudd_Ref(g);
    
                // D
                d = Cudd_Cofactor(_ddManager, _allBDD[ithCircuit][i][j], Cudd_Not(column_variable));
                Cudd_Ref(d);
    
                // Detect overflow
                if ((j == _r - 1) && !overflow_done)
                    if (overflow3(g, d, c))
                    {
                        _r += _inc;
                        for (int indexOfCircuit = 0; indexOfCircuit < _nCircuit; indexOfCircuit++)   // each BDD of each circuit should be extend together
                            allocBDD(_allBDD[indexOfCircuit], _inc, true); // add new BDDs
                        overflow_done = 1;
                    }
    
                // Sum
                Cudd_RecursiveDeref(_ddManager, _allBDD[ithCircuit][i][j]);
                _allBDD[ithCircuit][i][j] = Cudd_bddXor(_ddManager, g, d);
                Cudd_Ref(_allBDD[ithCircuit][i][j]);
                tmp = Cudd_bddXor(_ddManager, _allBDD[ithCircuit][i][j], c);
                Cudd_Ref(tmp);
                Cudd_RecursiveDeref(_ddManager, _allBDD[ithCircuit][i][j]);
                _allBDD[ithCircuit][i][j] = tmp;
    
                // Carry
                if (j == _r - 1)
                {
                    Cudd_RecursiveDeref(_ddManager, c);
                    Cudd_RecursiveDeref(_ddManager, g);
                    Cudd_RecursiveDeref(_ddManager, d);
                }
                else
                {
                    term1 = Cudd_bddAnd(_ddManager, g, d);
                    Cudd_Ref(term1);
                    term2 = Cudd_bddOr(_ddManager, g, d);
                    Cudd_Ref(term2);
                    Cudd_RecursiveDeref(_ddManager, g);
                    Cudd_RecursiveDeref(_ddManager, d);
                    tmp = Cudd_bddAnd(_ddManager, term2, c);
                    Cudd_Ref(tmp);
                    Cudd_RecursiveDeref(_ddManager, term2);
                    Cudd_RecursiveDeref(_ddManager, c);
                    term2 = tmp;
                    c = Cudd_bddOr(_ddManager, term1, term2);
                    Cudd_Ref(c);
                    Cudd_RecursiveDeref(_ddManager, term1);
                    Cudd_RecursiveDeref(_ddManager, term2);
                }
            }
        }
        Cudd_RecursiveDeref(_ddManager, column_variable);
    }
            
    if (_eqType == EqType::CWPeq) // need to remove some entries
    {
        DdNode *mask = Cudd_Not(_zeroNode);
        Cudd_Ref(mask);
        for (int ith_var = 0; ith_var < _nQm; ++ith_var)
        {
            DdNode *column_variable = Cudd_bddIthVar(_ddManager, ith_var + _nQd + _n);
            Cudd_Ref(column_variable);
            DdNode *copy_variable = Cudd_bddIthVar(_ddManager, ith_var + _n + _n + _nQextra);
            Cudd_Ref(copy_variable);
            
            DdNode *temp1 = Cudd_bddXnor(_ddManager, column_variable, copy_variable);
            Cudd_Ref(temp1);
            Cudd_RecursiveDeref(_ddManager, column_variable);
            Cudd_RecursiveDeref(_ddManager, copy_variable);
            
            DdNode *temp2 = Cudd_bddAnd(_ddManager, mask, temp1);
            Cudd_Ref(temp2);
            Cudd_RecursiveDeref(_ddManager, mask);
            Cudd_RecursiveDeref(_ddManager, temp1);
            
            mask = temp2;
        }
    
        for(int i = 0; i < _w; i++)
        {
            for(int j = 0; j < _r; j++)
            {
                DdNode *temp = Cudd_bddAnd(_ddManager, _allBDD[ithCircuit][i][j], mask);
                Cudd_Ref(temp);
                Cudd_RecursiveDeref(_ddManager, _allBDD[ithCircuit][i][j]);
                _allBDD[ithCircuit][i][j] = temp;
            }
        }
        
        Cudd_RecursiveDeref(_ddManager, mask);
    }
    
    if (_isReorder) Cudd_AutodynDisable(_ddManager);
}

/**Function*************************************************************

  Synopsis    [Check if two circuits are fully equivalent.]

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

  Synopsis    [Check if two circuits are partially equivalent.]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/

void EquivalenceChecker::checkPeq()
{
    assert((_k[0]-_k[1]) % 2 == 0); // _k[0] - _k[1] should be even because each gate (including H) will be applied pairwisely in (U^-1)U

    int small, large, dk;
    if (_k[0] >= _k[1])
    {
        small = 1;
        large = 0;
        dk = (_k[0] - _k[1]) / 2;
    }
    else
    {
        small = 0;
        large = 1;
        dk = (_k[1] - _k[0]) / 2;
    }   // need to consider different k

    for(int i = 0; i < _w; i++)
    {
        for (int j = 0; j < dk; j++)
        {
            // for 1's complement, the higher bits should be filled with the same bit
            if (_allBDD[small][i][_r - dk + j] != _allBDD[small][i][_r - dk - 1])
            {
                _isEq = false;
                return;
            }
            if (_allBDD[large][i][j] != _zeroNode)
            {
                _isEq = false;
                return;
            }
        }
        for (int j = 0; j < _r - dk; j++)
        {
            if (_allBDD[small][i][j] != _allBDD[large][i][j + dk])
            {
                _isEq = false;
                return;
            }
        }
    }

    _isEq = true;
    return;
}



/**Function*************************************************************

  Synopsis    [Build a BdNode* type from a function.   ]

  Description []

  SideEffects [The returned node is already referenced.]

  SeeAlso     []

***********************************************************************/


DdNode* EquivalenceChecker::func2node(std::vector<std::string>& func, int bias)
{
    assert(!func.empty());
    
    std::stack<DdNode*> waiting;
    for (int i = 0; i < func.size(); ++i)
    {
        if (func[i] == "not")
        {
            assert(waiting.size() >= 1);
            DdNode* temp = waiting.top();
            waiting.pop();
            
            DdNode* new_item = Cudd_Not(temp);
            Cudd_Ref(new_item);
            waiting.push(new_item);
            
            Cudd_RecursiveDeref(_ddManager, temp);
        }
        else if (func[i] == "and" || func[i] == "or" || func[i] == "xor")
        {
            assert(waiting.size() >= 2);
            DdNode* temp_1 = waiting.top();
            waiting.pop();
            DdNode* temp_2 = waiting.top();
            waiting.pop();
             
            DdNode* new_item;              
                 if (func[i] == "and") new_item = Cudd_bddAnd(_ddManager, temp_1, temp_2);
            else if (func[i] == "or")  new_item = Cudd_bddOr (_ddManager, temp_1, temp_2);
            else if (func[i] == "xor") new_item = Cudd_bddXor(_ddManager, temp_1, temp_2);
            Cudd_Ref(new_item);
            waiting.push(new_item);
            
            Cudd_RecursiveDeref(_ddManager, temp_1);
            Cudd_RecursiveDeref(_ddManager, temp_2);
        }
        else
        {
            int ith_var = stoi(func[i]);
            waiting.push(Cudd_bddIthVar(_ddManager, ith_var + bias));
            Cudd_Ref(waiting.top());
        }
    }
    
    assert(waiting.size() == 1);
    DdNode *result = waiting.top();
    waiting.pop();
    return result;
}



/**Function*************************************************************

  Synopsis    [Check if two circuits are partially equivalent in the special case.]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/

void EquivalenceChecker::checkPeqS()
{  
    // M <- [X] \notin S AND [Y] \notin S
    DdNode *mask;
    DdNode *temp_row, *temp_col;
    if (_careSet.empty())    // universal set
    {
        mask = _zeroNode;
        Cudd_Ref(mask);
        temp_row = _zeroNode;
        Cudd_Ref(temp_row);
        temp_col = _zeroNode;
        Cudd_Ref(temp_col);
    }
    else 
    {
        temp_row = func2node(_careSet, 0);    // already referenced
        temp_col = func2node(_careSet, _n);   // already referenced
                
        mask = Cudd_bddAnd(_ddManager, Cudd_Not(temp_row), Cudd_Not(temp_col)); 
        Cudd_Ref(mask);
        Cudd_RecursiveDeref(_ddManager, temp_row);
        Cudd_RecursiveDeref(_ddManager, temp_col);
    }
    
    // M <- M OR [X]=[Y]
    DdNode *temp_xeqy = Cudd_Not(_zeroNode);
    Cudd_Ref(temp_xeqy);
    for (int ith_var = 0; ith_var < _nQm; ++ith_var)
    {
        DdNode *row_variable = Cudd_bddIthVar(_ddManager, ith_var);
        Cudd_Ref(row_variable);
        DdNode *column_variable = Cudd_bddIthVar(_ddManager, ith_var + _n);
        Cudd_Ref(column_variable); 
        
        DdNode *temp1, *temp2;
        temp1 = Cudd_bddXnor(_ddManager, row_variable, column_variable);
        Cudd_Ref(temp1);
        temp2 = Cudd_bddAnd(_ddManager, temp_xeqy, temp1);
        Cudd_Ref(temp2);
        Cudd_RecursiveDeref(_ddManager, temp1);
        Cudd_RecursiveDeref(_ddManager, temp_xeqy);
        temp_xeqy = temp2;
        
        Cudd_RecursiveDeref(_ddManager, row_variable);
        Cudd_RecursiveDeref(_ddManager, column_variable);        
    }
    DdNode *temp_newmask = Cudd_bddOr(_ddManager, mask, temp_xeqy);
    Cudd_RecursiveDeref(_ddManager, mask);
    Cudd_RecursiveDeref(_ddManager, temp_xeqy);        
    mask = temp_newmask;
    
    
    
    // F AND M ?= F
    for(int i = 0; i < _w; i++)
    {
        for(int j = 0; j < _r; j++)
        {
            DdNode *temp;
            temp = Cudd_bddAnd(_ddManager, _allBDD[0][i][j], mask);
            Cudd_Ref(temp);
            
            if (temp != _allBDD[0][i][j])
            {
                Cudd_RecursiveDeref(_ddManager, mask);
                Cudd_RecursiveDeref(_ddManager, temp);
                _isEq = false;                
                return;
            }
            Cudd_RecursiveDeref(_ddManager, temp);
        }
    }
    
    // F <- F AND [X]=[Y]       // F <- F AND [X] \in S        // F00 ?= F11
    for(int i = 0; i < _w; i++)
    {
        for(int j = 0; j < _r; j++)
        {
            // main  
            // F <- F AND [X]=[Y]
            for (int ith_var = 0; ith_var < _nQm; ++ith_var)
            {
                DdNode *row_variable = Cudd_bddIthVar(_ddManager, ith_var);
                Cudd_Ref(row_variable);
                DdNode *column_variable = Cudd_bddIthVar(_ddManager, ith_var + _n);
                Cudd_Ref(column_variable); 
                
                DdNode *temp1, *temp2;
                temp1 = Cudd_bddXnor(_ddManager, row_variable, column_variable);
                Cudd_Ref(temp1);
                temp2 = Cudd_bddAnd(_ddManager, _allBDD[0][i][j], temp1);
                Cudd_Ref(temp2);
                Cudd_RecursiveDeref(_ddManager, temp1);
                Cudd_RecursiveDeref(_ddManager, _allBDD[0][i][j]);
                _allBDD[0][i][j] = temp2;
                
                Cudd_RecursiveDeref(_ddManager, row_variable);
                Cudd_RecursiveDeref(_ddManager, column_variable);        
            }
            
            // F <- F AND [X] \in S
            DdNode *temp;
            temp = Cudd_bddAnd(_ddManager, _allBDD[0][i][j], temp_row);
            Cudd_Ref(temp);
            Cudd_RecursiveDeref(_ddManager, _allBDD[0][i][j]);
            _allBDD[0][i][j] = temp;
            
            // F00 ?= F11
            for (int ith_var = _nQm; ith_var < _nQm + _nQp; ++ith_var)
            {
                DdNode *row_variable = Cudd_bddIthVar(_ddManager, ith_var);
                Cudd_Ref(row_variable);
                DdNode *column_variable = Cudd_bddIthVar(_ddManager, ith_var + _n);
                Cudd_Ref(column_variable); 
                
                DdNode *temp_00, *temp_11;
                temp_00 = Cudd_Cofactor(_ddManager, Cudd_Cofactor(_ddManager, _allBDD[0][i][j], Cudd_Not(row_variable)), Cudd_Not(column_variable));
                Cudd_Ref(temp_00);
                temp_11 = Cudd_Cofactor(_ddManager, Cudd_Cofactor(_ddManager, _allBDD[0][i][j], row_variable), column_variable);
                Cudd_Ref(temp_11);
                if (temp_00 != temp_11)
                {
                    Cudd_RecursiveDeref(_ddManager, temp_00);
                    Cudd_RecursiveDeref(_ddManager, temp_11); 
                    Cudd_RecursiveDeref(_ddManager, row_variable);
                    Cudd_RecursiveDeref(_ddManager, column_variable);
                    Cudd_RecursiveDeref(_ddManager, mask);
                    _isEq = false;
                    return;
                }
                Cudd_RecursiveDeref(_ddManager, temp_00);
                Cudd_RecursiveDeref(_ddManager, temp_11);
                
                Cudd_RecursiveDeref(_ddManager, row_variable);
                Cudd_RecursiveDeref(_ddManager, column_variable);        
            }
        }
    }
    
    Cudd_RecursiveDeref(_ddManager, mask);
    _isEq = true;
    return;
}

/**Function*************************************************************

  Synopsis    [Print the equivalence checking result.]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void EquivalenceChecker::printResult() const
{
    std::cout << "{\n"
              << "\t#Qubits (n): " << _n << std::endl;
    
    if (_eqType == EqType::Feq)
    {
        std::cout << "\tGatecount of circuit1: " << ((_isGatesSwap)? _gates[1].size() : _gates[0].size()) << std::endl;
        std::cout << "\tGatecount of circuit2: " << ((_isGatesSwap)? _gates[0].size() : _gates[1].size()) << std::endl;
        std::cout << "\tIs equivalent? " << (_isEq ? "Yes" : "No") << std::endl;
    }
    else
    {
        std::cout << "\t#Data qubits (n_d):                    "     << _nQd  << std::endl;
        std::cout << "\t#Clean ancilla qubits (n_kc):          "     << _nQkc << std::endl;
        std::cout << "\t#Dirty ancilla qubits (n_kd):          "     << _nQkd << std::endl;
        std::cout << "\t#Measured qubits (n_m):                "     << _nQm  << std::endl;
        std::cout << "\t#Weighted qubits (n_w):                "     << _nQw  << std::endl;
        std::cout << "\t#Propagating qubits (n_p):             "     << _nQp  << std::endl;
        std::cout << "\t#Garbage qubits (n_g):                 "     << _nQg  << std::endl;
        std::cout << "\t#Reverted clean ancilla qubits (n_kr): "     << _nQkr << std::endl;
        std::cout << "\tGatecount of circuit1: " << ((_isGatesSwap)? _gates[1].size() : _gates[0].size()) << std::endl;
        std::cout << "\tGatecount of circuit2: " << ((_isGatesSwap)? _gates[0].size() : _gates[1].size()) << std::endl;
        
             if (_eqType == EqType::SPeq || _eqType == EqType::SPeqS)  std::cout << "\tIs strongly partially equivalent? ";
        else if (_eqType == EqType::CWPeq)                             std::cout << "\tIs weakly partially equivalent? ";
        else                                                           std::cerr << "Unknown equivalence type\n";

             if (_isEq == 1)  std::cout << "Yes\n";
        else if (_isEq == 0)  std::cout << "No\n";
        else if (_isEq == -1) std::cout << "No (unreverted ancilla qubits)\n";
        else                  std::cerr << "Unknown equivalence result\n";
    }
    std::cout << "}" << std::endl;
}

/**Function*************************************************************

  Synopsis    [Print statistics.]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void EquivalenceChecker::printInfo(double runtime, size_t memPeak) const
{
    std::cout << '\n';
    std::cout << "Runtime: " << runtime << " seconds\n";
    std::cout << "Peak memory usage: " << memPeak << " bytes\n"; 
}
