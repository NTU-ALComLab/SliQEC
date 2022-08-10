#include "Simulator.h"
#include "util_sim.h"

PartialEquivalenceChecker::PartialEquivalenceChecker(std::string G, std::string G_p, int inputNumber, int outputNumber, int bitSize, bool reorder): BDDbased(2, bitSize, reorder)
{
    nInput = inputNumber;
    nOutput = outputNumber;

    gates.resize(2);
    qubits.resize(2);
    circuitParser(G, gates[0], qubits[0]);
    circuitParser(G_p, gates[1], qubits[1]);

    // the longer circuit (in gatecount) is stored in gates[1]
    if (gates[0].size() > gates[1].size())
    {
        std::vector<gateType> tmp;
        tmp = gates[0];
        gates[0] = gates[1];
        gates[1] = tmp;
        tmp.clear();
        std::vector<std::vector<int> > tmp1;
        tmp1 = qubits[0];
        qubits[0] = qubits[1];
        qubits[1] = tmp1;
        tmp1.clear();
    }

    // compute ratio
    //ratio = ceil(((double) gates[1].size()) / ((double) gates[0].size()));
    ratio = round(((double) gates[1].size()) / ((double) gates[0].size()));
}

/**Function*************************************************************

  Synopsis    [setup BDDs]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/

void PartialEquivalenceChecker::setupDD(bool special){
    if (!special){
        int nAncilla = n - nInput;
        if (nOutput > nAncilla){
            n += nOutput - nAncilla;
        }
    }

    manager = Cudd_Init(2*n, 2*n, CUDD_UNIQUE_SLOTS, CUDD_CACHE_SLOTS, 0); // 0~(n-1): 0-variables, n~(2n-1): 1-variables
    init_identity();
    initBaseBDD();
}

/**Function*************************************************************

  Synopsis    [parse and store two gateuits in vectors]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/

void PartialEquivalenceChecker::circuitParser(std::string G, std::vector<gateType> &gate, std::vector<std::vector<int> > &qubit)
{
    std::string inStr;
    std::stringstream inFile_ss(G);
    while (getline(inFile_ss, inStr))
    {
        inStr = inStr.substr(0, inStr.find("//"));
        if (inStr.find_first_not_of("\t\n ") != std::string::npos)
        {
            std::stringstream inStr_ss(inStr);
            getline(inStr_ss, inStr, ' ');
            if (inStr == "qreg")
            {
                getline(inStr_ss, inStr, '[');
                getline(inStr_ss, inStr, ']');

                assert(stoi(inStr) >= nInput);
                assert(stoi(inStr) >= nOutput);
                n = std::max(n, stoi(inStr));
            }
            else if (inStr == "creg"){;}
            else if (inStr == "OPENQASM"){;}
            else if (inStr == "include"){;}
            else if (inStr == "measure"){;}
            else
            {
                if (inStr == "x")
                {
                    gate.push_back(X);
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    int var = stoi(inStr);
                    assert(var < n);
                    qubit.push_back(std::vector<int>(1, var));
                }
                else if (inStr == "y")
                {
                    gate.push_back(Y);
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    int var = stoi(inStr);
                    assert(var < n);
                    qubit.push_back(std::vector<int>(1, var));
                }
                else if (inStr == "z")
                {
                    gate.push_back(Z);
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    int var = stoi(inStr);
                    assert(var < n);
                    qubit.push_back(std::vector<int>(1, var));
                }
                else if (inStr == "h")
                {
                    gate.push_back(H);
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    int var = stoi(inStr);
                    assert(var < n);
                    qubit.push_back(std::vector<int>(1, var));
                }
                else if (inStr == "s")
                {
                    gate.push_back(S);
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    int var = stoi(inStr);
                    assert(var < n);
                    qubit.push_back(std::vector<int>(1, var));
                }
                else if (inStr == "sdg")
                {
                    gate.push_back(SDG);
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    int var = stoi(inStr);
                    assert(var < n);
                    qubit.push_back(std::vector<int>(1, var));
                }
                else if (inStr == "t")
                {
                    gate.push_back(T);
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    int var = stoi(inStr);
                    assert(var < n);
                    qubit.push_back(std::vector<int>(1, var));
                }
                else if (inStr == "tdg")
                {
                    gate.push_back(TDG);
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    int var = stoi(inStr);
                    assert(var < n);
                    qubit.push_back(std::vector<int>(1, var));
                }
                else if (inStr == "rx(pi/2)")
                {
                    gate.push_back(RX_PI_2);
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    int var = stoi(inStr);
                    assert(var < n);
                    qubit.push_back(std::vector<int>(1, var));
                }
                else if (inStr == "rx(-pi/2)")
                {
                    gate.push_back(RX_PI_2_DG);
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    int var = stoi(inStr);
                    assert(var < n);
                    qubit.push_back(std::vector<int>(1, var));
                }
                else if (inStr == "ry(pi/2)")
                {
                    gate.push_back(RY_PI_2);
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    int var = stoi(inStr);
                    assert(var < n);
                    qubit.push_back(std::vector<int>(1, var));
                }
                else if (inStr == "ry(-pi/2)")
                {
                    gate.push_back(RY_PI_2_DG);
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    int var = stoi(inStr);
                    assert(var < n);
                    qubit.push_back(std::vector<int>(1, var));
                }
                else if (inStr == "cx")
                {
                    gate.push_back(CX);
                    std::vector<int> iqubit(2);
                    for (int i = 0; i < 2; i++)
                    {
                        getline(inStr_ss, inStr, '[');
                        getline(inStr_ss, inStr, ']');
                        iqubit[i] = stoi(inStr);
                        assert(iqubit[i] < n);
                    }
                    qubit.push_back(iqubit);
                    iqubit.clear();
                }
                else if (inStr == "cz")
                {
                    gate.push_back(CZ);
                    std::vector<int> iqubit(2);
                    for (int i = 0; i < 2; i++)
                    {
                        getline(inStr_ss, inStr, '[');
                        getline(inStr_ss, inStr, ']');
                        iqubit[i] = stoi(inStr);
                        assert(iqubit[i] < n);
                    }
                    qubit.push_back(iqubit);
                    iqubit.clear();
                }
                else if (inStr == "swap")
                {
                    gate.push_back(SWAP);
                    std::vector<int> iqubit(2);
                    for (int i = 0; i < 2; i++)
                    {
                        getline(inStr_ss, inStr, '[');
                        getline(inStr_ss, inStr, ']');
                        iqubit[i] = stoi(inStr);
                        assert(iqubit[i] < n);
                    }
                    qubit.push_back(iqubit);
                    iqubit.clear();
                }
                else if (inStr == "cswap")
                {
                    gate.push_back(CSWAP);
                    std::vector<int> iqubit(3);
                    for (int i = 0; i < 3; i++)
                    {
                        getline(inStr_ss, inStr, '[');
                        getline(inStr_ss, inStr, ']');
                        iqubit[i] = stoi(inStr);
                        assert(iqubit[i] < n);
                    }
                    qubit.push_back(iqubit);
                    iqubit.clear();
                }
                else if (inStr == "ccx" || inStr == "mcx")
                {
                    gate.push_back(CCX);
                    std::vector<int> iqubit(0);
                    getline(inStr_ss, inStr, '[');
                    while(getline(inStr_ss, inStr, ']'))
                    {
                        iqubit.push_back(stoi(inStr));
                        assert(iqubit.back() < n);
                        getline(inStr_ss, inStr, '[');
                    }
                    qubit.push_back(iqubit);
                    iqubit.clear();
                }
                else
                {
                    std::cerr << std::endl
                            // << "[warning]: Gate \'" << inStr << "\' is not supported in this simulator. The gate is ignored ..." << std::endl;
                            << "[warning]: Syntax \'" << inStr << "\' is not supported in this simulator. The line is ignored ..." << std::endl;
                }
            }
        }
    }
}

/**Function*************************************************************

  Synopsis    [invert the gates in the given circuit]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/

void PartialEquivalenceChecker::invertCircuit(std::vector<gateType> &gate)
{
    for (int i = 0; i < gate.size(); i++)
    {
        if (gate[i] == S) gate[i] = SDG;
        else if (gate[i] == SDG) gate[i] = S;
        else if (gate[i] == T) gate[i] = TDG;
        else if (gate[i] == TDG) gate[i] = T;
        else if (gate[i] == RX_PI_2) gate[i] = RX_PI_2_DG;
        else if (gate[i] == RX_PI_2_DG) gate[i] = RX_PI_2;
        else if (gate[i] == RY_PI_2) gate[i] = RY_PI_2_DG;
        else if (gate[i] == RY_PI_2_DG) gate[i] = RY_PI_2;
    }
}

/**Function*************************************************************

  Synopsis    [initilize zero BDD for checking]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void PartialEquivalenceChecker::initBaseBDD()
{
    zero = Cudd_Not(Cudd_ReadOne(manager));
    Cudd_Ref(zero);
}

/**Function*************************************************************

  Synopsis    [Simulate the miter]

  Description [
               Apply gates in G and G' to evolve matrx from identity interleavingly.
               The longer gateuit (gates[1]) is seen as G (applied to the left side of the matrix)
               to save computation overhead
              ]

  SideEffects []

  SeeAlso     []

***********************************************************************/
void PartialEquivalenceChecker::sim_miter()
{
    int cntCir0 = 0, cntCir1 = 0;

    if (isReorder) Cudd_AutodynEnable(manager, CUDD_REORDER_SYMM_SIFT);

    invertCircuit(gates[1]);
    while (cntCir0 < gates[0].size() || cntCir1 < gates[1].size())
    {
        // apply 1 gate from gates[0]
        if (cntCir0 < gates[0].size())
        {
            applyGate(0, gates[0][cntCir0], qubits[0][cntCir0], false);
            cntCir0++;
        }
        // apply ratio gate(s) from gates[1]
        while(  cntCir1 * gates[0].size() < cntCir0 * gates[1].size()  &&  cntCir1 < gates[1].size()   )  
        {
            applyGate(0, gates[1][cntCir1], qubits[1][cntCir1], true);
            cntCir1++;
        }
    }
    invertCircuit(gates[1]);

    if (isReorder) Cudd_AutodynDisable(manager);
}

/**Function*************************************************************

  Synopsis    [Apply a gate]

  Description [Apply a gate to the left side of the matrix if right = 0, the right side otherwise.]

  SideEffects []

  SeeAlso     []

***********************************************************************/
void PartialEquivalenceChecker::applyGate(int ithCircuit, gateType type, std::vector<int> qubit, bool right)
{
    /*
    std::cout<<type<<"========="<<gatecount<<"=============================\n";
    for(int i=0; i<w; i++){
        for(int j=0; j<r;j++){
            std::cout<<i<<" "<<j<<"\n";
            Cudd_PrintMinterm(manager, All_BDD[0][i][j]);
            std::cout<<"--\n";
            Cudd_PrintMinterm(manager, All_BDD[1][i][j]);
        }
    }*/
    if (right) for (int i = 0; i < qubit.size(); i++) qubit[i] += n;

    if (type == X) PauliX(ithCircuit, qubit[0]);
    else if (type == Y) PauliY(ithCircuit, qubit[0], right);
    else if (type == Z) PauliZ(ithCircuit, qubit);
    else if (type == H) Hadamard(ithCircuit, qubit[0]);
    else if (type == S) Phase_shift(ithCircuit, 2, qubit[0]);
    else if (type == SDG) Phase_shift_dagger(ithCircuit, -2, qubit[0]);
    else if (type == T) Phase_shift(ithCircuit, 4, qubit[0]);
    else if (type == TDG) Phase_shift_dagger(ithCircuit, -4, qubit[0]);
    else if (type == RX_PI_2) rx_pi_2(ithCircuit, qubit[0], false);
    else if (type == RX_PI_2_DG) rx_pi_2(ithCircuit, qubit[0], true);
    else if (type == RY_PI_2) ry_pi_2(ithCircuit, qubit[0], right^false);
    else if (type == RY_PI_2_DG) ry_pi_2(ithCircuit, qubit[0], right^true);
    else if (type == CX)
    {
        std::vector<int> ncont(0);
        int targ = qubit[1];
        qubit.pop_back();
        Toffoli(ithCircuit, targ, qubit, ncont);
        ncont.clear();
    }
    else if (type == CZ) PauliZ(ithCircuit, qubit);
    else if (type == SWAP)
    {
        std::vector<int> cont(0);
        Fredkin(ithCircuit, qubit[0], qubit[1], cont);
        cont.clear();
    }
    else if (type == CSWAP)
    {
        int swapA = qubit[1], swapB = qubit[2];
        qubit.pop_back();
        qubit.pop_back();
        Fredkin(ithCircuit, swapA, swapB, qubit);
    }
    else if (type == CCX)
    {
        std::vector<int> ncont(0);
        int targ = qubit.back();
        qubit.pop_back();
        Toffoli(ithCircuit, targ, qubit, ncont);
        ncont.clear();
    }

    if (manager != NULL)
        nodecount();
}

/**Function*************************************************************

  Synopsis    [Extract information]

  Description [Extract the information needed to compare between two circuits.]

  SideEffects []

  SeeAlso     []

***********************************************************************/

void PartialEquivalenceChecker::extract(int ithCircuit){
    if (isReorder) Cudd_AutodynEnable(manager, CUDD_REORDER_SYMM_SIFT);

    for(int cntCir = 0; cntCir < gates[ithCircuit].size(); cntCir ++){
        applyGate(ithCircuit, gates[ithCircuit][cntCir], qubits[ithCircuit][cntCir], 0);
    }
    // setup G

    for(int i = 0; i < w; i++){
        for(int j = 0; j < r; j++){
            for(int variable_1 = nInput + n; variable_1 < 2*n; variable_1++){      // index of 1-variable
                DdNode *temp1, *temp2;

                temp1 = Cudd_Not(Cudd_bddIthVar(manager, variable_1));
                Cudd_Ref(temp1);

                temp2 = Cudd_Cofactor(manager, All_BDD[ithCircuit][i][j], temp1);
                Cudd_Ref(temp2);
                Cudd_RecursiveDeref(manager, temp1);

                Cudd_RecursiveDeref(manager, All_BDD[ithCircuit][i][j]);
                All_BDD[ithCircuit][i][j] = temp2;
            }
        }
    }
    // cofactor

    for(int i = 0; i < w; i++){
        for(int j = 0; j < r; j++){
            for(int index = 0; index < nOutput; index++){
                int variable_0 = index;                                    // index of 0-variable
                int variable_1 = n + (n - nOutput + index);                // index of 1-variable
                DdNode *temp1, *temp2, *temp3;

                temp1 = Cudd_bddXor(manager, Cudd_bddIthVar(manager, variable_0), Cudd_bddIthVar(manager, variable_1));
                Cudd_Ref(temp1);

                temp2 = Cudd_Not(temp1);
                Cudd_Ref(temp2);
                Cudd_RecursiveDeref(manager, temp1);

                temp3 = Cudd_bddAnd(manager, All_BDD[ithCircuit][i][j], temp2);
                Cudd_Ref(temp3);
                Cudd_RecursiveDeref(manager, temp2);

                Cudd_RecursiveDeref(manager, All_BDD[ithCircuit][i][j]);
                All_BDD[ithCircuit][i][j] = temp3;
            }
        }
    }
    //

    for(int i = 0; i < w; i++){
        for(int j = 0; j < r; j++){
            for(int variable_1 = nInput + n; variable_1 < (n - nOutput) + n; variable_1++){     // index of 1-variable
                DdNode *temp1, *temp2;

                temp1 = Cudd_Not(Cudd_bddIthVar(manager, variable_1));
                Cudd_Ref(temp1);

                temp2 = Cudd_bddAnd(manager, All_BDD[ithCircuit][i][j], temp1);
                Cudd_Ref(temp2);
                Cudd_RecursiveDeref(manager, temp1);

                Cudd_RecursiveDeref(manager, All_BDD[ithCircuit][i][j]);
                All_BDD[ithCircuit][i][j] = temp2;
            }
        }
    }
    //

    invertCircuit(gates[ithCircuit]);
    for(int cntCir = gates[ithCircuit].size() - 1; cntCir >= 0 ; cntCir --){
        applyGate(ithCircuit, gates[ithCircuit][cntCir], qubits[ithCircuit][cntCir], 0);
    }
    invertCircuit(gates[ithCircuit]);
    // apply G^-1

    for(int i = 0; i < w; i++){
        for(int j = 0; j < r; j++){
            for(int variable_0 = nInput; variable_0 < n; variable_0++){
                DdNode *temp1, *temp2;

                temp1 = Cudd_Not( Cudd_bddIthVar(manager, variable_0) );
                Cudd_Ref(temp1);

                temp2 = Cudd_bddAnd(manager, All_BDD[ithCircuit][i][j], temp1);
                Cudd_Ref(temp2);
                Cudd_RecursiveDeref(manager, temp1);

                Cudd_RecursiveDeref(manager, All_BDD[ithCircuit][i][j]);
                All_BDD[ithCircuit][i][j] = temp2;
            }
        }
    }
    //

    if (isReorder) Cudd_AutodynDisable(manager);
}

/**Function*************************************************************

  Synopsis    [Check if two circuits are partially equivalent]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/

bool PartialEquivalenceChecker::checkPEC(){
    if ( (k[0]-k[1])%2 != 0){   // k[0] - k[1] must be even for two matrices to be equivalent
        assert(false);
        // this condition should not appear, because each gate will be applied pairwisely in (U^-1)U
    }

    int small, large, dk;
    if (k[0] >= k[1]){
        small = 1;
        large = 0;
        dk = (k[0] - k[1]) / 2;
    }
    else{
        small = 0;
        large = 1;
        dk = (k[1] - k[0]) / 2;
    }                               // need to consider different k

    for(int i = 0; i < w; i++){
        for (int j = 0; j < dk; j++){
            if (All_BDD[small][i][r - dk + j] != All_BDD[small][i][r - dk - 1]){     // for 1's complement, the higher bits should be filled with the same bit
                return false;
            }
            if (All_BDD[large][i][j] != zero){
                return false;
            }
        }
        for (int j = 0; j < r - dk; j++){
            if (All_BDD[small][i][j] != All_BDD[large][i][j + dk]){
                return false;
            }
        }
    }
    return true;
}

/**Function*************************************************************

  Synopsis    [Check if two circuits are partially equivalent -- in the special case]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/

bool PartialEquivalenceChecker::checkPEC_special(){
    DdNode *mask = zero;
    Cudd_Ref(mask);
    for (int index = 0; index < nOutput; index++){
        int variable_0 = index;
        int variable_1 = n + index;

        DdNode *temp1, *temp2;

        temp1 = Cudd_bddXor(manager, Cudd_bddIthVar(manager, variable_0), Cudd_bddIthVar(manager, variable_1));
        Cudd_Ref(temp1);

        temp2 = Cudd_bddOr(manager, mask, temp1);
        Cudd_RecursiveDeref(manager, temp1);
        Cudd_Ref(temp2);

        Cudd_RecursiveDeref(manager, mask);
        mask = temp2;
    }

    for(int i = 0; i < w; i++){
        for(int j = 0; j < r; j++){
            DdNode *temp;
            temp = Cudd_bddAnd(manager, All_BDD[0][i][j], mask);
            Cudd_Ref(temp);

            if (temp != zero){
                Cudd_RecursiveDeref(manager, mask);
                Cudd_RecursiveDeref(manager, temp);
                return false;
            }
            Cudd_RecursiveDeref(manager, temp);
        }
    }
    Cudd_RecursiveDeref(manager, mask);
    return true;
}

/**Function*************************************************************

  Synopsis    [Output final result]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void PartialEquivalenceChecker::getResult(bool special)
{
    std::cout << std::endl;
    std::cout << "  #Qubits: " << n << std::endl;
    std::cout << "  #Input Qubits: " << nInput << std::endl;
    std::cout << "  #Output Qubits: " << nOutput << std::endl;
    std::cout << "  Gatecount of G : " << gates[0].size() << std::endl;
    std::cout << "  Gatecount of G': " << gates[1].size() << std::endl;
    printf("  |G'|/|G|: %.2f\n", ((double) gates[1].size()) / ((double) gates[0].size()));
    if (special) std::cout << "  Is partially equivalent? (special case) ";
    else std::cout << "  Is partially equivalent? ";
    if (isPEC) std::cout << "Yes" << std::endl;
    else std::cout << "No" << std::endl;
    std::cout << std::endl;
}


/**Function*************************************************************

  Synopsis    [Run PEC procedure for the special case]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void PartialEquivalenceChecker::runPEC()
{
    setupDD(false);
    extract(0);
    extract(1);
    isPEC = checkPEC();
    getResult(0);
}

/**Function*************************************************************

  Synopsis    [Run PEC procedure for the special case]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void PartialEquivalenceChecker::runPEC_special()
{
    if (n != nInput)
        std::cout << "all qubits should be input qubits in the special case\n";
    assert(n == nInput);

    setupDD(true);
    sim_miter();
    isPEC = checkPEC_special();
    getResult(1);
}


/**Function*************************************************************

  Synopsis    [print statistics]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void PartialEquivalenceChecker::print_info(double runtime, size_t memPeak)
{
    std::cout << "  Runtime: " << runtime << " seconds" << std::endl;
    std::cout << "  Peak memory usage: " << memPeak << " bytes" << std::endl; //unit in bytes
    std::cout << "  Max #nodes: " << NodeCount << std::endl;
    std::cout << "  Integer bit size: " << r << std::endl;
}

//Cudd_PrintMinterm(manager, node);
