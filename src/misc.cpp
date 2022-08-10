#include "Simulator.h"
#include "util_sim.h"


/**Function*************************************************************

  Synopsis    [initialize an identity matrix represented by BDDs]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void BDDbased::init_identity()
{
    k = new int [nCircuit];
    All_BDD = new DdNode ***[nCircuit];
    for (int ithCircuit = 0; ithCircuit < nCircuit; ithCircuit++){
        k[ithCircuit] = 0;

        DdNode *tmp, *tmp1, *tmp2;
        All_BDD[ithCircuit] = new DdNode **[w];
        for (int i = 0; i < w; i++)
            All_BDD[ithCircuit][i] = new DdNode *[r];

        // reorder into a better order for identity
        int * permut = new int[2*n];
        for (int i = 0; i < n; i++)
        {
            permut[2*i] = i;
            permut[2*i+1] = i + n;
        }
        Cudd_ShuffleHeap(manager, permut);
        delete [] permut;

        for (int i = 0; i < r; i++)
        {
            if (i == 0)
            {
                for (int j = 0; j < w - 1; j++)
                {
                    All_BDD[ithCircuit][j][i] = Cudd_Not(Cudd_ReadOne(manager));
                    Cudd_Ref(All_BDD[ithCircuit][j][i]);
                }
                All_BDD[ithCircuit][w - 1][i] = Cudd_ReadOne(manager);
                Cudd_Ref(All_BDD[ithCircuit][w - 1][i]);
                for (int j = n - 1; j >= 0; j--)
                {
                    tmp1 = Cudd_bddAnd(manager, Cudd_Not(Cudd_bddIthVar(manager, j)), Cudd_Not(Cudd_bddIthVar(manager, j+n)));
                    Cudd_Ref(tmp1);
                    tmp2 = Cudd_bddAnd(manager, Cudd_bddIthVar(manager, j), Cudd_bddIthVar(manager, j+n));
                    Cudd_Ref(tmp2);
                    tmp = Cudd_bddOr(manager, tmp1, tmp2);
                    Cudd_Ref(tmp);
                    Cudd_RecursiveDeref(manager, tmp1);
                    Cudd_RecursiveDeref(manager, tmp2);
                    tmp1 = Cudd_bddAnd(manager, All_BDD[ithCircuit][w - 1][i], tmp);
                    Cudd_Ref(tmp1);
                    Cudd_RecursiveDeref(manager, All_BDD[ithCircuit][w - 1][i]);
                    All_BDD[ithCircuit][w - 1][i] = tmp1;
                }
            }
            else
            {
                for (int j = 0; j < w; j++)
                {
                    All_BDD[ithCircuit][j][i] = Cudd_Not(Cudd_ReadOne(manager));
                    Cudd_Ref(All_BDD[ithCircuit][j][i]);
                }
            }
        }
    }
}

/**Function*************************************************************

  Synopsis    [allocate new BDDs for each integer vector]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void BDDbased::alloc_BDD(DdNode ***Bdd, bool extend)
{
    DdNode *tmp;

    DdNode ***W = new DdNode **[w];
    for (int i = 0; i < w; i++)
        W[i] = new DdNode *[r];

    for (int i = 0; i < r - inc; i++)
        for (int j = 0; j < w; j++)
            W[j][i] = Bdd[j][i];

    for (int i = 0; i < w; i++)
        delete[] Bdd[i];

    for (int i = 0; i < w; i++)
        Bdd[i] = W[i];

    if (extend)
    {
        for (int i = r - inc; i < r; i++)
        {
            for (int j = 0; j < w; j++)
            {
                Bdd[j][i] = Cudd_ReadOne(manager);
                Cudd_Ref(Bdd[j][i]);
                tmp = Cudd_bddAnd(manager, Bdd[j][r - inc - 1], Bdd[j][i]);
                Cudd_Ref(tmp);
                Cudd_RecursiveDeref(manager, Bdd[j][i]);
                Bdd[j][i] = tmp;
            }
        }
    }
}

/**Function*************************************************************

  Synopsis    [detect overflow in integer vectors]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
int BDDbased::overflow3(DdNode *g, DdNode *h, DdNode *crin){
    DdNode *tmp, *dd1, *dd2;
    int overflow;

    dd1 = Cudd_bddXor(manager, g, crin);
    Cudd_Ref(dd1);

    dd2 = Cudd_bddXnor(manager, g, h);
    Cudd_Ref(dd2);

    tmp = Cudd_bddAnd(manager, dd1, dd2);
    Cudd_Ref(tmp);
    Cudd_RecursiveDeref(manager, dd1);
    Cudd_RecursiveDeref(manager, dd2);

    if (Cudd_CountPathsToNonZero(tmp))
        overflow = 1;
    else
        overflow = 0;
    Cudd_RecursiveDeref(manager, tmp);

    return overflow;
}

/**Function*************************************************************

  Synopsis    [detect overflow in integer vectors -- for the case that h is 0]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
int BDDbased::overflow2(DdNode *g, DdNode *crin){
    DdNode *tmp;
    int overflow;

    tmp = Cudd_bddAnd(manager, Cudd_Not(g), crin);
    Cudd_Ref(tmp);

    if (Cudd_CountPathsToNonZero(tmp))
        overflow = 1;
    else
        overflow = 0;
    Cudd_RecursiveDeref(manager, tmp);

    return overflow;
}

/**Function*************************************************************

  Synopsis    [reorder BDDs]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void BDDbased::reorder()
{
    int reorder_signal = Cudd_ReduceHeap(manager, CUDD_REORDER_SYMM_SIFT, 0);
    if (!reorder_signal)
        std::cout << "reorder fails" << std::endl;
}

/**Function*************************************************************

  Synopsis    [update max #nodes]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void BDDbased::nodecount()
{
    unsigned long NodeCount_new = Cudd_ReadNodeCount(manager);
    if (NodeCount_new > NodeCount)
         NodeCount = NodeCount_new;
}
