#include "bddSystem.h"

/**Function*************************************************************

  Synopsis    [Initialize BDD manager/_zeroNode/_identityNode.]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void BDDSystem::ddInitialize()
{
    _ddManager = Cudd_Init(2*_n, 2*_n, CUDD_UNIQUE_SLOTS, CUDD_CACHE_SLOTS, 0); // 0~(n-1): 0-variables, n~(2n-1): 1-variables

    _zeroNode = Cudd_Not(Cudd_ReadOne(_ddManager));
    Cudd_Ref(_zeroNode);
}

/**Function*************************************************************

  Synopsis    [Initialize an identity matrix represented by BDDs.]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void BDDSystem::initIdentity()
{
    _k = new int [_nCircuit];
    _allBDD = new DdNode ***[_nCircuit];
    for (int ithCircuit = 0; ithCircuit < _nCircuit; ithCircuit++){
        _k[ithCircuit] = 0;

        DdNode *tmp, *tmp1, *tmp2;
        _allBDD[ithCircuit] = new DdNode **[_w];
        for (int i = 0; i < _w; i++)
            _allBDD[ithCircuit][i] = new DdNode *[_r];

        // reorder into a better order for identity
        int * permut = new int[2*_n];
        for (int i = 0; i < _n; i++)
        {
            permut[2*i] = i;
            permut[2*i+1] = i + _n;
        }
        Cudd_ShuffleHeap(_ddManager, permut);
        delete [] permut;

        for (int i = 0; i < _r; i++)
        {
            if (i == 0)
            {
                for (int j = 0; j < _w - 1; j++)
                {
                    _allBDD[ithCircuit][j][i] = Cudd_Not(Cudd_ReadOne(_ddManager));
                    Cudd_Ref(_allBDD[ithCircuit][j][i]);
                }
                _allBDD[ithCircuit][_w - 1][i] = Cudd_ReadOne(_ddManager);
                Cudd_Ref(_allBDD[ithCircuit][_w - 1][i]);
                for (int j = _n - 1; j >= 0; j--)
                {
                    tmp1 = Cudd_bddAnd(_ddManager, Cudd_Not(Cudd_bddIthVar(_ddManager, j)), Cudd_Not(Cudd_bddIthVar(_ddManager, j+_n)));
                    Cudd_Ref(tmp1);
                    tmp2 = Cudd_bddAnd(_ddManager, Cudd_bddIthVar(_ddManager, j), Cudd_bddIthVar(_ddManager, j+_n));
                    Cudd_Ref(tmp2);
                    tmp = Cudd_bddOr(_ddManager, tmp1, tmp2);
                    Cudd_Ref(tmp);
                    Cudd_RecursiveDeref(_ddManager, tmp1);
                    Cudd_RecursiveDeref(_ddManager, tmp2);
                    tmp1 = Cudd_bddAnd(_ddManager, _allBDD[ithCircuit][_w - 1][i], tmp);
                    Cudd_Ref(tmp1);
                    Cudd_RecursiveDeref(_ddManager, _allBDD[ithCircuit][_w - 1][i]);
                    _allBDD[ithCircuit][_w - 1][i] = tmp1;
                }

                if(ithCircuit == 0)
                {
                    _identityNode = _allBDD[ithCircuit][_w - 1][i];
                    Cudd_Ref(_identityNode);
                }
             }
            else
            {
                for (int j = 0; j < _w; j++)
                {
                    _allBDD[ithCircuit][j][i] = Cudd_Not(Cudd_ReadOne(_ddManager));
                    Cudd_Ref(_allBDD[ithCircuit][j][i]);
                }
            }
        }
    }
}

/**Function*************************************************************

  Synopsis    [Allocate new BDDs for each integer vector.]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void BDDSystem::allocBDD(DdNode ***Bdd, bool extend)
{
    DdNode *tmp;

    DdNode ***W = new DdNode **[_w];
    for (int i = 0; i < _w; i++)
        W[i] = new DdNode *[_r];

    for (int i = 0; i < _r - _inc; i++)
        for (int j = 0; j < _w; j++)
            W[j][i] = Bdd[j][i];

    for (int i = 0; i < _w; i++)
        delete[] Bdd[i];

    for (int i = 0; i < _w; i++)
        Bdd[i] = W[i];

    if (extend)
    {
        for (int i = _r - _inc; i < _r; i++)
        {
            for (int j = 0; j < _w; j++)
            {
                Bdd[j][i] = Cudd_ReadOne(_ddManager);
                Cudd_Ref(Bdd[j][i]);
                tmp = Cudd_bddAnd(_ddManager, Bdd[j][_r - _inc - 1], Bdd[j][i]);
                Cudd_Ref(tmp);
                Cudd_RecursiveDeref(_ddManager, Bdd[j][i]);
                Bdd[j][i] = tmp;
            }
        }
    }
}

/**Function*************************************************************

  Synopsis    [Detect overflow in integer vectors.]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
int BDDSystem::overflow3(DdNode *g, DdNode *h, DdNode *crin) const
{
    DdNode *tmp, *dd1, *dd2;
    int overflow;

    dd1 = Cudd_bddXor(_ddManager, g, crin);
    Cudd_Ref(dd1);

    dd2 = Cudd_bddXnor(_ddManager, g, h);
    Cudd_Ref(dd2);

    tmp = Cudd_bddAnd(_ddManager, dd1, dd2);
    Cudd_Ref(tmp);
    Cudd_RecursiveDeref(_ddManager, dd1);
    Cudd_RecursiveDeref(_ddManager, dd2);

    if (Cudd_CountPathsToNonZero(tmp))
        overflow = 1;
    else
        overflow = 0;
    Cudd_RecursiveDeref(_ddManager, tmp);

    return overflow;
}

/**Function*************************************************************

  Synopsis    [Detect overflow in integer vectors -- for the case that h is 0.]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
int BDDSystem::overflow2(DdNode *g, DdNode *crin) const
{
    DdNode *tmp;
    int overflow;

    tmp = Cudd_bddAnd(_ddManager, Cudd_Not(g), crin);
    Cudd_Ref(tmp);

    if (Cudd_CountPathsToNonZero(tmp))
        overflow = 1;
    else
        overflow = 0;
    Cudd_RecursiveDeref(_ddManager, tmp);

    return overflow;
}

/**Function*************************************************************

  Synopsis    [Update max #nodes.]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void BDDSystem::updateNodeCount()
{
    _nodeCount = std::max(_nodeCount, static_cast<unsigned long>(Cudd_ReadNodeCount(_ddManager)));
}
