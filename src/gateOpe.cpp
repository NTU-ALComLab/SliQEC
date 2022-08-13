#include "bddSystem.h"

/*
Apply gate on the bit slicing BDD system.
*/

void BDDSystem::Toffoli(int ithCircuit, int targ, std::vector<int> cont, std::vector<int> ncont)
{
    assert((cont.size() + ncont.size()) < _n);
    int IsBadtarg = 0;
    int cont_tot = cont.size() + ncont.size();
    for (int i = 0; i < cont_tot; i++)
    {
        if (i < cont.size())
        {
            if (targ == cont[i])
            {
                IsBadtarg = 1;
                break;
            }
        }
        else
        {
            if (targ == ncont[i - cont.size()])
            {
                IsBadtarg = 1;
                break;
            }
        }
    }
    assert(!IsBadtarg);

    DdNode *term1, *term2, *term3, *g, *tmp;

    g = Cudd_ReadOne(_ddManager);
    Cudd_Ref(g);
    for (int h = cont.size() - 1; h >= 0; h--)
    {
        tmp = Cudd_bddAnd(_ddManager, Cudd_bddIthVar(_ddManager, cont[h]), g);
        Cudd_Ref(tmp);
        Cudd_RecursiveDeref(_ddManager, g);
        g = tmp;
    }
    for (int h = ncont.size() - 1; h >= 0; h--)
    {
        tmp = Cudd_bddAnd(_ddManager, Cudd_Not(Cudd_bddIthVar(_ddManager, ncont[h])), g);
        Cudd_Ref(tmp);
        Cudd_RecursiveDeref(_ddManager, g);
        g = tmp;
    }

    for (int i = 0; i < _w; i++) // F = _allBDD[ithCircuit][i][j]
    {
        for (int j = 0; j < _r; j++)
        {
            // term1
            term1 = Cudd_ReadOne(_ddManager);
            Cudd_Ref(term1);
            tmp = Cudd_bddAnd(_ddManager, _allBDD[ithCircuit][i][j], term1);
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, term1);
            term1 = tmp;
            tmp = Cudd_bddAnd(_ddManager, Cudd_Not(g), term1);
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, term1);
            term1 = tmp;

            // term2
            term2 = Cudd_Cofactor(_ddManager, _allBDD[ithCircuit][i][j], Cudd_Not(Cudd_bddIthVar(_ddManager, targ)));
            Cudd_Ref(term2);

            tmp = Cudd_Cofactor(_ddManager, term2, g);
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, term2);
            term2 = tmp;

            tmp = Cudd_bddAnd(_ddManager, term2, g);
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, term2);
            term2 = tmp;

            tmp = Cudd_bddAnd(_ddManager, term2, Cudd_bddIthVar(_ddManager, targ));
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, term2);
            term2 = tmp;

            // term3
            term3 = Cudd_Cofactor(_ddManager, _allBDD[ithCircuit][i][j], Cudd_bddIthVar(_ddManager, targ));
            Cudd_Ref(term3);
            Cudd_RecursiveDeref(_ddManager, _allBDD[ithCircuit][i][j]);

            tmp = Cudd_Cofactor(_ddManager, term3, g);
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, term3);
            term3 = tmp;

            tmp = Cudd_bddAnd(_ddManager, term3, g);
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, term3);
            term3 = tmp;

            tmp = Cudd_bddAnd(_ddManager, term3, Cudd_Not(Cudd_bddIthVar(_ddManager, targ)));
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, term3);
            term3 = tmp;

            // OR
            _allBDD[ithCircuit][i][j] = Cudd_bddOr(_ddManager, term1, term2);
            Cudd_Ref(_allBDD[ithCircuit][i][j]);
            Cudd_RecursiveDeref(_ddManager, term1);
            Cudd_RecursiveDeref(_ddManager, term2);
            tmp = Cudd_bddOr(_ddManager, term3, _allBDD[ithCircuit][i][j]);
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, term3);
            Cudd_RecursiveDeref(_ddManager, _allBDD[ithCircuit][i][j]);
            _allBDD[ithCircuit][i][j] = tmp;
        }
    }
    Cudd_RecursiveDeref(_ddManager, g);
    updateNodeCount();
}

void BDDSystem::Fredkin(int ithCircuit, int swapA , int swapB, std::vector<int> cont)
{
    assert(cont.size() < _n);
    int IsBadtarg = 0;
    for (int i = 0; i < cont.size(); i++)
    {
        if ((swapA == cont[i]) || (swapB == cont[i]))
        {
            IsBadtarg = 1;
            break;
        }
    }
    assert(!IsBadtarg);

    DdNode *term1, *term2, *term3, *g, *tmp, *tmp0;

    g = Cudd_ReadOne(_ddManager);
    Cudd_Ref(g);
    for (int h = cont.size() - 1; h >= 0; h--)
    {
        tmp = Cudd_bddAnd(_ddManager, Cudd_bddIthVar(_ddManager, cont[h]), g);
        Cudd_Ref(tmp);
        Cudd_RecursiveDeref(_ddManager, g);
        g = tmp;
    }

    for (int i = 0; i < _w; i++) // F = _allBDD[ithCircuit][i][j]
    {
        for (int j = 0; j < _r; j++)
        {
            // term1
            term1 = Cudd_ReadOne(_ddManager);
            Cudd_Ref(term1);
            tmp = Cudd_bddAnd(_ddManager, _allBDD[ithCircuit][i][j], term1);
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, term1);
            term1 = tmp;
            tmp = Cudd_bddXor(_ddManager, Cudd_bddIthVar(_ddManager, swapA), Cudd_bddIthVar(_ddManager, swapB));
            Cudd_Ref(tmp);
            tmp0 = Cudd_Not(Cudd_bddAnd(_ddManager, g, tmp));
            Cudd_Ref(tmp0);
            Cudd_RecursiveDeref(_ddManager, tmp);
            tmp = Cudd_bddAnd(_ddManager, term1, tmp0);
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, term1);
            Cudd_RecursiveDeref(_ddManager, tmp0);
            term1 = tmp;

            // term2
            term2 = Cudd_Cofactor(_ddManager, _allBDD[ithCircuit][i][j], g);
            Cudd_Ref(term2);

            tmp = Cudd_Cofactor(_ddManager, term2, Cudd_Not(Cudd_bddIthVar(_ddManager, swapA)));
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, term2);
            term2 = tmp;

            tmp = Cudd_Cofactor(_ddManager, term2, Cudd_bddIthVar(_ddManager, swapB));
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, term2);
            term2 = tmp;

            tmp = Cudd_bddAnd(_ddManager, term2, g);
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, term2);
            term2 = tmp;

            tmp = Cudd_bddAnd(_ddManager, term2, Cudd_bddIthVar(_ddManager, swapA));
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, term2);
            term2 = tmp;

            tmp = Cudd_bddAnd(_ddManager, term2, Cudd_Not(Cudd_bddIthVar(_ddManager, swapB)));
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, term2);
            term2 = tmp;

            // term3
            term3 = Cudd_Cofactor(_ddManager, _allBDD[ithCircuit][i][j], g);
            Cudd_Ref(term3);

            tmp = Cudd_Cofactor(_ddManager, term3, Cudd_bddIthVar(_ddManager, swapA));
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, term3);
            term3 = tmp;

            tmp = Cudd_Cofactor(_ddManager, term3, Cudd_Not(Cudd_bddIthVar(_ddManager, swapB)));
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, term3);
            term3 = tmp;

            tmp = Cudd_bddAnd(_ddManager, term3, g);
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, term3);
            term3 = tmp;

            tmp = Cudd_bddAnd(_ddManager, term3, Cudd_Not(Cudd_bddIthVar(_ddManager, swapA)));
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, term3);
            term3 = tmp;

            tmp = Cudd_bddAnd(_ddManager, term3, Cudd_bddIthVar(_ddManager, swapB));
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, term3);
            term3 = tmp;

            // OR
            _allBDD[ithCircuit][i][j] = Cudd_bddOr(_ddManager, term1, term2);
            Cudd_Ref(_allBDD[ithCircuit][i][j]);
            Cudd_RecursiveDeref(_ddManager, term1);
            Cudd_RecursiveDeref(_ddManager, term2);
            tmp = Cudd_bddOr(_ddManager, term3, _allBDD[ithCircuit][i][j]);
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, term3);
            Cudd_RecursiveDeref(_ddManager, _allBDD[ithCircuit][i][j]);
            _allBDD[ithCircuit][i][j] = tmp;
        }
    }
    Cudd_RecursiveDeref(_ddManager, g);
    updateNodeCount();
}

void BDDSystem::Hadamard(int ithCircuit, int iqubit)
{
    assert((iqubit >= 0) & (iqubit < 2*_n));

    _k[ithCircuit] = _k[ithCircuit] + 1;

    DdNode *g, *d, *c, *tmp, *term1, *term2;

    int overflow_done = 0;

    for (int i = 0; i < _w; i++) // F = _allBDD[ithCircuit][i][j]
    {
        c = Cudd_ReadOne(_ddManager); 
        Cudd_Ref(c);
        tmp = Cudd_bddAnd(_ddManager, c, Cudd_bddIthVar(_ddManager, iqubit));
        Cudd_Ref(tmp);
        Cudd_RecursiveDeref(_ddManager, c);
        c = tmp;
        for (int j = 0; j < _r; j++)
        {
            // G
            g = Cudd_Cofactor(_ddManager, _allBDD[ithCircuit][i][j], Cudd_Not(Cudd_bddIthVar(_ddManager, iqubit)));
            Cudd_Ref(g);

            // D
            term1 = Cudd_Cofactor(_ddManager, _allBDD[ithCircuit][i][j], Cudd_bddIthVar(_ddManager, iqubit));
            Cudd_Ref(term1);
            tmp = Cudd_bddAnd(_ddManager, term1, Cudd_Not(Cudd_bddIthVar(_ddManager, iqubit)));
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, term1);
            term1 = tmp;
            term2 = Cudd_Not(_allBDD[ithCircuit][i][j]);
            Cudd_Ref(term2);
            tmp = Cudd_bddAnd(_ddManager, term2, Cudd_bddIthVar(_ddManager, iqubit));
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, term2);
            term2 = tmp;
            d = Cudd_bddOr(_ddManager, term1, term2);
            Cudd_Ref(d);
            Cudd_RecursiveDeref(_ddManager, term1);
            Cudd_RecursiveDeref(_ddManager, term2);

            // Detect overflow
            if ((j == _r - 1) && !overflow_done)
                if (overflow3(g, d, c))
                {
                    _r += _inc;
                    for (int indexOfCircuit = 0; indexOfCircuit < _nCircuit; indexOfCircuit++)   // each BDD of each circuit should be extend together
                        allocBDD(_allBDD[indexOfCircuit], true); // add new BDDs
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
    updateNodeCount();
}

void BDDSystem::rx_pi_2(int ithCircuit, int iqubit, bool dagger)
{
    assert((iqubit >= 0) & (iqubit < 2*_n));

    _k[ithCircuit] = _k[ithCircuit] + 1;

    int nshift = _w / 2;
    int overflow_done = 0;

    DdNode *g, *d, *c, *tmp, *term1, *term2;
    DdNode ***copy = new DdNode **[_w];
    for (int i = 0; i < _w; i++)
        copy[i] = new DdNode *[_r];

    // Copy
    for (int i = 0; i < _w; i++)
    {
         for (int j = 0; j < _r; j++)
        {
            copy[i][j] = Cudd_ReadOne(_ddManager);
            Cudd_Ref(copy[i][j]);
            tmp = Cudd_bddAnd(_ddManager, copy[i][j], _allBDD[ithCircuit][i][j]);
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, copy[i][j]);
            copy[i][j] = tmp;
        }
    }

    for (int i = 0; i < _w; i++)
    {
        // Init C
        if ((i < nshift) ^ dagger)
            c = Cudd_ReadOne(_ddManager);
        else
            c = Cudd_Not(Cudd_ReadOne(_ddManager));
        Cudd_Ref(c);
        for (int j = 0; j < _r; j++)
        {
            // D
            term1 = Cudd_Cofactor(_ddManager, copy[(i + nshift) % _w][j], Cudd_Not(Cudd_bddIthVar(_ddManager, iqubit)));
            Cudd_Ref(term1);
            tmp = Cudd_bddAnd(_ddManager, term1, Cudd_bddIthVar(_ddManager, iqubit));
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, term1);
            term1 = tmp;
            term2 = Cudd_Cofactor(_ddManager, copy[(i + nshift) % _w][j], Cudd_bddIthVar(_ddManager, iqubit));
            Cudd_Ref(term2);
            tmp = Cudd_bddAnd(_ddManager, term2, Cudd_Not(Cudd_bddIthVar(_ddManager, iqubit)));
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, term2);
            term2 = tmp;
            if ((i < nshift) ^ dagger) d = Cudd_Not(Cudd_bddOr(_ddManager, term1, term2));
            else d = Cudd_bddOr(_ddManager, term1, term2);
            Cudd_Ref(d);
            Cudd_RecursiveDeref(_ddManager, term1);
            Cudd_RecursiveDeref(_ddManager, term2);
            // Detect overflow
            if ((j == _r - 1) && !overflow_done)
                if (overflow3(copy[i][j], d, c))
                {
                    _r += _inc;
                    for (int indexOfCircuit = 0; indexOfCircuit < _nCircuit; indexOfCircuit++)   // each BDD of each circuit should be extend together
                        allocBDD(_allBDD[indexOfCircuit], true); // add new BDDs
                    allocBDD(copy, true);
                    overflow_done = 1;
                }
            // Sum
            Cudd_RecursiveDeref(_ddManager, _allBDD[ithCircuit][i][j]);
            _allBDD[ithCircuit][i][j] = Cudd_bddXor(_ddManager, copy[i][j], d);
            Cudd_Ref(_allBDD[ithCircuit][i][j]);
            tmp = Cudd_bddXor(_ddManager, _allBDD[ithCircuit][i][j], c);
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, _allBDD[ithCircuit][i][j]);
            _allBDD[ithCircuit][i][j] = tmp;
            // Carry
            if (j == _r - 1)
            {
                Cudd_RecursiveDeref(_ddManager, c);
                Cudd_RecursiveDeref(_ddManager, d);
            }
            else
            {
                term1 = Cudd_bddAnd(_ddManager, copy[i][j], d);
                Cudd_Ref(term1);
                term2 = Cudd_bddOr(_ddManager, copy[i][j], d);
                Cudd_Ref(term2);
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
    for (int i = 0; i < _w; i++)
    {
        for (int j = 0; j < _r; j++)
            Cudd_RecursiveDeref(_ddManager, copy[i][j]);
        delete[] copy[i];
    }
    updateNodeCount();
}

void BDDSystem::ry_pi_2(int ithCircuit, int iqubit, bool transpose)
{

    assert((iqubit >= 0) & (iqubit < 2*_n));

    _k[ithCircuit] = _k[ithCircuit] + 1;

    int overflow_done = 0;

    DdNode *g, *d, *c, *tmp, *term1, *term2, *var;

    if (transpose) var = Cudd_Not(Cudd_bddIthVar(_ddManager, iqubit));
    else var = Cudd_bddIthVar(_ddManager, iqubit);

    for (int i = 0; i < _w; i++)
    {
        // Init C
        c = Cudd_ReadOne(_ddManager); 
        Cudd_Ref(c);
        tmp = Cudd_bddAnd(_ddManager, c, Cudd_Not(var));
        Cudd_Ref(tmp);
        Cudd_RecursiveDeref(_ddManager, c);
        c = tmp;
        for (int j = 0; j < _r; j++)
        {
            // G
            g = Cudd_Cofactor(_ddManager, _allBDD[ithCircuit][i][j], Cudd_Not(var));
            Cudd_Ref(g);
            // D
            term1 = Cudd_bddAnd(_ddManager, _allBDD[ithCircuit][i][j], var);
            Cudd_Ref(term1);
            term2 = Cudd_Not(Cudd_Cofactor(_ddManager, _allBDD[ithCircuit][i][j], var));
            Cudd_Ref(term2);
            tmp = Cudd_bddAnd(_ddManager, term2, Cudd_Not(var));
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, term2);
            term2 = tmp;
            d = Cudd_bddOr(_ddManager, term1, term2);
            Cudd_Ref(d);
            Cudd_RecursiveDeref(_ddManager, term1);
            Cudd_RecursiveDeref(_ddManager, term2);

            // Detect overflow
            if ((j == _r - 1) && !overflow_done)
                if (overflow3(g, d, c))
                {
                    _r += _inc;
                    for (int indexOfCircuit = 0; indexOfCircuit < _nCircuit; indexOfCircuit++)   // each BDD of each circuit should be extend together
                        allocBDD(_allBDD[indexOfCircuit], true); // add new BDDs
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
    updateNodeCount();
}

void BDDSystem::Phase_shift(int ithCircuit, int phase, int iqubit)
{
    assert((iqubit >= 0) & (iqubit < 2*_n));

    int nshift = _w / phase;
    int overflow_done = 0;

    DdNode *g, *c, *tmp, *term1, *term2;

    // Copy
    DdNode **copy[_w];
    for (int i = 0; i < _w; i++)
        copy[i] = new DdNode *[_r];
    for (int i = 0; i < _w; i++)
    {
         for (int j = 0; j < _r; j++)
        {
            copy[i][j] = Cudd_ReadOne(_ddManager);
            Cudd_Ref(copy[i][j]);
            tmp = Cudd_bddAnd(_ddManager, copy[i][j], _allBDD[ithCircuit][i][j]);
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, copy[i][j]);
            copy[i][j] = tmp;
        }
    }

    for (int i = 0; i < _w; i++)
    {
        // Init C
        if (i >= _w - nshift)
        {
            c = Cudd_bddIthVar(_ddManager, iqubit);
            Cudd_Ref(c);
        }

        for (int j = 0; j < _r; j++)
        {
            if (i >= _w - nshift)
            {
                term1 = Cudd_bddAnd(_ddManager, copy[i][j], Cudd_Not(Cudd_bddIthVar(_ddManager, iqubit)));
                Cudd_Ref(term1);
                term2 = Cudd_bddAnd(_ddManager, Cudd_Not(copy[i - (_w - nshift)][j]), Cudd_bddIthVar(_ddManager, iqubit));
                Cudd_Ref(term2);
                g = Cudd_bddOr(_ddManager, term1, term2);
                Cudd_Ref(g);
                Cudd_RecursiveDeref(_ddManager, term1);
                Cudd_RecursiveDeref(_ddManager, term2);

                // Detect overflow
                if ((j == _r - 1) && !overflow_done)
                    if (overflow2(g, c))
                    {
                        _r += _inc;
                        for (int indexOfCircuit = 0; indexOfCircuit < _nCircuit; indexOfCircuit++)   // each BDD of each circuit should be extend together
                            allocBDD(_allBDD[indexOfCircuit], true); // add new BDDs
                        allocBDD(copy, true);      // add new BDDs
                        overflow_done = 1;
                    }

                // Plus
                Cudd_RecursiveDeref(_ddManager, _allBDD[ithCircuit][i][j]);
                if (Cudd_IsConstant(c))     // must be constant 0
                    _allBDD[ithCircuit][i][j] = g;
                else
                {
                    // Sum
                    _allBDD[ithCircuit][i][j] = Cudd_bddXor(_ddManager, g, c);
                    Cudd_Ref(_allBDD[ithCircuit][i][j]);
                    // Carry
                    if (j == _r - 1)
                    {
                        Cudd_RecursiveDeref(_ddManager, g);
                        Cudd_RecursiveDeref(_ddManager, c);
                    }
                    else
                    {
                        tmp = Cudd_bddAnd(_ddManager, g, c);
                        Cudd_Ref(tmp);
                        Cudd_RecursiveDeref(_ddManager, g);
                        Cudd_RecursiveDeref(_ddManager, c);
                        c = tmp;
                    }
                }
            }
            else
            {
                term1 = Cudd_bddAnd(_ddManager, copy[i][j], Cudd_Not(Cudd_bddIthVar(_ddManager, iqubit)));
                Cudd_Ref(term1);
                term2 = Cudd_bddAnd(_ddManager, copy[i + nshift][j], Cudd_bddIthVar(_ddManager, iqubit));
                Cudd_Ref(term2);

                Cudd_RecursiveDeref(_ddManager, _allBDD[ithCircuit][i][j]);
                _allBDD[ithCircuit][i][j] = Cudd_bddOr(_ddManager, term1, term2);
                Cudd_Ref(_allBDD[ithCircuit][i][j]);

                Cudd_RecursiveDeref(_ddManager, term1);
                Cudd_RecursiveDeref(_ddManager, term2);
            }
        }
    }

    for (int i = 0; i < _w; i++)
    {
        for (int j = 0; j < _r; j++)
            Cudd_RecursiveDeref(_ddManager, copy[i][j]);
        delete[] copy[i];
    }
    updateNodeCount();
}

void BDDSystem::Phase_shift_dagger(int ithCircuit, int phase, int iqubit)
{
    assert((iqubit >= 0) & (iqubit < 2*_n));

    int nshift = _w / abs(phase);
    int overflow_done = 0;

    DdNode *g, *c, *tmp, *term1, *term2;

    // Copy
    DdNode **copy[_w];
    for (int i = 0; i < _w; i++)
        copy[i] = new DdNode *[_r];
    for (int i = 0; i < _w; i++)
    {
         for (int j = 0; j < _r; j++)
        {
            copy[i][j] = Cudd_ReadOne(_ddManager);
            Cudd_Ref(copy[i][j]);
            tmp = Cudd_bddAnd(_ddManager, copy[i][j], _allBDD[ithCircuit][i][j]);
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, copy[i][j]);
            copy[i][j] = tmp;
        }
    }

    for (int i = 0; i < _w; i++)
    {
        // Init C
        if (i < nshift)
        {
            c = Cudd_bddIthVar(_ddManager, iqubit);
            Cudd_Ref(c);
        }

        for (int j = 0; j < _r; j++)
        {
            if (i < nshift)
            {
                term1 = Cudd_bddAnd(_ddManager, copy[i][j], Cudd_Not(Cudd_bddIthVar(_ddManager, iqubit)));
                Cudd_Ref(term1);
                term2 = Cudd_bddAnd(_ddManager, Cudd_Not(copy[_w - nshift + i][j]), Cudd_bddIthVar(_ddManager, iqubit));
                Cudd_Ref(term2);
                g = Cudd_bddOr(_ddManager, term1, term2);
                Cudd_Ref(g);
                Cudd_RecursiveDeref(_ddManager, term1);
                Cudd_RecursiveDeref(_ddManager, term2);

                // Detect overflow
                if ((j == _r - 1) && !overflow_done)
                    if (overflow2(g, c))
                    {
                        _r += _inc;
                        for (int indexOfCircuit = 0; indexOfCircuit < _nCircuit; indexOfCircuit++)   // each BDD of each circuit should be extend together
                            allocBDD(_allBDD[indexOfCircuit], true); // add new BDDs
                        allocBDD(copy, true);      // add new BDD
                        overflow_done = 1;
                    }

                // Plus
                Cudd_RecursiveDeref(_ddManager, _allBDD[ithCircuit][i][j]);
                if (Cudd_IsConstant(c))     // must be constant 0
                    _allBDD[ithCircuit][i][j] = g;
                else
                {
                    // Sum
                    _allBDD[ithCircuit][i][j] = Cudd_bddXor(_ddManager, g, c);
                    Cudd_Ref(_allBDD[ithCircuit][i][j]);
                    // Carry
                    if (j == _r - 1)
                    {
                        Cudd_RecursiveDeref(_ddManager, g);
                        Cudd_RecursiveDeref(_ddManager, c);
                    }
                    else
                    {
                        tmp = Cudd_bddAnd(_ddManager, g, c);
                        Cudd_Ref(tmp);
                        Cudd_RecursiveDeref(_ddManager, g);
                        Cudd_RecursiveDeref(_ddManager, c);
                        c = tmp;
                    }
                }
            }
            else
            {
                term1 = Cudd_bddAnd(_ddManager, copy[i][j], Cudd_Not(Cudd_bddIthVar(_ddManager, iqubit)));
                Cudd_Ref(term1);
                term2 = Cudd_bddAnd(_ddManager, copy[i - nshift][j], Cudd_bddIthVar(_ddManager, iqubit));
                Cudd_Ref(term2);

                Cudd_RecursiveDeref(_ddManager, _allBDD[ithCircuit][i][j]);
                _allBDD[ithCircuit][i][j] = Cudd_bddOr(_ddManager, term1, term2);
                Cudd_Ref(_allBDD[ithCircuit][i][j]);

                Cudd_RecursiveDeref(_ddManager, term1);
                Cudd_RecursiveDeref(_ddManager, term2);
            }
        }
    }

    for (int i = 0; i < _w; i++)
    {
        for (int j = 0; j < _r; j++)
            Cudd_RecursiveDeref(_ddManager, copy[i][j]);
        delete[] copy[i];
    }
    updateNodeCount();
}

void BDDSystem::PauliX(int ithCircuit, int iqubit)
{
    assert((iqubit >= 0) & (iqubit < 2*_n));

    DdNode *tmp, *term1, *term2;

    for (int i = 0; i < _w; i++) // F = _allBDD[ithCircuit][i][j]
    {
        for (int j = 0; j < _r; j++)
        {
            // term1
            term1 = Cudd_Cofactor(_ddManager, _allBDD[ithCircuit][i][j], Cudd_Not(Cudd_bddIthVar(_ddManager, iqubit)));
            Cudd_Ref(term1);

            tmp = Cudd_bddAnd(_ddManager, term1, Cudd_bddIthVar(_ddManager, iqubit));
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, term1);
            term1 = tmp;

            // term2
            term2 = Cudd_Cofactor(_ddManager, _allBDD[ithCircuit][i][j], Cudd_bddIthVar(_ddManager, iqubit));
            Cudd_Ref(term2);
            Cudd_RecursiveDeref(_ddManager, _allBDD[ithCircuit][i][j]);

            tmp = Cudd_bddAnd(_ddManager, term2, Cudd_Not(Cudd_bddIthVar(_ddManager, iqubit)));
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, term2);
            term2 = tmp;

            // OR
            _allBDD[ithCircuit][i][j] = Cudd_bddOr(_ddManager, term1, term2);
            Cudd_Ref(_allBDD[ithCircuit][i][j]);
            Cudd_RecursiveDeref(_ddManager, term1);
            Cudd_RecursiveDeref(_ddManager, term2);
        }
    }
    updateNodeCount();
}

void BDDSystem::PauliY(int ithCircuit, int iqubit, bool transpose)
{
    assert((iqubit >= 0) & (iqubit < 2*_n));

    int nshift = _w / 2;

    DdNode *g, *c, *tmp, *term1, *term2, *var;
    int overflow_done = 0;

    if (transpose) var = Cudd_Not(Cudd_bddIthVar(_ddManager, iqubit));
    else var = Cudd_bddIthVar(_ddManager, iqubit);

    // PauliX(iqubit);
    for (int i = 0; i < _w; i++) // F = _allBDD[ithCircuit][i][j]
    {
        for (int j = 0; j < _r; j++)
        {
            // term1
            term1 = Cudd_Cofactor(_ddManager, _allBDD[ithCircuit][i][j], Cudd_Not(var));
            Cudd_Ref(term1);

            tmp = Cudd_bddAnd(_ddManager, term1, var);
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, term1);
            term1 = tmp;

            // term2
            term2 = Cudd_Cofactor(_ddManager, _allBDD[ithCircuit][i][j], var);
            Cudd_Ref(term2);
            Cudd_RecursiveDeref(_ddManager, _allBDD[ithCircuit][i][j]);

            tmp = Cudd_bddAnd(_ddManager, term2, Cudd_Not(var));
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, term2);
            term2 = tmp;

            // OR
            _allBDD[ithCircuit][i][j] = Cudd_bddOr(_ddManager, term1, term2);
            Cudd_Ref(_allBDD[ithCircuit][i][j]);
            Cudd_RecursiveDeref(_ddManager, term1);
            Cudd_RecursiveDeref(_ddManager, term2);
        }
    }

    // Copy
    DdNode **copy[_w];
    for (int i = 0; i < _w; i++)
        copy[i] = new DdNode *[_r];
    for (int i = 0; i < _w; i++)
    {
         for (int j = 0; j < _r; j++)
        {
            copy[i][j] = Cudd_ReadOne(_ddManager);
            Cudd_Ref(copy[i][j]);
            tmp = Cudd_bddAnd(_ddManager, copy[i][j], _allBDD[ithCircuit][i][j]);
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, copy[i][j]);
            copy[i][j] = tmp;
        }
    }

    for (int i = 0; i < _w; i++)
    {
        // Init C
        if (i < nshift)
            c = Cudd_Not(var);
        else
            c = Cudd_bddAnd(_ddManager, Cudd_ReadOne(_ddManager), var);
        Cudd_Ref(c);

        for (int j = 0; j < _r; j++)
        {
            if (i < nshift)
            {
                term1 = Cudd_bddAnd(_ddManager, copy[i + nshift][j], var);
                Cudd_Ref(term1);
                term2 = Cudd_bddAnd(_ddManager, Cudd_Not(copy[i + nshift][j]), Cudd_Not(var));
                Cudd_Ref(term2);
            }
            else
            {
                term1 = Cudd_bddAnd(_ddManager, Cudd_Not(copy[i - nshift][j]), var);
                Cudd_Ref(term1);
                term2 = Cudd_bddAnd(_ddManager, copy[i - nshift][j], Cudd_Not(var));
                Cudd_Ref(term2);
            }
            g = Cudd_bddOr(_ddManager, term1, term2);
            Cudd_Ref(g);
            Cudd_RecursiveDeref(_ddManager, term1);
            Cudd_RecursiveDeref(_ddManager, term2);

            // Detect overflow
            if ((j == _r - 1) && !overflow_done)
                if (overflow2(g, c))
                {
                    _r += _inc;
                    for (int indexOfCircuit = 0; indexOfCircuit < _nCircuit; indexOfCircuit++)   // each BDD of each circuit should be extend together
                        allocBDD(_allBDD[indexOfCircuit], true); // add new BDDs
                    allocBDD(copy, true);
                    overflow_done = 1;
                }

            // Plus 1
            Cudd_RecursiveDeref(_ddManager, _allBDD[ithCircuit][i][j]);
            if (Cudd_IsConstant(c))
                _allBDD[ithCircuit][i][j] = g;
            else
            {
                // Sum
                _allBDD[ithCircuit][i][j] = Cudd_bddXor(_ddManager, g, c);
                Cudd_Ref(_allBDD[ithCircuit][i][j]);
                // Carry
                if (j == _r - 1)
                {
                    Cudd_RecursiveDeref(_ddManager, g);
                    Cudd_RecursiveDeref(_ddManager, c);
                }
                else
                {
                    tmp = Cudd_bddAnd(_ddManager, g, c);
                    Cudd_Ref(tmp);
                    Cudd_RecursiveDeref(_ddManager, g);
                    Cudd_RecursiveDeref(_ddManager, c);
                    c = tmp;
                }
            }
        }
    }

    for (int i = 0; i < _w; i++)
    {
        for (int j = 0; j < _r; j++)
            Cudd_RecursiveDeref(_ddManager, copy[i][j]);
        delete[] copy[i];
    }
    updateNodeCount();
}

void BDDSystem::PauliZ(int ithCircuit, std::vector<int> iqubit)
{
    for (int i = 0; i < iqubit.size(); i++)
    {
        assert((iqubit[i] >= 0) & (iqubit[i] < 2*_n));
    }
    assert((iqubit.size() == 1) || (iqubit.size() == 2));

    DdNode *c, *tmp, *term1, *term2, *inter, *qubit_and;

    // Init qubit and
    qubit_and = Cudd_ReadOne(_ddManager); 
    Cudd_Ref(qubit_and);
    for (int i = iqubit.size() - 1; i >= 0; i--)
    {
        tmp = Cudd_bddAnd(_ddManager, qubit_and, Cudd_bddIthVar(_ddManager, iqubit[i]));
        Cudd_Ref(tmp);
        Cudd_RecursiveDeref(_ddManager, qubit_and);
        qubit_and = tmp;
    }
    int overflow_done = 0;
    for (int i = 0; i < _w; i++)
    {
        // Init C
        c = Cudd_ReadOne(_ddManager); 
        Cudd_Ref(c);
        tmp = Cudd_bddAnd(_ddManager, c, qubit_and);
        Cudd_Ref(tmp);
        Cudd_RecursiveDeref(_ddManager, c);
        c = tmp;
        for (int j = 0; j < _r; j++)
        {
            term1 = Cudd_bddAnd(_ddManager, _allBDD[ithCircuit][i][j], Cudd_Not(qubit_and));
            Cudd_Ref(term1);
            term2 = Cudd_Not(_allBDD[ithCircuit][i][j]);
            Cudd_Ref(term2);
            Cudd_RecursiveDeref(_ddManager, _allBDD[ithCircuit][i][j]);
            tmp = Cudd_bddAnd(_ddManager, term2, qubit_and);
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, term2);
            term2 = tmp;
            inter = Cudd_bddOr(_ddManager, term1, term2);
            Cudd_Ref(inter);
            Cudd_RecursiveDeref(_ddManager, term1);
            Cudd_RecursiveDeref(_ddManager, term2);

            // Plus 1
            if (Cudd_IsConstant(c))
                _allBDD[ithCircuit][i][j] = inter;
            else
            {
                // Detect overflow
                if ((i == _r - 1) && !overflow_done)
                    if (overflow2(inter, c))
                    {
                        _r += _inc;
                        for (int indexOfCircuit = 0; indexOfCircuit < _nCircuit; indexOfCircuit++)   // each BDD of each circuit should be extend together
                            allocBDD(_allBDD[indexOfCircuit], true); // add new BDDs
                        overflow_done = 1;
                    }

                // Sum
                _allBDD[ithCircuit][i][j] = Cudd_bddXor(_ddManager, inter, c);
                Cudd_Ref(_allBDD[ithCircuit][i][j]);
                // Carry
                if (i == _r - 1)
                    Cudd_RecursiveDeref(_ddManager, inter);
                else
                {
                    tmp = Cudd_bddAnd(_ddManager, inter, c);
                    Cudd_Ref(tmp);
                    Cudd_RecursiveDeref(_ddManager, c);
                    Cudd_RecursiveDeref(_ddManager, inter);
                    c = tmp;
                }
            }
        }
        Cudd_RecursiveDeref(_ddManager, c);
    }
    Cudd_RecursiveDeref(_ddManager, qubit_and);
    updateNodeCount();
}
