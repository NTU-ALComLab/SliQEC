#include "util_sim.h"
#include "Simulator.h" // exp

std::ofstream outFile; // exp

/**Function*************************************************************

  Synopsis    []

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void full_adder_plus_1(int length, int *reg)
{
    int one = 1, carry = 0;
    for (int i = 0; i < length; i++)
    {
        carry = reg[i] & one;
        reg[i] = reg[i] ^ one;

        if (carry == 0)
            break;
    }
}

/**Function*************************************************************

  Synopsis    []

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void full_adder_plus_1_start(int length, int *reg, int start)
{
    int one = 1, carry = 0;
    for (int i = start; i < length; i++)
    {
        carry = reg[i] & one;
        reg[i] = reg[i] ^ one;

        if (carry == 0)
            break;
    }
}

/**Function*************************************************************

  Synopsis    []

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void full_adder_plus_1_measure(int length, int *reg, int *order)
{
    int one = 1, carry = 0;
    for (int i = length - 1; i >= 0; i--)
    {
        carry = reg[order[i]] & one;
        reg[order[i]] = reg[order[i]] ^ one;

        if (carry == 0)
            break;
    }
}

/**Function*************************************************************

  Synopsis    []

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
int int_array_full_check(int length, int *reg)
{
    int check = 1;
    for (int i = 0; i < length; i++)
        check *= reg[i];

    return check;
}


/**Function*************************************************************

  Synopsis    []

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void signalHandler(int signum)
{
    std::cout << "Interrupt signal (" << signum << ") received.\n";

    outFile  << "TO/MO" << std::endl;
    outFile.close();

    // terminate program
    exit(signum);
}
