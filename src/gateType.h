#ifndef _GATETYPE_H_
#define _GATETYPE_H_

enum class GateType
{
    X,              // Pauli-X
    Y,              // Pauli-Y
    Z,              // Pauli-Z
    H,              // Hadamard
    S,              // Phase 
    SDG,            // Phase inverse
    T,              // Pi/8
    TDG,            // Pi/8 inverse
    RX_PI_2,        // Rotation-X Pi/2
    RX_PI_2_DG,     // Rotation-X Pi/2 inverse
    RY_PI_2,        // Rotation-Y Pi/2 
    RY_PI_2_DG,     // Rotation-Y Pi/2 inverse
    CX,             // Controlled-NOT
    CZ,             // Controlled-Z
    CCX,            // Toffoli
    SWAP,           // SWAP
    CSWAP           // Fredkin
};

#endif
