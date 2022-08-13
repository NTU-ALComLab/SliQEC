#include <string>
#include <sstream>
#include <vector>
#include <cassert>
#include <iostream>
#include <fstream>

#include "gateType.h"

/**Function*************************************************************

  Synopsis    [parse QASM file and store information into gate/qubit/n]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/

void qasmParser(std::ifstream &inFile, std::vector<GateType> &gates, std::vector<std::vector<int> > &qubits, int &n)
{
    std::stringstream strStream;
    strStream << inFile.rdbuf();
    std::string fileStr = strStream.str(); 
    std::stringstream file_ss(fileStr);

    std::string inStr;
    
    while(getline(file_ss, inStr))
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
                
                n = stoi(inStr); 
            }
            else if (inStr == "creg"){;}
            else if (inStr == "OPENQASM"){;}
            else if (inStr == "include"){;}
            else if (inStr == "measure"){;}
            else
            {
                if (inStr == "x")
                {
                    gates.push_back(GateType::X);
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    int iQubit = stoi(inStr);
                    assert(iQubit < n);
                    qubits.push_back(std::vector<int>(1, iQubit));
                }
                else if (inStr == "y")
                {
                    gates.push_back(GateType::Y);
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    int iQubit = stoi(inStr);
                    assert(iQubit < n);
                    qubits.push_back(std::vector<int>(1, iQubit));
                }
                else if (inStr == "z")
                {
                    gates.push_back(GateType::Z);
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    int iQubit = stoi(inStr);
                    assert(iQubit < n);
                    qubits.push_back(std::vector<int>(1, iQubit));
                }
                else if (inStr == "h")
                {
                    gates.push_back(GateType::H);
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    int iQubit = stoi(inStr);
                    assert(iQubit < n);
                    qubits.push_back(std::vector<int>(1, iQubit));
                }
                else if (inStr == "s")
                {
                    gates.push_back(GateType::S);
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    int iQubit = stoi(inStr);
                    assert(iQubit < n);
                    qubits.push_back(std::vector<int>(1, iQubit));
                }
                else if (inStr == "sdg")
                {
                    gates.push_back(GateType::SDG);
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    int iQubit = stoi(inStr);
                    assert(iQubit < n);
                    qubits.push_back(std::vector<int>(1, iQubit));
                }
                else if (inStr == "t")
                {
                    gates.push_back(GateType::T);
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    int iQubit = stoi(inStr);
                    assert(iQubit < n);
                    qubits.push_back(std::vector<int>(1, iQubit));
                }
                else if (inStr == "tdg")
                {
                    gates.push_back(GateType::TDG);
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    int iQubit = stoi(inStr);
                    assert(iQubit < n);
                    qubits.push_back(std::vector<int>(1, iQubit));
                }
                else if (inStr == "rx(pi/2)")
                {
                    gates.push_back(GateType::RX_PI_2);
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    int iQubit = stoi(inStr);
                    assert(iQubit < n);
                    qubits.push_back(std::vector<int>(1, iQubit));
                }
                else if (inStr == "rx(-pi/2)")
                {
                    gates.push_back(GateType::RX_PI_2_DG);
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    int iQubit = stoi(inStr);
                    assert(iQubit < n);
                    qubits.push_back(std::vector<int>(1, iQubit));
                }
                else if (inStr == "ry(pi/2)")
                {
                    gates.push_back(GateType::RY_PI_2);
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    int iQubit = stoi(inStr);
                    assert(iQubit < n);
                    qubits.push_back(std::vector<int>(1, iQubit));
                }
                else if (inStr == "ry(-pi/2)")
                {
                    gates.push_back(GateType::RY_PI_2_DG);
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    int iQubit = stoi(inStr);
                    assert(iQubit < n);
                    qubits.push_back(std::vector<int>(1, iQubit));
                }
                else if (inStr == "cx")
                {
                    gates.push_back(GateType::CX);
                    std::vector<int> iQubitsList(2);
                    for (int i = 0; i < 2; i++)
                    {
                        getline(inStr_ss, inStr, '[');
                        getline(inStr_ss, inStr, ']');
                        iQubitsList[i] = stoi(inStr);
                        assert(iQubitsList[i] < n);
                    }
                    qubits.push_back(iQubitsList);
                    iQubitsList.clear();
                }
                else if (inStr == "cz")
                {
                    gates.push_back(GateType::CZ);
                    std::vector<int> iQubitsList(2);
                    for (int i = 0; i < 2; i++)
                    {
                        getline(inStr_ss, inStr, '[');
                        getline(inStr_ss, inStr, ']');
                        iQubitsList[i] = stoi(inStr);
                        assert(iQubitsList[i] < n);
                    }
                    qubits.push_back(iQubitsList);
                    iQubitsList.clear();
                }
                else if (inStr == "swap")
                {
                    gates.push_back(GateType::SWAP);
                    std::vector<int> iQubitsList(2);
                    for (int i = 0; i < 2; i++)
                    {
                        getline(inStr_ss, inStr, '[');
                        getline(inStr_ss, inStr, ']');
                        iQubitsList[i] = stoi(inStr);
                        assert(iQubitsList[i] < n);
                    }
                    qubits.push_back(iQubitsList);
                    iQubitsList.clear();
                }
                else if (inStr == "cswap")
                {
                    gates.push_back(GateType::CSWAP);
                    std::vector<int> iQubitsList(3);
                    for (int i = 0; i < 3; i++)
                    {
                        getline(inStr_ss, inStr, '[');
                        getline(inStr_ss, inStr, ']');
                        iQubitsList[i] = stoi(inStr);
                        assert(iQubitsList[i] < n);
                    }
                    qubits.push_back(iQubitsList);
                    iQubitsList.clear();
                }
                else if (inStr == "ccx" || inStr == "mcx")
                {
                    gates.push_back(GateType::CCX);
                    std::vector<int> iQubitsList(0);
                    getline(inStr_ss, inStr, '[');
                    while(getline(inStr_ss, inStr, ']'))
                    {
                        iQubitsList.push_back(stoi(inStr));
                        assert(iQubitsList.back() < n);
                        getline(inStr_ss, inStr, '[');
                    }
                    qubits.push_back(iQubitsList);
                    iQubitsList.clear();
                }
                else
                {
                    std::cerr << std::endl
                            << "[Warning]: Syntax \'" << inStr << "\' is not supported. The line is ignored ..." << std::endl;
                }
            }
        }
    }
}
