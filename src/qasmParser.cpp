#include <string>
#include <sstream>
#include <vector>
#include <cassert>
#include <iostream>
#include <fstream>
#include <map>
#include <stack>

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
        if (inStr.find_first_not_of("\t\n\r ") != std::string::npos)
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
                    gates.push_back(GateType::MCS);
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    int iQubit = stoi(inStr);
                    assert(iQubit < n);
                    qubits.push_back(std::vector<int>(1, iQubit));
                }
                else if (inStr == "sdg")
                {
                    gates.push_back(GateType::MCSDG);
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    int iQubit = stoi(inStr);
                    assert(iQubit < n);
                    qubits.push_back(std::vector<int>(1, iQubit));
                }
                else if (inStr == "t")
                {
                    gates.push_back(GateType::MCT);
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    int iQubit = stoi(inStr);
                    assert(iQubit < n);
                    qubits.push_back(std::vector<int>(1, iQubit));
                }
                else if (inStr == "tdg")
                {
                    gates.push_back(GateType::MCTDG);
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
                else if (inStr == "cs" || inStr == "mcs")
                {
                    gates.push_back(GateType::MCS);
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
                else if (inStr == "csdg" || inStr == "mcsdg")
                {
                    gates.push_back(GateType::MCSDG);
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
                else if (inStr == "ct" || inStr == "mct")
                {
                    gates.push_back(GateType::MCT);
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
                else if (inStr == "ctdg" || inStr == "mctdg")
                {
                    gates.push_back(GateType::MCTDG);
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


/**Function*************************************************************

  Synopsis    [parse an infix Boolean function to postfix]

  Description [use shunting yard algorithm]

  SideEffects []

  SeeAlso     []

***********************************************************************/

std::vector<std::string> booleanParser(std::string& inStr)
{
    std::vector<std::string> result;
    std::map<std::string, int> operator_priority = { {"not", 4},
                                                     {"xor", 3},
                                                     {"and", 2},
                                                     {"or",  1} };
        
    std::stringstream inStr_ss(inStr);
    std::stack<std::string> waiting_operators;
    int nLeftTuple = 0;
    while(getline(inStr_ss, inStr, ' '))
    {
        if (inStr == "(")
        {
            waiting_operators.push(inStr);
            nLeftTuple++;
        }
        else if (inStr == ")")
        {
            assert(nLeftTuple > 0);
            while (waiting_operators.top() != "(") 
            {
                result.emplace_back(waiting_operators.top());
                waiting_operators.pop();
            }
            waiting_operators.pop();
            nLeftTuple--;
        }
        else if (operator_priority.count(inStr) != 0)
        {
            while (!waiting_operators.empty() && waiting_operators.top() != "(" && operator_priority[inStr] <= operator_priority[waiting_operators.top()])
            {
                result.emplace_back(waiting_operators.top());
                waiting_operators.pop();
            }
            waiting_operators.push(inStr);
        }
        else
        {
            std::stringstream inStr_sss(inStr);
            getline(inStr_sss, inStr, '[');
            getline(inStr_sss, inStr, ']');
            assert(inStr != "");
            result.emplace_back(inStr);
        }
    }
    assert(nLeftTuple == 0);
    
    while (!waiting_operators.empty())
    {
        result.emplace_back(waiting_operators.top());
        waiting_operators.pop();
    }
    return result;
}


/**Function*************************************************************

  Synopsis    [parse careSet file and store information in postfix notation]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/

void careSetParser(std::ifstream &inFile, std::vector<std::string>& careSet)  
{
    std::stringstream strStream;
    strStream << inFile.rdbuf();
    std::string fileStr = strStream.str(); 
    std::stringstream file_ss(fileStr);
    std::string inStr;
    
    getline(file_ss, inStr);    // only a single line is allowed
    inStr = inStr.substr(0, inStr.find("//"));
    careSet = booleanParser(inStr);
}

/**Function*************************************************************

  Synopsis    [parse weight function file and store information in postfix notation]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/

void weightFunParser(std::ifstream &inFile, std::vector<std::pair<std::vector<std::string>, std::pair<int, int>>>& weightFun)
{
    std::stringstream strStream;
    strStream << inFile.rdbuf();
    std::string fileStr = strStream.str(); 
    std::stringstream file_ss(fileStr);
    std::string inStr;
    
    while (getline(file_ss, inStr))
    {  
        // <Boolean function>
        inStr = inStr.substr(0, inStr.find("//"));
        std::vector<std::string> function = booleanParser(inStr);
        
        // <int, int>
        getline(file_ss, inStr);
        inStr = inStr.substr(0, inStr.find("//"));
        std::stringstream inStr_ss(inStr);
        
        getline(inStr_ss, inStr, ' ');
        int sign = stoi(inStr);            // 0 for 0, 1 for +, -1 for -
        getline(inStr_ss, inStr, ' ');
        int power = stoi(inStr);           // i for 2^i
        
        // combine
        weightFun.emplace_back( std::make_pair(function, std::make_pair(sign, power)) );
    }
}
