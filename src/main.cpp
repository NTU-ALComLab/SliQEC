#include <boost/program_options.hpp>
#include <sys/time.h> 
#include <fstream>

#include "eqChecker.h"
#include "memMeasure.h"

extern void qasmParser(std::ifstream &inFile, std::vector<GateType> &gates, std::vector<std::vector<int> > &qubits, int &n);
extern void careSetParser(std::ifstream &inFile, std::vector<std::string>& careSet);
extern void weightFunParser(std::ifstream &inFile, std::vector<std::pair<std::vector<std::string>, std::pair<int, int>>>& weightFun);


int main(int argc, char **argv)
{
    // Program options
    namespace po = boost::program_options;
    po::options_description description("Options");
    description.add_options()
    ("help", "produce help message.")
    ("reorder", po::value<bool>()->default_value(true), "allow variable reordering or not.\n"
                                                        "0: disable 1: enable") 
    ("circuit1", po::value<std::string>()->default_value(""), "1st circuit for equivalence checking.")
    ("circuit2", po::value<std::string>()->default_value(""), "2nd circuit for equivalence checking.")
    ("p", po::value<int>()->default_value(0), "conduct full or partial equivalence checking.\n"
                                              "0: full 1: strong partial\n"
                                              "2: unweighted zero-clean-ancilla strong partial\n"
                                              "3: constant-probability weak partial")
    ("nQd",  po::value<int>()->default_value(0), "(only for --p 1/2/3) #data qubits.")
    ("nQkc", po::value<int>()->default_value(0), "(only for --p 1/3)   #clean ancilla qubits.")
    ("nQkd", po::value<int>()->default_value(0), "(only for --p 1/2/3) #dirty ancilla qubits.")
    ("nQm",  po::value<int>()->default_value(0), "(only for --p 1/2/3) #measured qubits.")
    ("nQw",  po::value<int>()->default_value(0), "(only for --p 1/3)   #weighted qubits.")
    ("nQp",  po::value<int>()->default_value(0), "(only for --p 1/2/3) #propagating qubits.")
    ("nQg",  po::value<int>()->default_value(0), "(only for --p 1/2/3) #garbage qubits.")
    ("nQkr", po::value<int>()->default_value(0), "(only for --p 1/3)   #reverted clean ancilla qubits.")
    ("careSet", po::value<std::string>()->default_value(""), "(only for --p 1/2/3) the care set of the circuits.\n"
                                                               "Can be omitted if the care set is the universal set.")
    ("weightFun1", po::value<std::string>()->default_value(""), "(only for --p 1/2/3 and nQw > 0)"
                                                                 "the weight function of the 1st circuit")
    ("weightFun2", po::value<std::string>()->default_value(""), "(only for --p 1/2/3 and nQw > 0)"
                                                                 "the weight function of the 2nd circuit")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, description), vm);
    po::notify(vm);

    if (vm.count("help") || argc == 1)
    {
	    std::cout << description;
	    return 0;
	}


    bool isReorder = vm["reorder"].as<bool>();

    // Parse QASM files
    std::vector<std::vector<GateType>> gates(2);
    std::vector<std::vector<std::vector<int>>> qubits(2);
    int nQ1, nQ2, n; 

    std::ifstream inFile;

    inFile.open(vm["circuit1"].as<std::string>()); 
    if (!inFile)
    {
        std::cerr << "Circuit1 file doesn't exist." << std::endl;
        return 0;
    }
    qasmParser(inFile, gates[0], qubits[0], nQ1);
    inFile.close();

    inFile.open(vm["circuit2"].as<std::string>());
    if (!inFile)
    {
        std::cerr << "Circuit2 file doesn't exist." << std::endl;
        return 0;
    }
    qasmParser(inFile, gates[1], qubits[1], nQ2);
    inFile.close();
   

    // parse parameters
    n = std::max(nQ1, nQ2);

    int p = vm["p"].as<int>();
    
    int nQd, nQkc, nQkd, nQm, nQw, nQp, nQg, nQkr; 
    
    EqType eqType;

    if (p == 0)
    {
        if (nQ1 != nQ2)
        {
            std::cerr << "The two circuits have different number of qubits." << std::endl;
            return 0;
        }

        eqType = EqType::Feq;

        nQd   = -1;
        nQkc  = -1;
        nQkd  = -1;
        nQm   = -1;
        nQw   = -1;
        nQp   = -1;
        nQg   = -1;
        nQkr  = -1;
    }
    
    else
    {
        if      (p == 1)  eqType = EqType::SPeq;
        else if (p == 2)  eqType = EqType::SPeqS;
        else if (p == 3)  eqType = EqType::CWPeq;
        else 
        {
            std::cerr << "--p should be either 0, 1, 2, or 3;" << std::endl;
            return 0;
        }
        
        nQd   = vm["nQd"].as<int>();
        nQkc  = vm["nQkc"].as<int>();
        nQkd  = vm["nQkd"].as<int>();
        nQm   = vm["nQm"].as<int>();
        nQw   = vm["nQw"].as<int>();
        nQp   = vm["nQp"].as<int>();
        nQg   = vm["nQg"].as<int>();
        nQkr  = vm["nQkr"].as<int>();
        
        if (nQd < 0 || nQkc < 0 || nQkd < 0 || nQm < 0 || nQw < 0 || nQp < 0 || nQg < 0 || nQkr < 0)
        {
            std::cerr << "{nQd, nQkc, nQkd, nQm, nQw, nQp, nQg, nQkr} should be non-negative integers." << std::endl;
            return 0;
        }
        if (n != nQd + nQkc + nQkd || n != nQm + nQw + nQp + nQg +  nQkr + nQkd)
        {
            std::cerr << "(n = nQd + nQkc + nQkd = nQm + nQw + nQp + nQg +  nQkr + nQkd) did not hold." << std::endl;
            return 0;
        }
        if (eqType == EqType::SPeqS) 
        {
            if (nQkc != 0 || nQkr != 0 || nQw != 0) 
            {
                std::cerr << "For --p 2, (nQw = nQkc = nQkr) should hold." << std::endl;
                return 0;
            }
        }
        if (nQkr > nQkc)
        {
            std::cerr << "(nQkr <= nQkc) did not hold." << std::endl;
            return 0;
        }
    }

    // parse care set    // can be omitted, default: universal set
    std::vector<std::string> careSet;
    
    inFile.open(vm["careSet"].as<std::string>()); 
    if (inFile)
    {
        careSetParser(inFile, careSet);
        inFile.close();
    }
    
    // parse weight function
    std::vector<std::vector<std::pair<std::vector<std::string>, std::pair<int, int>>>> weightFuns(2);    // weight functions of circuits
    
    if ((p == 0 || p == 1 || p == 2) && nQw > 0) 
    {
        inFile.open(vm["weightFun1"].as<std::string>()); 
        if (!inFile)
        {
            std::cerr << "Weight function 1 file doesn't exist." << std::endl;
            return 0;
        }
        weightFunParser(inFile, weightFuns[0]);
        inFile.close();
        
        inFile.open(vm["weightFun2"].as<std::string>()); 
        if (!inFile)
        {
            std::cerr << "Weight function 2 file doesn't exist." << std::endl;
            return 0;
        }
        weightFunParser(inFile, weightFuns[1]);
        inFile.close();
    }

    // run
    struct timeval tStart, tFinish;
    double elapsedTime;
    double runtime;
    size_t memPeak;

    gettimeofday(&tStart, NULL);

    EquivalenceChecker checker(gates, qubits, n, nQd, nQkc, nQkd, nQm, nQw, nQp, nQg, nQkr, careSet, weightFuns, isReorder, eqType);

    checker.check();

    gettimeofday(&tFinish, NULL);
    elapsedTime = (tFinish.tv_sec - tStart.tv_sec) * 1000.0;
    elapsedTime += (tFinish.tv_usec - tStart.tv_usec) / 1000.0;

    runtime = elapsedTime / 1000.0;
    memPeak = getPeakRSS();
    
    checker.printInfo(runtime, memPeak);

    return 0;
}
