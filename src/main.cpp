#include <boost/program_options.hpp>
#include <sys/time.h> 
#include <fstream>

#include "eqChecker.h"
#include "memMeasure.h"

extern void qasmParser(std::ifstream &inFile, std::vector<GateType> &gates, std::vector<std::vector<int> > &qubits, int &n);

int main(int argc, char **argv)
{
    // Program options
    namespace po = boost::program_options;
    po::options_description description("Options");
    description.add_options()
    ("help", "produce help message.")
    ("reorder", po::value<bool>()->default_value(true), "allow variable reordering or not.\n"
                                                             "0: disable 1: enable") 
    ("circuit1", po::value<std::string>()->implicit_value(""), "1st circuit for equivalence checking.")
    ("circuit2", po::value<std::string>()->implicit_value(""), "2nd circuit for equivalence checking.")
    ("p", po::value<bool>()->default_value(false), "conduct full or partial equivalence checking.\n"
                                                    "0: full 1: partial")
    ("nQd", po::value<int>()->default_value(0), "(only for --p 1) the number of data qubits.")
    ("nQm", po::value<int>()->default_value(0), "(only for --p 1) the number of measured qubits.")
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

    n = std::max(nQ1, nQ2);

    bool p = vm["p"].as<bool>();
    
    int nQd, nQm; 
    
    EqType eqType;

    if(!p)
    {
        if(nQ1 != nQ2)
        {
            std::cerr << "The two circuits have different number of qubits." << std::endl;
            return 0;
        }

        eqType = EqType::Feq;

        nQd = -1;
        nQm = -1;
    }
    else
    {
        nQd = vm["nQd"].as<int>();
        
        if(nQd < 0)
        {
            std::cerr << "nQd should be a positive integer." << std::endl;
            return 0;
        }
        else if(nQd > nQ1 || nQd > nQ2)
        {
            std::cerr << "nQd cannot be larger than #qubit in circuit1 and circuit2." << std::endl;
            return 0;
        }

        nQm = vm["nQm"].as<int>();

        if(nQm < 0)
        {
            std::cerr << "nQm should be a positive integer." << std::endl;
            return 0;
        }
        else if(nQm > nQ1 || nQm > nQ2)
        {
            std::cerr << "nQm cannot be larger than #qubit in circuit1 and circuit2." << std::endl;
            return 0;
        }

        if(nQd == n) 
            eqType = EqType::PeqS;
        else
            eqType = EqType::Peq;
    }

    struct timeval tStart, tFinish;
    double elapsedTime;
    double runtime;
    size_t memPeak;

    gettimeofday(&tStart, NULL);

    EquivalenceChecker checker(gates, qubits, n, nQd, nQm, isReorder, eqType);

    checker.check();

    gettimeofday(&tFinish, NULL);
    elapsedTime = (tFinish.tv_sec - tStart.tv_sec) * 1000.0;
    elapsedTime += (tFinish.tv_usec - tStart.tv_usec) / 1000.0;

    runtime = elapsedTime / 1000.0;
    memPeak = getPeakRSS();
    
    checker.printInfo(runtime, memPeak);

    return 0;
}
