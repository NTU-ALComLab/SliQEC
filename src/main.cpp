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
    ("print_info", "print statistics such as runtime, memory, etc.")
    ("r", po::value<int>()->default_value(32), "integer bit size.")
    ("reorder", po::value<bool>()->default_value(true), "allow variable reordering or not.\n"
                                                             "0: disable reordering.\n"
                                                             "1: enable reordering (default option).")
    ("p", "toggle conducting partial equivalence checking.")
    ("circuit1", po::value<std::string>()->implicit_value(""), "1st circuit for equivalence checking.")
    ("circuit2", po::value<std::string>()->implicit_value(""), "2nd circuit for equivalence checking.")
    ("nQin", po::value<int>()->default_value(0), "the number of input qubits.")
    ("nQout", po::value<int>()->default_value(0), "the number of output qubits.")
    ("s", "toggle using algorithm2 for partial equivalence checking.")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, description), vm);
    po::notify(vm);

    if (vm.count("help") || argc == 1)
    {
	    std::cout << description << std::endl;
	    return 1;
	}


    bool isReorder = vm["reorder"].as<bool>();

    int r = vm["r"].as<int>();
    if(r <= 0)
    {
        std::cerr << "r should be a positive integer." << std::endl;
        return 0;
    }

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

    int nQin, nQout; 
    

    EqType eqType;

    if(!vm.count("p"))
    {
        if(nQ1 != nQ2)
        {
            std::cerr << "The two circuits have different number of qubits." << std::endl;
            return 0;
        }

        eqType = EqType::Feq;

        nQin = -1;
        nQout = -1;
    }
    else
    {
        nQin = vm["nQin"].as<int>();
        
        if(nQin < 0)
        {
            std::cerr << "nQin should be a positive integer." << std::endl;
            return 0;
        }
        else if(nQin > nQ1 || nQin > nQ2)
        {
            std::cerr << "nQin cannot be larger than #qubit in circuit1 and circuit2." << std::endl;
            return 0;
        }

        nQout = vm["nQout"].as<int>();

        if(nQout < 0)
        {
            std::cerr << "nQout should be a positive integer." << std::endl;
            return 0;
        }
        else if(nQout > nQ1 || nQout > nQ2)
        {
            std::cerr << "nQout cannot be larger than #qubit in circuit1 and circuit2." << std::endl;
            return 0;
        }

        if(nQin == n && !vm.count("s"))
            eqType = EqType::PeqS;
        else
            eqType = EqType::Peq;
    }

    struct timeval tStart, tFinish;
    double elapsedTime;
    double runtime;
    size_t memPeak;

    gettimeofday(&tStart, NULL);

    EquivalenceChecker checker(gates, qubits, n, nQin, nQout, r, isReorder, eqType);

    checker.check();

    gettimeofday(&tFinish, NULL);
    elapsedTime = (tFinish.tv_sec - tStart.tv_sec) * 1000.0;
    elapsedTime += (tFinish.tv_usec - tStart.tv_usec) / 1000.0;

    runtime = elapsedTime / 1000.0;
    memPeak = getPeakRSS();
    if (vm.count("print_info"))
        checker.printInfo(runtime, memPeak);

    return 0;
}
