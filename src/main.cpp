#include <boost/program_options.hpp>

#include "peqChecker.h"
#include "memMeasure.h"

int main(int argc, char **argv)
{
    namespace po = boost::program_options;
    po::options_description description("Options");
    description.add_options()
    ("help", "produce help message.")
    ("print_info", "print statistics such as runtime, memory, etc.")
    ("r", po::value<unsigned int>()->default_value(32), "integer bit size.")
    ("reorder", po::value<bool>()->default_value(1), "allow variable reordering or not.\n"
                                                             "0: disable reordering.\n"
                                                             "1: enable reordering (default option).")
    ("pec", "conduct partial equivalence checking")
    ("circuit1", po::value<std::string>()->implicit_value(""), "1st circuit for equivalence checking.")
    ("circuit2", po::value<std::string>()->implicit_value(""), "2nd circuit for equivalence checking.")
    ("nQin", po::value<unsigned int>()->default_value(0), "the number of input qubits.")
    ("nQout", po::value<unsigned int>()->default_value(0), "the number of output qubits.")
    ("s", "use the method for special cases")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, description), vm);
    po::notify(vm);

    if (vm.count("help") || argc == 1)
    {
	    std::cout << description << std::endl;
	    return 1;
	}

    struct timeval tStart, tFinish;
    double elapsedTime;
    double runtime;
    size_t memPeak;

    int r = vm["r"].as<unsigned int>();
    bool isReorder = vm["reorder"].as<bool>();

    int nQubitIn = vm["nQin"].as<unsigned int>();
    int nQubitOut = vm["nQout"].as<unsigned int>();
    assert(nQubitIn > 0);
    assert(nQubitOut > 0);
    assert(r > 1);

    if (vm.count("pec"))
    {
        // read in file into a string
        std::stringstream strStream;
        if (vm["circuit1"].as<std::string>() == "")
        {
            strStream << std::cin.rdbuf();
        }
        else
        {
            std::ifstream inFile;
            inFile.open(vm["circuit1"].as<std::string>()); //open the input file
            if (!inFile)
            {
                std::cerr << "No such file : " << vm["circuit1"].as<std::string>() << "\n";
                return 0;
            }
            strStream << inFile.rdbuf(); //read the file
        }
        std::string circuit1Str = strStream.str(); //str holds the content of the file

        strStream.str("");
        strStream.clear();
        if (vm["circuit2"].as<std::string>() == "")
        {
            strStream << std::cin.rdbuf();
        }
        else
        {
            std::ifstream inFile;
            inFile.open(vm["circuit2"].as<std::string>()); //open the input file
            if (!inFile)
            {
                std::cerr << "No such file : " << vm["circuit2"].as<std::string>() << "\n";
                return 0;
            }
            strStream << inFile.rdbuf(); //read the file
        }
        std::string circuit2Str = strStream.str(); //str holds the content of the file

        // start timer
        gettimeofday(&tStart, NULL);

        PartialEquivalenceChecker checker(circuit1Str, circuit2Str, nQubitIn, nQubitOut, r, isReorder);

        if (vm.count("s"))
            checker.runPECSpecial();
        else
            checker.runPEC();

        //end timer
        gettimeofday(&tFinish, NULL);
        elapsedTime = (tFinish.tv_sec - tStart.tv_sec) * 1000.0;
        elapsedTime += (tFinish.tv_usec - tStart.tv_usec) / 1000.0;

        runtime = elapsedTime / 1000.0;
        memPeak = getPeakRSS();
        if (vm.count("print_info"))
            checker.printInfo(runtime, memPeak);
    }

    return 0;
}
