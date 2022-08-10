#include <boost/program_options.hpp>
#include "Simulator.h"
#include "util_sim.h"

int main(int argc, char **argv)
{
    namespace po = boost::program_options;
    po::options_description description("Options");
    description.add_options()
    ("help", "produce help message")
    ("print_info", "print statistics such as runtime, memory, etc.")
    ("r", po::value<unsigned int>()->default_value(32), "integer bit size.")
    ("reorder", po::value<bool>()->default_value(1), "allow variable reordering or not.\n"
                                                             "0: disable reordering.\n"
                                                             "1: enable reordering (default option).")
    ("pec", "conduct partial equivalence checking")
    ("circuit1", po::value<std::string>()->implicit_value(""), "1st circuit for equivalence checking")
    ("circuit2", po::value<std::string>()->implicit_value(""), "2nd circuit for equivalence checking")
    ("nQin", po::value<unsigned int>()->default_value(0), "the number of input qubits")
    ("nQout", po::value<unsigned int>()->default_value(0), "the number of output qubits")
    ("s", "use the method for special cases")
    ("output", po::value<std::string>()->implicit_value("exp.csv"), "output experimental data to the file.") // exp
    ("f", "calculate fidelity.")
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

    int inputNumber = vm["nQin"].as<unsigned int>();
    int outputNumber = vm["nQout"].as<unsigned int>();
    assert(inputNumber>0);
    assert(outputNumber>0);
    assert(r>1);

    if (vm.count("output")) // exp
    {
        outFile.open(vm["output"].as<std::string>(), std::ios::app);
        signal(SIGTERM, signalHandler);
    }

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
            if (!inFile){
                std::cout << "no such file : " << vm["circuit1"].as<std::string>() << "\n";
                assert(inFile);
            }
            strStream << inFile.rdbuf(); //read the file
        }
        std::string G = strStream.str(); //str holds the content of the file

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
            if (!inFile){
                std::cout << "no such file : " << vm["circuit2"].as<std::string>() << "\n";
                assert(inFile);
            }
            strStream << inFile.rdbuf(); //read the file
        }
        std::string G_p = strStream.str(); //str holds the content of the file

        // start timer
        gettimeofday(&tStart, NULL);

        PartialEquivalenceChecker checker(G, G_p, inputNumber, outputNumber, r, isReorder);
        if (vm.count("s"))
            checker.runPEC_special();
        else
            checker.runPEC();

        //end timer
        gettimeofday(&tFinish, NULL);
        elapsedTime = (tFinish.tv_sec - tStart.tv_sec) * 1000.0;
        elapsedTime += (tFinish.tv_usec - tStart.tv_usec) / 1000.0;

        runtime = elapsedTime / 1000;
        memPeak = getPeakRSS();
        if (vm.count("print_info"))
            checker.print_info(runtime, memPeak);
    }

    if (vm.count("output")) // exp
    {
        outFile  << "," << runtime << "," << memPeak;
        outFile << std::endl;
        outFile.close();
    }

    return 0;
}
