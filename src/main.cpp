#include "FASTAParse.h"
#include <iostream>
#include "DNATranslator.h"
#include "SequenceAlignment.h"
#include <pthread.h>

int main(int argc, const char* argv[])
{
    // TODO
        int threads = atoi(argv[3]);
       std::string f1 = argv[1];
       std::string f2 = argv[2];
    // std::string f1 = "input/Zaire_ebolavirus.fasta";
    // std::string f2 = "input/Reston_ebolavirus.fasta";
        SequenceAlignment test(threads,f1,f2);
        test.processGenes();
        test.createFile();
    return 1;

    }

