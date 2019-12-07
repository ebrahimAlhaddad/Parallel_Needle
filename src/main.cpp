#include "FASTAParse.h"
#include <iostream>
#include "DNATranslator.h"
#include "SequenceAlignment.h"
int main(int argc, const char* argv[])
{
    // TODO
        std::string f1 = argv[1];
        std::string f2 = argv[2];
        int threads = atoi(argv[3]);
        SequenceAlignment test(f1,f2,threads);
        test.processGenes();
        test.createFile();
    return 1;

    }

