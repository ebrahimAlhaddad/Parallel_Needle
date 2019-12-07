#pragma once

#include <stdio.h>
#include "DNATranslator.h"
#include "FASTAParse.h"
#include <string>



class SequenceAlignment{
private:
    int mThreads;
    FASTAParse* mFastaFile1;
    FASTAParse* mFastaFile2;
    int mGridLength;
    int mGridWidth;
    short** mScoregrid;
    char** mChargrid;
    std::string mResultA;
    std::string mResultB;
    short mBestscore;
public:
    static void* matrixInit(void *threadArg);
    SequenceAlignment(std::string &file1, std::string &file2,int threads);
    void processGenes();
    void createFile();


};
