#pragma once

#include <stdio.h>
#include "DNATranslator.h"
#include "FASTAParse.h"
#include <string>



class SequenceAlignment{
private:
//    enum direction:char {
//        aboveLeft,
//        left,
//        above
//    };
    FASTAParse* mFastaFile1;
    FASTAParse* mFastaFile2;
    int mGridLength;
    int mGridWidth;
    short** mScoregrid;
    char** mChargrid;
    std::string mResultA;
    std::string mResultB;
    short mBestscore;
    int mThreads;
public:
    SequenceAlignment(int threads, std::string &file1, std::string &file2);
    void processGenes();
    void createFile();
    
    
};
