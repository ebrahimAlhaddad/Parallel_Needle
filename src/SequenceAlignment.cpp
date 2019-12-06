#include "SequenceAlignment.h"
#include <algorithm>
#include <time.h>


//creates comparison file (.result)
//input:none
//output:none
void SequenceAlignment::createFile(){
    //process indentical pairs between the two sequences
    char* bonds = new char[mResultA.length()];
    for(int i = 0; i < mResultA.length();i++){
        if(mResultA[i] == mResultB[i])
            bonds[i] = '|';
        else
            bonds[i] = ' ';
    }
    std::ofstream oFile("output/match.result");
    if(oFile.is_open()){
        //print headers
        oFile << "A: " << mFastaFile1->mHeader << std::endl;
        oFile << "B: " << mFastaFile2->mHeader << std::endl;
        oFile << "Score: " << mBestscore <<std::endl << std::endl;
        int tempInd;
        int i = 0;
        //for loops are to limit each line to 70 characters
        while(i < mResultB.length()){
            tempInd = i;
            for(int j = 0; j < 70; j++){
                oFile << mResultA[i];
                i++;
                if(i >= mResultA.length())
                    break;
            }
            oFile << std::endl;
            i = tempInd;
            for(int j = 0; j <70; j++){
                oFile << bonds[i];
                i++;
                if(i >= mResultA.length())
                    break;
            }
            oFile << std::endl;
            i = tempInd;
            for(int j = 0; j <70; j++){
                oFile << mResultB[i];
                i++;
                if(i >= mResultA.length())
                    break;
            }
            oFile << std::endl << std::endl;

        }
    }
    oFile.close();
    delete[] bonds;
};

void SequenceAlignment::processGenes(){
       struct timespec start, stop;
    double time;
    if( clock_gettime(CLOCK_REALTIME, &start) == -1) { perror("clock gettime");}
    
    //initialization
    short initScore = 0;
    printf("%d",left);
    //direction initVal = left;
    for(int i = 0; i < mGridWidth + 1; i++){
        mScoregrid[0][i] = initScore;
        initScore--;
        if(i != 0){
            mChargrid[0][i] = left;
        }
    }
    initScore = 0;
    //initVal = above;
    for(int i = 0; i < mGridLength + 1; i++){
        mScoregrid[i][0] = initScore;
        initScore--;
        if(i != 0){
            mChargrid[i][0] = above;
        }
    }
    
    //main loop
    short scoreA = 0;
    short scoreB = 0;
    short scoreC = 0;
    short max;
    for(int i = 1; i < mGridLength + 1; i++){
        for(int j = 1; j < mGridWidth + 1; j++){
            //compute score
            //scoreA
            scoreA = mScoregrid[i-1][j-1];
            if(mFastaFile2->mSequence[i-1] == mFastaFile1->mSequence[j-1]){
                scoreA++;
            } else {
                scoreA--;
            }
            //scoreB
            scoreB = mScoregrid[i][j-1];
            scoreB--;
            //scoreC
            scoreC = mScoregrid[i-1][j];
            scoreC--;
            //compare scores
            if(scoreB > scoreA)
                max = scoreB;
            else
                max = scoreA;
            if(scoreC > max)
                max = scoreC;
            //fill grid
            if(max == scoreA){
                mScoregrid[i][j] = max;
                mChargrid[i][j] = aboveLeft;
            } else if(max == scoreB){
                mScoregrid[i][j] = max;
                mChargrid[i][j] = left;
            } else if(max == scoreC){
                mScoregrid[i][j] = max;
                mChargrid[i][j] = above;
            }
        }
    }

    // for(int i = 0; i < mGridLength; ++i){
    //     for(int j = 0; j < mGridWidth; ++j){
    //         printf("|%d|",mChargrid[i][j]);
    //     }
    //     printf("\n");
    // }
    // std::ofstream oFile("scoreMatrix");
    // for(int i = 0; i <= mGridLength; ++i){
    //     for(int j = 0; j <= mGridWidth; ++j){
    //         oFile << "|" << mScoregrid[i][j] << "|";
    //     }
    //     oFile << std::endl;
    // }
    // oFile.close();

    //picking the perfect alignment
    int finalI = mGridLength;
    int finalJ = mGridWidth;
    printf("\n%d,%d",finalI,finalJ);
    mBestscore = mScoregrid[finalI][finalJ];
    while(finalI != 0 || finalJ != 0){
        if(mChargrid[finalI][finalJ] == aboveLeft){
            mResultB += mFastaFile2->mSequence[finalI - 1];
            mResultA += mFastaFile1->mSequence[finalJ - 1];
            finalI--;
            finalJ--;
        } else if(mChargrid[finalI][finalJ] == left){
            mResultB += '_';
            mResultA += mFastaFile1->mSequence[finalJ - 1];
            finalJ--;
        } else if(mChargrid[finalI][finalJ] == above){
            mResultB += mFastaFile2->mSequence[finalI - 1];
            mResultA += '_';
            finalI--;
        }

    }
    //delete grids to improve performance
    for(int i = 0; i < mGridLength; i++){
        delete[] mScoregrid[i];
    }
    delete[] mScoregrid;
    for(int j = 0; j < mGridLength; j++){
        delete[] mChargrid[j];
    }
    delete[] mChargrid;
    //reverse result sequence
    std::reverse(mResultA.begin(),mResultA.end());
    std::reverse(mResultB.begin(),mResultB.end());
    if( clock_gettime( CLOCK_REALTIME, &stop) == -1 ) { perror("clock gettime");}
    time = (stop.tv_sec - start.tv_sec)+ (double)(stop.tv_nsec - start.tv_nsec)/1e9;
    printf("Execution time = %f sec.\n", time);
};

//initialize grids for the algorithm and parse fasta files
//input: directory into comparison fasta files
SequenceAlignment::SequenceAlignment(std::string &file1, std::string &file2){
    //initialize Fasta parser and DNA Translator
    mFastaFile1 = new FASTAParse(file1);
    mFastaFile2 = new FASTAParse(file2);
    mGridLength = mFastaFile2->mSequence.length();
    mGridWidth = mFastaFile1->mSequence.length();
    if(mGridLength > mGridWidth){
        mResultA.reserve(mGridLength);
        mResultB.reserve(mGridLength);
    } else {
        mResultA.reserve(mGridWidth);
        mResultB.reserve(mGridWidth);
    }
    mScoregrid = new short*[mGridLength + 1];
    mChargrid = new char*[mGridLength + 1];
    for(int i = 0; i < mGridLength + 1; i++){
        mScoregrid[i] = new short[mGridWidth + 1];
        mChargrid[i] = new char[mGridWidth + 1];
    }
    
    
};
