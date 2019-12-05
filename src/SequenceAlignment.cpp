#include "SequenceAlignment.h"
#include <algorithm>
#include <cmath>
#include <time.h>
#include <pthread.h>

enum direction:char
{
    aboveLeft,
    left,
    above
};

// thread data for each thread when initializing matrix
pthread_mutex_t lock_Init;
struct  thread_data_Init
{
    int thread_id_Init;
    short** pScoreGrid;
    char** pCharGrid;
    int pGridLength;
    int pGridWidth;
    int pInitScore;


};

void* SequenceAlignment::matrixInit(void *threadArg)
{
  enum direction:char
  {
      aboveLeft,
      left,
      above
  };
  struct thread_data_Init *my_data_Init;
  my_data_Init = (struct thread_data_Init *) threadArg;
  short** myScoreGrid_Init = my_data_Init->pScoreGrid;
  char** myCharGrid_Init = my_data_Init->pCharGrid;
  int myThreadId_Init = my_data_Init -> thread_id_Init;
  int myGridLength_Init = my_data_Init -> pGridLength;
  int myGridWith_Init = my_data_Init -> pGridWidth;
  int myScore_Init = my_data_Init -> pInitScore;
  // char myLeft = my_data_Init -> DirectionLeft;
  // char myAbove = my_data_Init -> DirectionAbove;
  if(myThreadId_Init == 0)
  {
    for(int i = 0; i < myGridWith_Init + 1; i++)
    {
        myScoreGrid_Init[0][i] = myScore_Init;
        myScore_Init--;
        if(i != 0)
        {
            myCharGrid_Init[0][i] = left;

        }
    }
  }
  else
  {
    //initVal = above;
    for(int i = 0; i < myGridLength_Init + 1; i++)
    {
        myScoreGrid_Init[i][0] = myScore_Init;
        myScore_Init--;
        if(i != 0)
        {
            myCharGrid_Init[i][0] = above;
        }
    }
  }

  pthread_exit(NULL);
}



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

void SequenceAlignment::processGenes()
{
  struct timespec start, stop;
   double time;

//measure time starts now
   if( clock_gettime(CLOCK_REALTIME, &start) == -1) { perror("clock gettime");}
  //**********begin of paralell section below ************
  //return val when creating thread
    int return_Val_Init;
    //array with all thread data needed to be passed when creating active thread
    struct  thread_data_Init  thread_data_arrayJP[2];
    //initialization
    short initScore = 0;
    //declaring threads:
   pthread_t threads_Init[2];
   for(int i = 0; i < 2; i++)
   {
     thread_data_arrayJP[i].thread_id_Init = i;
     thread_data_arrayJP[i].pScoreGrid = mScoregrid;
     thread_data_arrayJP[i].pCharGrid = mChargrid;
     thread_data_arrayJP[i].pGridLength = mGridLength;
     thread_data_arrayJP[i].pGridWidth = mGridWidth;
     thread_data_arrayJP[i].pInitScore = initScore;
     thread_data_arrayJP[i].pInitScore = initScore;
     // thread_data_arrayJP[i].DirectionLeft = left;
     // thread_data_arrayJP[i].DirectionAbove= above;
     return_Val_Init = pthread_create(&threads_Init[i],NULL,&SequenceAlignment::matrixInit, (void *) &thread_data_arrayJP[i]);
       //in case pthreate create fails and returns 0
       if (return_Val_Init) { printf("ERROR; return code from pthread_create() is %d\n", return_Val_Init); exit(-1);}
   }

   for(int i=0; i < 2; i++)
   {
     (void) pthread_join(threads_Init[i], NULL);
   }

   // for(int i = 0; i < mGridWidth + 1; i++)
   // {
   //    std::cout << "score Grid[0][" << i << "]: " << mScoregrid[0][i];
   //    std::cout << "  char Grid[0][" << i << "]: ";
   //    if(mChargrid[0][i] == left)
   //    {
   //      std::cout << "== left \n";
   //    }
   //    else
   //    {
   //      std::cout << "\n";
   //    }
   // }
   // for(int i = 0; i < mGridLength + 1; i++)
   // {
   //    std::cout << "score Grid[" << i << "][0]: " << mScoregrid[i][0];
   //    std::cout << "  char Grid[" << i << "][0]: ";
   //    if(mChargrid[i][0] == above)
   //    {
   //      std::cout << "== above \n";
   //    }
   //    else
   //    {
   //      std::cout << "\n";
   //    }
   // }



      //**********end of paralell section below ************


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
    //picking the perfect alignment
    int finalI = mGridLength;
    int finalJ = mGridWidth;
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

    // measure the end time here
      if( clock_gettime( CLOCK_REALTIME, &stop) == -1 ) { perror("clock gettime");}
      time = (stop.tv_sec - start.tv_sec)+ (double)(stop.tv_nsec - start.tv_nsec)/1e9;

      // print out the execution time here
  printf("Execution time = %f sec\n", time);

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

    //**********begin of paralell section below ************
    for(int i = 0; i < mGridLength + 1; i++){
        mScoregrid[i] = new short[mGridWidth + 1];
        mChargrid[i] = new char[mGridWidth + 1];
    }
        //**********end of paralell section below ************


};
