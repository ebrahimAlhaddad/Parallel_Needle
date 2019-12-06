#include "SequenceAlignment.h"
#include <algorithm>
#include <pthread.h>
#include <time.h>
#define aboveLeft 0
#define left 1
#define above 2

//length is ROWS
//width is mGridWidthS
// A utility function to find min of two integers
int minu(int a, int b)
{ return (a < b)? a: b; } 
// A utility function to find min of three integers
int min(int a, int b, int c)
{ return minu(minu(a, b), c);}
// A utility function to find max of two integers
int max(int a, int b)
{ return (a > b)? a: b; }



struct thread_data{
    int thread_id;
    short **p_scoreGrid;
    char **p_charGrid;
    int p_gridLength;
    int p_gridWidth;
    FASTAParse *p_fastaFile1;
    FASTAParse *p_fastaFile2;
    int *p_i_indices;
    int *p_j_indices;
    int start;
    int end;
};

void *threadFunc(void *threadarg){
    struct thread_data *my_data;
    my_data = (struct thread_data *) threadarg;
    int thread_id = my_data->thread_id;
    short **mScoregrid = my_data->p_scoreGrid;
    char **mChargrid = my_data->p_charGrid;
    int mGridLength = my_data->p_gridLength;
    int mGridWidth = my_data->p_gridWidth;
    FASTAParse *mFastaFile1 = my_data->p_fastaFile1;
    FASTAParse *mFastaFile2 = my_data->p_fastaFile2;
    int *i_indices = my_data->p_i_indices;
    int *j_indices = my_data->p_j_indices;
    int start = my_data->start;
    int end = my_data->end;
    //main loop
    short scoreA = 0;
    short scoreB = 0;
    short scoreC = 0;
    short max;
    //printf("|start:%d, end:%d|",start,end);
    for(int k = start; k < end; k++){
        //compute score
            //scoreA
            //printf("T(%d,%d)",i_indices[k],j_indices[k]);
            
            scoreA = mScoregrid[i_indices[k]-1][j_indices[k]-1];
            if(mFastaFile2->mSequence[i_indices[k]-1] == mFastaFile1->mSequence[j_indices[k]-1]){
                scoreA++;
            } else {
                scoreA--;
            }
            //scoreB
            scoreB = mScoregrid[i_indices[k]][j_indices[k]-1];
            scoreB--;
            //scoreC
            scoreC = mScoregrid[i_indices[k]-1][j_indices[k]];
            scoreC--;
            //compare scores
            if(scoreB > scoreA)
                max = scoreB;
            else
                max = scoreA;
            if(scoreC > max)
                max = scoreC;
            //fill grid
                        //printf("T(%d) ",max);

            if(max == scoreA){
                mScoregrid[i_indices[k]][j_indices[k]] = max;
                mChargrid[i_indices[k]][j_indices[k]] = aboveLeft;
            } else if(max == scoreB){
                mScoregrid[i_indices[k]][j_indices[k]] = max;
                mChargrid[i_indices[k]][j_indices[k]] = left;
            } else if(max == scoreC){
                mScoregrid[i_indices[k]][j_indices[k]] = max;
                mChargrid[i_indices[k]][j_indices[k]] = above;
            }
    }
    
    
    //printf("Thread %d completed execution \n", thread_id);
    pthread_exit(NULL);
}
//creates comparison file (.result)
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

    //Diagonal Traversal
    const int indices_size = minu(mGridLength,mGridWidth) + 1;
    //printf("length:%d, width:%d:",mGridLength,mGridWidth);
    int *i_indices = new int[indices_size];
    int *j_indices = new int[indices_size];
    //pThread setup
    int rc;
    int NUM_THREADS = mThreads;
    struct thread_data thread_data_array[NUM_THREADS];
    pthread_t threads[NUM_THREADS];
        for(int i = 1; i < mGridWidth; i++){
            
            int r = 0;
            int count = 0;
            for(int j = i; j >= 1; --j){
                r = r + 1;
                int c = j;
                if(r > mGridLength){
                    break;
                }else{
                    i_indices[count] = r;
                    j_indices[count] = c;
                    count++;
                    //printf("(%d,%d)",r,c);
                    }
            }
            
        // for(int k = 0; k < count; k++){
        //     printf("(%d,%d) ", i_indices[k],j_indices[k]);
        // }
        // printf("\n");
        for(int i = 0; i < NUM_THREADS; i++){
                            thread_data_array[i].thread_id = i;
                            thread_data_array[i].p_charGrid = mChargrid;
                            thread_data_array[i].p_scoreGrid = mScoregrid;
                            thread_data_array[i].p_gridWidth = mGridWidth;
                            thread_data_array[i].p_gridLength = mGridLength;
                            thread_data_array[i].p_fastaFile1 = mFastaFile1;
                            thread_data_array[i].p_fastaFile2 = mFastaFile2;
                            thread_data_array[i].p_i_indices = i_indices;
                            thread_data_array[i].p_j_indices = j_indices;
                            thread_data_array[i].start = ((count)/NUM_THREADS)*i;
                            if(i == NUM_THREADS -1 ){
                                thread_data_array[i].end = count;
                            } else {
                                thread_data_array[i].end = ((count)/NUM_THREADS)*(i+1);
                            }
                            //printf("(start:%d,end:%d)",thread_data_array[i].start,thread_data_array[i].end);
                            rc = pthread_create(&threads[i],NULL,threadFunc,(void*)&thread_data_array[i]);
                            if (rc) { printf("ERROR; return code from pthread_create() is %d\n", rc); exit(-1);}
                        }

                        for(int i=0; i <NUM_THREADS; i++){
                        (void) pthread_join(threads[i], NULL);
                        } 
        }

    for(int k = 1; k <= mGridLength; k++){

        int j = mGridWidth;
        int count = 0;
        for(int i = k; i <= mGridLength; ++i){
            int r = i;
            int c = j;
            --j;
            i_indices[count] = r;
            j_indices[count] = c;
            ++count;
        }

        for(int i = 0; i < NUM_THREADS; i++){
                            thread_data_array[i].thread_id = i;
                            thread_data_array[i].p_charGrid = mChargrid;
                            thread_data_array[i].p_scoreGrid = mScoregrid;
                            thread_data_array[i].p_gridWidth = mGridWidth;
                            thread_data_array[i].p_gridLength = mGridLength;
                            thread_data_array[i].p_fastaFile1 = mFastaFile1;
                            thread_data_array[i].p_fastaFile2 = mFastaFile2;
                            thread_data_array[i].p_i_indices = i_indices;
                            thread_data_array[i].p_j_indices = j_indices;
                            thread_data_array[i].start = ((count)/NUM_THREADS)*i;
                            if(i == NUM_THREADS -1 ){
                                thread_data_array[i].end = count;
                            } else {
                                thread_data_array[i].end = ((count)/NUM_THREADS)*(i+1);
                            }
                            //printf("(start:%d,end:%d)",thread_data_array[i].start,thread_data_array[i].end);
                            rc = pthread_create(&threads[i],NULL,threadFunc,(void*)&thread_data_array[i]);
                            if (rc) { printf("ERROR; return code from pthread_create() is %d\n", rc); exit(-1);}
                        }
                        for(int i=0; i <NUM_THREADS; i++){
                        (void) pthread_join(threads[i], NULL);
                        } 
                                                //printf("\n");

    }
            
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
SequenceAlignment::SequenceAlignment(int threads, std::string &file1, std::string &file2){
    //initialize Fasta parser and DNA Translator
    mThreads = threads;
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