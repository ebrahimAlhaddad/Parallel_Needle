#include <iostream>

using namespace std;
int minu(int a, int b)
{ return (a < b)? a: b; } 
// A utility function to find min of three integers
int min(int a, int b, int c)
{ return minu(minu(a, b), c);}
// A utility function to find max of two integers
int max(int a, int b)
{ return (a > b)? a: b; }
int main()
{
    // int mGridLength = 4;
    // int mGridWidth =  6;
    // for (int line=2; line<=(mGridLength + mGridWidth -1); line++)
    // {
    //     int start_col =  max(0, line-mGridLength);
    //     int count = min(line, (mGridLength-start_col), mGridLength);
    //     for (int j=1; j<count; j++)
    //         {
    //          printf("(%d,%d)",minu(mGridLength, line)-j,start_col+j);
    //         }
    //         printf("start_col:%d, count:%d line:%d\n",start_col,count,line);
    // }

    int ROW = 7;
    int COL = 11;
    for(int i = 0; i < ROW; ++i){
        for(int j = 0; j < COL; ++j){
            printf("|(%d,%d)|",i,j);
        }
        printf("\n");
    }

    for(int i = 1; i < COL; i++){
        int r = 0;
        for(int j = i; j >= 1; --j){
            r = r + 1;
            int c = j;
            if(r > ROW-1){
                break;
            }else{
                printf("(%d,%d)",r,c);
            }
        }
        printf("\n");
    }

    for(int k = 2; k < ROW; k++){
        int j = COL-1;
        for(int i = k; i < ROW; ++i){
            int r = i;
            int c = j;
            --j;
                printf("(%d,%d)",r,c);
        }
      printf("\n");

    }
}


