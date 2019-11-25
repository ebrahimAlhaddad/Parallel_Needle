#include <stdio.h>
#include "DNATranslator.h"


//Initialize state machine array and amino acid map
//Input: A parsed fasta file object pointer
DNATranslator::DNATranslator(FASTAParse *inData){
    /**I = 9
     M = 13
     T = 20
     N = 14
     K = 11
     S = 19
     R = 18
     F = 6
     L = 12
     Y = 25
     C = 3
     W = 23
     P = 16
     Q = 17
     H = 8
     V = 22
     A = 1
     D = 4
     G = 7
     E = 5
     **/
    mFASTAData = inData;
    mAminoCount['A'] = 0;
    mAminoCount['C'] = 0;
    mAminoCount['D'] = 0;
    mAminoCount['E'] = 0;
    mAminoCount['F'] = 0;
    mAminoCount['G'] = 0;
    mAminoCount['H'] = 0;
    mAminoCount['I'] = 0;
    mAminoCount['K'] = 0;
    mAminoCount['L'] = 0;
    mAminoCount['M'] = 0;
    mAminoCount['N'] = 0;
    mAminoCount['P'] = 0;
    mAminoCount['Q'] = 0;
    mAminoCount['R'] = 0;
    mAminoCount['S'] = 0;
    mAminoCount['T'] = 0;
    mAminoCount['V'] = 0;
    mAminoCount['W'] = 0;
    mAminoCount['Y'] = 0;
    mStateGraph[0][0] = 6;
    mStateGraph[0][1] = 11;
    mStateGraph[0][2] = 1;
    mStateGraph[0][3] = 16;
    mStateGraph[1][0] = 2;
    mStateGraph[1][1] = 3;
    mStateGraph[1][2] = 4;
    mStateGraph[1][3] = 5;
    mStateGraph[2][0] = 73;
    mStateGraph[2][1] = 73;
    mStateGraph[2][2] = 73;
    mStateGraph[2][3] = 77; //start codon
    mStateGraph[3][0] = 84;
    mStateGraph[3][1] = 84;
    mStateGraph[3][2] = 84;
    mStateGraph[3][3] = 84;
    mStateGraph[4][0] = 78;
    mStateGraph[4][1] = 78;
    mStateGraph[4][2] = 75;
    mStateGraph[4][3] = 75;
    mStateGraph[5][0] = 83;
    mStateGraph[5][1] = 83;
    mStateGraph[5][2] = 82;
    mStateGraph[5][3] = 82;
    mStateGraph[6][0] = 7;
    mStateGraph[6][1] = 8;
    mStateGraph[6][2] = 9;
    mStateGraph[6][3] = 10;
    mStateGraph[7][0] = 70;
    mStateGraph[7][1] = 70;
    mStateGraph[7][2] = 76;
    mStateGraph[7][3] = 76;
    mStateGraph[8][0] = 83;
    mStateGraph[8][1] = 83;
    mStateGraph[8][2] = 83;
    mStateGraph[8][3] = 83;
    mStateGraph[9][0] = 89;
    mStateGraph[9][1] = 89;
    mStateGraph[9][2] = -1; //stop codon
    mStateGraph[9][3] = -1; //stop codon
    mStateGraph[10][0] = 67;
    mStateGraph[10][1] = 67;
    mStateGraph[10][2] = -1;
    mStateGraph[10][3] = 87;
    mStateGraph[11][0] = 12;
    mStateGraph[11][1] = 13;
    mStateGraph[11][2] = 14;
    mStateGraph[11][3] = 15;
    mStateGraph[12][0] = 76;
    mStateGraph[12][1] = 76;
    mStateGraph[12][2] = 76;
    mStateGraph[12][3] = 76;
    mStateGraph[13][0] = 80;
    mStateGraph[13][1] = 80;
    mStateGraph[13][2] = 80;
    mStateGraph[13][3] = 80;
    mStateGraph[14][0] = 72;
    mStateGraph[14][1] = 72;
    mStateGraph[14][2] = 81;
    mStateGraph[14][3] = 81;
    mStateGraph[15][0] = 82;
    mStateGraph[15][1] = 82;
    mStateGraph[15][2] = 82;
    mStateGraph[15][3] = 82;
    mStateGraph[16][0] = 17;
    mStateGraph[16][1] = 18;
    mStateGraph[16][2] = 19;
    mStateGraph[16][3] = 20;
    mStateGraph[17][0] = 86;
    mStateGraph[17][1] = 86;
    mStateGraph[17][2] = 86;
    mStateGraph[17][3] = 86;
    mStateGraph[18][0] = 65;
    mStateGraph[18][1] = 65;
    mStateGraph[18][2] = 65;
    mStateGraph[18][3] = 65;
    mStateGraph[19][0] = 68;
    mStateGraph[19][1] = 68;
    mStateGraph[19][2] = 69;
    mStateGraph[19][3] = 69;
    mStateGraph[20][0] = 71;
    mStateGraph[20][1] = 71;
    mStateGraph[20][2] = 71;
    mStateGraph[20][3] = 71;
    
};


//prints the content of amino acid count hash map into amino.txt
//input: none
//output: none
void DNATranslator::printAmino(){
    
    std::ofstream oFile("amino.txt");
    if(oFile.is_open()){
        oFile << mFASTAData->mHeader << std::endl;
        oFile << "Total amino acids produced: " << totalAmino << std::endl;
        oFile << "(A) Alanine: " << mAminoCount['A'] << std::endl;
        oFile << "(C) Cysteine: " << mAminoCount['C'] << std::endl;
        oFile << "(D) Aspartic acid: " << mAminoCount['D'] << std::endl;
        oFile << "(E) Glutamic acid: " << mAminoCount['E'] << std::endl;
        oFile << "(F) Phenylalanine: " << mAminoCount['F'] << std::endl;
        oFile << "(G) Glycine: " << mAminoCount['G'] << std::endl;
        oFile << "(H) Histidine: " << mAminoCount['H'] << std::endl;
        oFile << "(I) Isoleucine: " << mAminoCount['I'] << std::endl;
        oFile << "(K) Lysine: " << mAminoCount['K'] << std::endl;
        oFile << "(L) Leucine: " << mAminoCount['L'] << std::endl;
        oFile << "(M) Methionine: " << mAminoCount['M'] << std::endl;
        oFile << "(N) Asparagine: " << mAminoCount['N'] << std::endl;
        oFile << "(P) Proline: " << mAminoCount['P'] << std::endl;
        oFile << "(Q) Glutamine: " << mAminoCount['Q'] << std::endl;
        oFile << "(R) Arginine: " << mAminoCount['R'] << std::endl;
        oFile << "(S) Serine: " << mAminoCount['S'] << std::endl;
        oFile << "(T) Threonine: " << mAminoCount['T'] << std::endl;
        oFile << "(V) Valine: " << mAminoCount['V'] << std::endl;
        oFile << "(W) Tryptophan: " << mAminoCount['W'] << std::endl;
        oFile << "(Y) Tyrosine: " << mAminoCount['Y'] << std::endl;
    }
    oFile.close();
};


//converts A, T, C, G into index value
//input: character from sequence
//output: index for the state machine
int elemToInd(char element){
    if(element == 'T')
        return 0;
    else if(element == 'C')
        return 1;
    else if (element == 'A')
        return 2;
    else
        return 3;
};

//counts amino acids produced by DNA sequence
//input: none
//output: none
void DNATranslator::translateDNA(){
    /**
     
     //str | T | C | A | G |
     //0   | 6 | 11 | 1 | 16 |  -> looking for ATG
     //1   | 2 | 3 | 4 | 5 |  -> A**
     //2   | ‘I’| ‘I’| ‘I’|’M’/start| -> AT*
     //3   |’T’|’T’|’T’|’T’| -> AC*
     //4   |’N’|’N’|’K’|’K’| ->AA*
     //5   |’S’|’S’|’R’|’R’| -> AG*
     
     //6   | 7 | 8 | 9 | 10 | -> T**
     //7   |’F’|’F’|’L’|’L’| -> TT*
     //8   |’S’|’S’|’S’|’S’| -> TC*
     //9   |’Y’|’Y’|sto|sto| -> TA*
     //10   |’C’|’C’|sto|’W’| -> TG*
     
     //11   |12|13|14|15| -> C**
     //12   |’L’|’L’| ’L’|’L’| -> CT*
     //13   |’P’|’P’|’P’|’P’| -> CC*
     //14   |’H’|’H’|’Q’|’Q’| -> CA*
     //15   |’R’|’R’|’R’|’R’| -> CG*
     
     //16  |17|18|19|20| ->G**
     //17  | V | V | V | V | ->GT*
     //18  | A | A | A | A | ->GC*
     //19  | D | D | E | E | ->GA*
     //20  |G |G |G | G| ->GG*
     **/
    
    int i = 0;
    //read first element
    char elem = mFASTAData->mSequence[i];
    int stateInd = 0;
    int inStateInd;
    char aminoProduced;
    bool bufferSeq = true;
    
    int startState = 0;
    while(elem != '\0'){
        if(bufferSeq){
            //check starter codon
            switch (startState) {
                case 0:
                    if(elem == 'A')
                        startState = 1;
                    break;
                    
                case 1:
                    if(elem == 'T')
                        startState = 2;
                    else if(elem != 'A')
                        startState = 0;
                    break;
                    
                case 2:
                    if(elem == 'G'){
                        bufferSeq = false;
                        mAminoCount['M']++;
                        startState = 0;
                    } else if(elem == 'A'){
                        startState = 1;
                    } else
                        startState = 0;
                    break;
            }
            
        } else {
            inStateInd = elemToInd(elem);
            stateInd = mStateGraph[stateInd][inStateInd];
            //produce amino acid and count it if the array element is an ASCII value
            if(stateInd > 20){
                aminoProduced = char(stateInd);
                mAminoCount[aminoProduced]++;
                stateInd = 0;
            } else if(stateInd == -1){
                //stop producing amino acids and start buffering characters
                bufferSeq = true;
                stateInd = 0;
            }
            
        }
        i++;
        elem = mFASTAData->mSequence[i];
    }
    
    //count total amino acids
    totalAmino = 0;
    for(auto elem : mAminoCount)
    {
        totalAmino += elem.second;
    }
    
    //print amino map
    printAmino();
    
};
