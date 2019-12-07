# Diagonal Traversal Parallel Computation

### Steps:
- Run make
- To execute program in commandline :
    - `./test input/small_test1.fasta input/small_test2.fasta [#ofthreads]`
    - `./test input/Zaire_ebolavirus.fasta input/Reston_ebolavirus.fasta [#ofthreads]`
    - `./test input/TAS2R16_Homo_sapiens.fasta input/TAS2R16_Pan_troglodytes.fasta [#ofthreads]`
    - `./test input/A1C1_Homo_sapiens.fasta input/AKAP1_Homo_sapiens.fasta [#ofthreads]`

### Output is in the output folder as 'match.result' compare it with pre-checked sampleoutput
Note: input files must be called in the order mentioned. This is to make sure the matrix has larger width than length. Necessary for diagonal traversal without running into issues