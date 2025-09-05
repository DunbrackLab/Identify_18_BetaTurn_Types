# Identify_18_BetaTurn_Types
Identifies beta turns in proteins according to Shapovalov, Vucetic, and Dunbrack, PLOSCompBio 2019

Requires mkdssp from DSSP4.5 : https://github.com/PDB-REDO/dssp
It uses the command line version of DSSP4.5 not the python module. The way to install it in /usr/local/bin/mkdssp is:
    git clone https://github.com/PDB-REDO/dssp.git
    cd dssp
    cmake -S . -B build
    cmake --build build
    cmake --install build
