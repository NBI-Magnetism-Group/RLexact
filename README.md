# RLexact
A program to exactly diagonalize hamiltons for small spin systems.

The source code for compiling the binary is found in src. 
The few examples are given in RLexamples.
The manual contains most of the documentation for this project.

## Installation
### Linux:
1. First update apt: ```sudo apt-get update```
2. Then install openmpi and make: ```sudo apt-get install -y make clang pkg-config openmpi-bin libopenmpi-dev```
3. Clone this repository with `git clone https://github.com/NBI-Magnetism-Group/RLexact`
4. (Optional) Change flags in `src/RLexact.h` to correspond to the wished mode of operation according to the manual.
5. Move into the directory of RLexact, i.e Run ```cd RLexact```
6. Run `make -C src` from the main folder to compile the binary.
7. Run ```src/RLexact RLexamples/test/test.h``` To execute your first RLexact program
   
### Macos
1. Install brew, such that you can install make and libopenmpi.
2. Run the following commands in your terminal:
   
   ```brew install llvm make```
   
   ```brew install open-mpi```
   
3. Clone this repository with `git clone https://github.com/NBI-Magnetism-Group/RLexact`
4. (Optional) Change flags in `src/RLexact.h` to correspond to the wished mode of operation according to the manual.
5. Move into the directory of RLexact, i.e Run ```cd RLexact```
6. Run ```make -C src``` from the main folder to compile the binary.
7. Run ```src/RLexact RLexamples/test/test.h``` To execute your first RLexact program

### Windows
1. Install Linux with WSL and then follow Linux instructions.
