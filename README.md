# RLExact
A program to exactly diagonalize hamiltons for small spin systems.

## Installation
## Linux:

1. Clone this repository with `git clone https://github.com/NBI-Magnetism-Group/RLexact`
2. (Optional) Change flags in `src/RLexact.h` to correspond to the wished mode of operation according to the manual.
3. Run `make -C src` from the main folder to compile the binary.

## Macos tahoe
1. Install brew, such that you can install make and libopenmpi.
2. Run the following commands in your terminal:
   
   ```brew install llvm make```
   
   ```brew install open-mpi```
   
4. Clone this repository with `git clone https://github.com/NBI-Magnetism-Group/RLexact`
5. (Optional) Change flags in `src/RLexact.h` to correspond to the wished mode of operation according to the manual.
3. Run ```make -C src``` from the main folder to compile the binary.
4. Run ```src/RLexact RLexamples/test/test.h$ To execute your first RLexact program```
