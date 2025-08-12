Implementation of Chevallier-Mames et al.'s (MANA) atomic pattern kP-algorithm.

This project is an implementation of the Left-to-Right (L2R), Left-to-Right with projective coordinate randomization (L2RwithPCR) & Right-to-Left (R2L) scalar multiplication algorithm based on the Chevallier-Mames et al.'s atomic pattern, targeting Elliptic Curve Cryptosystems (ECC).
The implementation is done for a Texas Instruments microcontroller environment, using Code Composer Studio (CCS).

The aim is to investigate the distinguishability of Chevallier-Mames et al.'s atomic blocks on embedded systems.
Main Files:  L2R_main.c, L2RwithPCR_main.c, and R2L_main.c
Thess file contains the core implementation of the the Left-to-Right (L2R), Left-to-Right with projective coordinate randomization (L2RwithPCR) & Right-to-Left (R2L) scalar multiplication algorithm using Chevallier-Mames et al.'s atomic blocks. 
The algorithm is written in C and implemented using the FLECC_IN_C cryptographic library, which provides essential field and ECC operations.
