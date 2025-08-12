# kP-algorithm-for-atomic-pattern

This project is an implementation of the Left-to-Right (L2R), Left-to-Right with projective coordinate randomazation (L2RwithPCR) and the Right-to-Left (R2L) scalar multiplication algorithm based on the Chevallier-Mames et al.'s atomic pattern, targeting Elliptic Curve Cryptosystems (ECC). The implementation is done for a Texas Instruments microcontroller environment, using Code Composer Studio (CCS).

The aim is to investigate the distinguishability of Chevallier-Mames et al.'s atomic blocks on embedded systems.

#Main Files: R2L_main.c, L2R_main.c and L2RwithPRC_main.c

These files contains the core implementation of the the scalar multiplication algorithms using Chevallier-Mames et al.'s atomic blocks. The algorithm is written in C and implemented using the FLECC_IN_C cryptographic library, which provides essential field and ECC operations.

If you’re looking for the core logic of the algorithm, start here.
