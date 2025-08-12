#include "F28x_Project.h"
#include "driverlib.h"
#include "device.h"


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>


#include <flecc_in_c/gfp/gfp_const_runtime.h>
#include <flecc_in_c/gfp/gfp.h>
#include <flecc_in_c/io/io.h>
#include <flecc_in_c/bi/bi.h>
#include <flecc_in_c/utils/param.h>
#include <flecc_in_c/utils/parse.h>
#include <time.h>

// Curve parameters (P-256)
const char *curve_type_str = "secp256r1";
// R squared
const char *Rsq_str = "4FFFFFFFDFFFFFFFFFFFFFFFEFFFFFFFBFFFFFFFF0000000000000003";

// Coordinates of P
const char *p0x_str = "6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296";
const char *p0y_str = "4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5";
const char *p0z_str = "0000000000000000000000000000000000000000000000000000000000000001";
const char *q0T0 =  "0xffffffff00000001000000000000000000000000fffffffffffffffffffffffc";


// Coordinates of Q
const char *q0x_str = "6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296";
const char *q0y_str = "4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5";
const char *q0z_str = "0000000000000000000000000000000000000000000000000000000000000001";



// scalar value
//const char *KB =  "1001101101011111110111"; //26D7F7 this should change according to your data
const char *KB =  "11111"; //31 this should change according to your data
//const char *KB =  "100111111100011010000101110001011111110000110100111111110011011"; //4FE342E2FE1A7F9B this Long Key used to see why we had to use 5-bit key

// Field elements &
gfp_t T0, T1, T2, T3, T4, T5, T6, T7, T8, T9;
gfp_t Q0x, Q0y, Q0z, Px, Py, Pz, r_sq;
gfp_t DummyOp;

//Affinate Coordiantes
gfp_t X_A, Y_A;


eccp_parameters_t curve_params;
gfp_prime_data_t prime_data;

int counter = 0;

/**void perform_dummy_nops(int count) {
    for (int i = 0; i < count; i++) {

    counter++;

    // DELAY_US(10);
    }

}*/

void perform_dummy_nops(int microseconds) {
    DEVICE_DELAY_US(microseconds);
}


// Parse a bigint string into a gfp_t type
void parse_bigint(const char *string, uint_t *big_int, const int bi_length) {
    int len = strlen(string);
    bigint_parse_hex_var(big_int, bi_length, string, len);
}

// Modular multiplication function for M-A-N-A
void modular_multiply(gfp_t res,const gfp_t a, const gfp_t b, const gfp_prime_data_t *prime_data) {
    gfp_cr_mont_multiply_sos(res, a, b, prime_data); // Perform Montgomery multiplication
    gfp_cr_mont_multiply_sos(res, res, r_sq, prime_data); // Adjust result using R-squared
}

#define EXT_PIN 32  // Example: GPIO32, from J1

void init_ext_pin() {
    GPIO_setDirectionMode(EXT_PIN, GPIO_DIR_MODE_OUT);
    GPIO_setPadConfig(EXT_PIN, GPIO_PIN_TYPE_STD);
    GPIO_writePin(EXT_PIN, 0);  // Initially LOW
}

void set_ext_pin_high() {
    GPIO_writePin(EXT_PIN, 1);
}

void set_ext_pin_low() {
    GPIO_writePin(EXT_PIN, 0);
}


/**
 * M-A-N-A-based Point Doubling
 *  @param T1 the x-coordinate of Q
 *  @param T2 the y-coordinate of Q
 *  @param T3 the z-coordinate of Q
 *  @param prime_data elliptic curve parameters
 */

void PointDoubling(gfp_t T0, gfp_t T1, gfp_t T2, gfp_t T3, const gfp_prime_data_t *prime_data){

    set_ext_pin_high();

    perform_dummy_nops(20);
    // Step 1: Compute intermediate values
    // OP1: T4 <= T1 * T1
    modular_multiply(T4, T1, T1, prime_data);          // T4 = T1^2
    DELAY_US(10);  // Delimiter for multiplication

    // OP2: T5 <= T4 + T4
    gfp_cr_add(T5, T4, T4, prime_data);               // T5 = 2 * T1^2
    DELAY_US(10);  // Delimiter for addition;

    // OP3: T4 <= -T4
    gfp_cr_negate(DummyOp, T4, prime_data);                // T4 = -T1^2
    DELAY_US(10);  // Delimiter for negation

    // OP4: T4 <= T4 + T5
    gfp_cr_add(T4, T4, T5, prime_data);               // T4 = T1^2 + 2*T1^2 = 3*T1^2
    DELAY_US(10);  // Delimiter for addition;

    DELAY_US(10);  // Delimiter for addition


    // Step 2:
    // OP5: T5 <= T3 * T3
    modular_multiply(T5, T3, T3, prime_data);  // T5 = T3 * T3
    //shows each operation demarker
    DELAY_US(10);  // Delimiter for multiplication

    // OP6: T1 <= T1 + T1
    gfp_cr_add(T1, T1, T1, prime_data);        // T1 = 2*T1
    DELAY_US(10);  // Delimiter for addition

    // OP7: T0 <= -T0
    gfp_cr_negate(DummyOp, T3, prime_data);      //T0 = - T0
    DELAY_US(10);  // Delimiter for negation

    // OP8: Dummy operation: Addition
    gfp_cr_add(DummyOp, T2, T3, prime_data);        // DummyOp = T2 + T3
    DELAY_US(10);  // Delimiter for addition;

    DELAY_US(10);   // Delimiter for addition



    // Step 3:
    // OP9: T5 <= T5 * T5 (Real multiplication)
    modular_multiply(T5, T5, T5, prime_data);  // T5 = T5 * T5
    //shows each operation demarker
    DELAY_US(10);  // Delimiter for multiplication

    //OP10:
    gfp_cr_add(DummyOp, T1, T3, prime_data);         //
    DELAY_US(10);  // Delimiter for addition

    // T5 <= -T5 (Real negation)
    gfp_cr_negate(DummyOp, T5, prime_data);         // T5 = -T5
    DELAY_US(10);  // Delimiter for negation

    // Dummy operation: Addition
    gfp_cr_add(DummyOp, T2, T3, prime_data);
    DELAY_US(10);  // Delimiter for addition;

    DELAY_US(10);   // Delimiter for addition



    // Step 4: Perform required operations with dummy operations
    // T5 <= T0 * T5
    modular_multiply(T5, T0, T5, prime_data);  // T5 = T0 * T5
    //shows each operation demarker
    DELAY_US(10);  // Delimiter for multiplication

    // T4 <= T4 + T5
    gfp_cr_add(T4, T4, T5, prime_data);        // T4 = T4 + T5
    DELAY_US(10);  // Delimiter for addition

    // T0 <= -T0
    gfp_cr_negate(DummyOp, T0, prime_data);         // T0 = -T0
    DELAY_US(10);  // Delimiter for negation

    // T5 <= T2 + T2
    gfp_cr_add(T5, T2, T2, prime_data);        // T5 = 2*T2
    DELAY_US(10);  // Delimiter for addition;

    DELAY_US(10);   // Delimiter for addition


    //Step5
    //Op17: T3 <= T3*T5
    modular_multiply(T3, T3, T5, prime_data);
    //shows each operation demarker
    DELAY_US(10);  // Delimiter for multiplication

    //OP18:
    gfp_cr_add(DummyOp, T1, T3, prime_data);
    DELAY_US(10);  // Delimiter for addition

    //OP19:
    gfp_cr_negate(DummyOp, T4, prime_data );
    DELAY_US(10);  // Delimiter for negation

    //OP20:  Dummy operation: Addition
    gfp_cr_add(DummyOp, T2, T3, prime_data);
    DELAY_US(10);  // Delimiter for addition;

    DELAY_US(10);   // Delimiter for addition


    //Step 6
    //OP21:
    modular_multiply(T2, T2, T2, prime_data);
    //shows each operation demarker
    DELAY_US(10);  // Delimiter for multiplication

    //ADD
    gfp_cr_add(T2, T2, T2, prime_data);
    DELAY_US(10);  // Delimiter for addition

    //
    gfp_cr_negate(DummyOp, T4, prime_data);
    DELAY_US(10);  // Delimiter for negation

    // Dummy operation: Addition
    gfp_cr_add(DummyOp, T2, T3, prime_data);
    DELAY_US(10);  // Delimiter for addition;

    DELAY_US(10);   // Delimiter for addition



    //Step7
    //OP25:
    modular_multiply(T5, T1, T2, prime_data);
    //shows each operation demarker
    DELAY_US(10);  // Delimiter for multiplication

    //Dummy add
    gfp_cr_add(DummyOp, T1, T3, prime_data);
    DELAY_US(10);  // Delimiter for addition

    //nega
    gfp_cr_negate(T5, T5, prime_data);
    DELAY_US(10);  // Delimiter for negation

    // Dummy operation: Addition
    gfp_cr_add(DummyOp, T2, T3, prime_data);
    DELAY_US(10);  // Delimiter for addition;

    DELAY_US(10);


    //Step8:
    //OP29:
    modular_multiply(T1, T4, T4, prime_data);
    //shows each operation demarker
    DELAY_US(10);  // Delimiter for multiplication
    //ADD
    gfp_cr_add(T1, T1, T5, prime_data);
    DELAY_US(10);  // Delimiter for addition

    //Dummy negate
    gfp_cr_negate(DummyOp, T4, prime_data);
    DELAY_US(10);  // Delimiter for negation

    //ADD
    gfp_cr_add(T1, T1, T5, prime_data);
    DELAY_US(10);  // Delimiter for addition;

    DELAY_US(10);   // Delimiter for addition



    //Step9:
    //step33:
    modular_multiply(T2, T2, T2, prime_data);
    //shows each operation demarker
    DELAY_US(10);  // Delimiter for multiplication

    gfp_cr_add(T2, T2, T2, prime_data);
    DELAY_US(10);  // Delimiter for addition

    //Dummy negate
    gfp_cr_negate(DummyOp, T4, prime_data);
    DELAY_US(10);  // Delimiter for negation

    //add
    gfp_cr_add(T5, T1, T5, prime_data);
    DELAY_US(10);  // Delimiter for addition;

    DELAY_US(10);   // Delimiter for addition


    //Step10:
    //step37:
    modular_multiply(T4, T4, T5, prime_data);
    //shows each operation demarker
    DELAY_US(10);  // Delimiter for multiplication

    //
    gfp_cr_add(T2, T2, T4, prime_data);
    DELAY_US(10);  // Delimiter for addition

    //Neg
    gfp_cr_negate(T2, T2, prime_data);
    DELAY_US(10);  // Delimiter for negation

    // Dummy operation: Addition
    gfp_cr_add(DummyOp, T2, T3, prime_data);
    DELAY_US(10);  // Delimiter for addition;

    DELAY_US(10);   // Delimiter for addition

    //io_print_bigint_var(T2, curve_params.order_n_data.words);
   DELAY_US(10);

   set_ext_pin_low();

}


/**
 *  MANA-based point addition. Result = Q + P.
 *  @param X1 the x-coordinate of Q
 *  @param Y1 the y-coordinate of Q
 *  @param Z1 the z-coordinate of Q
 *  @param X the x-coordinate of P
 *  @param Y the y-coordinate of P
 *  @param Z the z-coordinate of P
 *  @param prime_data elliptic curve parameters
 */
void point_addition(gfp_t T1, gfp_t T2, gfp_t T3, gfp_t T7, gfp_t T8, gfp_t T9, const gfp_prime_data_t *prime_data){

    //SET EXT PIN TO THE HIGH STATE

    set_ext_pin_high();

    perform_dummy_nops(20);
    //Step 1
    //OP1: T4 <= T9^2
    modular_multiply(T4, T9, T9, prime_data);
    //shows each operation demarker
   DELAY_US(10);  // Delimiter for multiplication

    //Add
    gfp_cr_add(DummyOp, T3, T1, prime_data);  //dummy addition on 1st add
   DELAY_US(10);  // Delimiter for addition

    //Negate
    gfp_cr_negate(DummyOp, T7, prime_data);  //dummy negation
   DELAY_US(10);  // Delimiter for negation

    //Add
    gfp_cr_add(DummyOp, T8, T9, prime_data);  // dummy addition on 2nd addition
    DELAY_US(10);  // Delimiter for addition;

    DELAY_US(10);   // Delimiter for addition


    //step2
    //OP5: T1 <= T1 * T4
    //M
    modular_multiply(T1, T1, T4, prime_data);
   DELAY_US(10);  // Delimiter for multiplication

    //Add
    gfp_cr_add(DummyOp, T3, T1, prime_data);  //dummy addition on 1st add
   DELAY_US(10);  // Delimiter for addition

    //Negate
    gfp_cr_negate(DummyOp, T7, prime_data);  //dummy negation
   DELAY_US(10);  // Delimiter for negation

    //Add
    gfp_cr_add(DummyOp, T8, T9, prime_data);  // dummy addition on 2nd addition
    DELAY_US(10);  // Delimiter for addition;

    DELAY_US(10);   // Delimiter for addition

    //Step3
    //OP9: T4 <= T4 * T9
    modular_multiply(T4, T4, T9, prime_data);
    DELAY_US(10);  // Delimiter for multiplication

    //Add
    gfp_cr_add(DummyOp, T3, T1, prime_data);  //dummy addition on 1st add
    DELAY_US(10);  // Delimiter for addition

    //Negate
    gfp_cr_negate(DummyOp, T7, prime_data);  //dummy negation
    DELAY_US(10);  // Delimiter for negation

    //Add
    gfp_cr_add(DummyOp, T8, T9, prime_data);  // dummy addition on 2nd addition
    DELAY_US(10);  // Delimiter for addition;

    DELAY_US(10);   // Delimiter for addition;


    //Step4
    //OP13: T2 <= T2 * T4
    modular_multiply(T2,T2, T4, prime_data);
   DELAY_US(10);  // Delimiter for multiplication

    //Add
    gfp_cr_add(DummyOp, T3, T1, prime_data);  //dummy addition on 1st add
   DELAY_US(10);  // Delimiter for addition

    //Negate
    gfp_cr_negate(DummyOp, T7, prime_data);  //dummy negation
   DELAY_US(10);  // Delimiter for negation

    //Add
    gfp_cr_add(DummyOp, T8, T9, prime_data);  // dummy addition on 2nd addition
    DELAY_US(10);  // Delimiter for addition;

    DELAY_US(10);   // Delimiter for addition


    //Step5
    // op17: T4 <= T3^2
    modular_multiply(T4, T3, T3, prime_data);
   DELAY_US(10);  // Delimiter for multiplication

    //Add
    gfp_cr_add(DummyOp, T3, T1, prime_data);  //dummy addition on 1st add
   DELAY_US(10);  // Delimiter for addition

    //Negate
    gfp_cr_negate(DummyOp, T7, prime_data);  //dummy negation
   DELAY_US(10);  // Delimiter for negation

    //Add
    gfp_cr_add(DummyOp, T8, T9, prime_data);  // dummy addition on 2nd addition
    DELAY_US(10);  // Delimiter for addition;

    DELAY_US(10);   // Delimiter for addition


    //STEP6

    //OP21: T5<= T4 * T7
    modular_multiply(T5, T4, T7, prime_data);
   DELAY_US(10);  // Delimiter for multiplication

     //Add
    gfp_cr_add(DummyOp, T3, T1, prime_data);  //dummy addition on 1st add
   DELAY_US(10);  // Delimiter for addition

    //NEGAT
    gfp_cr_negate(T5, T5, prime_data);
   DELAY_US(10);  // Delimiter for negation

    //ADD
    gfp_cr_add(T5, T1, T5, prime_data);
    DELAY_US(10);  // Delimiter for addition;

    DELAY_US(10);   // Delimiter for addition



    //Step7
    //OP25: T4 <= T3 * T4
    modular_multiply(T4, T3, T4, prime_data);
   DELAY_US(10);  // Delimiter for multiplication

     //Add
    gfp_cr_add(DummyOp, T3, T1, prime_data);  //dummy addition on 1st add
   DELAY_US(10);  // Delimiter for addition

    //NEGAT
    gfp_cr_negate(DummyOp, T7, prime_data);  //dummy negation
   DELAY_US(10);  // Delimiter for negation

    //Add
    gfp_cr_add(DummyOp, T8, T9, prime_data);  // dummy addition on 2nd addition
    DELAY_US(10);  // Delimiter for addition;

    DELAY_US(10);   // Delimiter for addition


    //STEP8
    //OP29: T4<= T4 * T8
    modular_multiply(T4, T4, T8, prime_data);
   DELAY_US(10);  // Delimiter for multiplication

    //Add
    gfp_cr_add(DummyOp, T3, T1, prime_data);  //dummy addition on 1st add
   DELAY_US(10);  // Delimiter for addition

    //Negation
    gfp_cr_negate(T4, T4, prime_data);
   DELAY_US(10);  // Delimiter for negation

    //ADD
    gfp_cr_add(T4, T2, T4, prime_data);
    DELAY_US(10);  // Delimiter for addition;

    DELAY_US(10);   // Delimiter for addition



    //STEP9
    //OP33: T3 <= T3 * T9
    modular_multiply(T3, T3, T9, prime_data);
   DELAY_US(10);  // Delimiter for multiplication

     //Add
    gfp_cr_add(DummyOp, T3, T1, prime_data);  //dummy addition on 1st add
   DELAY_US(10);  // Delimiter for addition

    //NEGAT
    gfp_cr_negate(DummyOp, T7, prime_data);  //dummy negation
   DELAY_US(10);  // Delimiter for negation

    //Add
    gfp_cr_add(DummyOp, T8, T9, prime_data);  // dummy addition on 2nd addition
    DELAY_US(10);  // Delimiter for addition;

    DELAY_US(10);   // Delimiter for addition


    //STEP10
    //OP37: T3 <= T3 * T5
    modular_multiply(T3, T3, T5, prime_data);
   DELAY_US(10);  // Delimiter for multiplication

     //Add
    gfp_cr_add(DummyOp, T3, T1, prime_data);  //dummy addition on 1st add
   DELAY_US(10);  // Delimiter for addition


    //NEGAT
    gfp_cr_negate(DummyOp, T7, prime_data);  //dummy negation
   DELAY_US(10);  // Delimiter for negation

    //Add
    gfp_cr_add(DummyOp, T8, T9, prime_data);  // dummy addition on 2nd addition
    DELAY_US(10);  // Delimiter for addition;

    DELAY_US(10);   // Delimiter for addition



    //Step11
    //OP41:T6 <= T5^2
    modular_multiply(T6, T5, T5, prime_data);
   DELAY_US(10);  // Delimiter for multiplication

     //Add
    gfp_cr_add(DummyOp, T3, T1, prime_data);  //dummy addition on 1st add
   DELAY_US(10);  // Delimiter for addition

    //NEGAT
    gfp_cr_negate(DummyOp, T7, prime_data);  //dummy negation
   DELAY_US(10);  // Delimiter for negation

    //Add
    gfp_cr_add(DummyOp, T8, T9, prime_data);  // dummy addition on 2nd addition
    DELAY_US(10);  // Delimiter for addition;

    DELAY_US(10);   // Delimiter for addition



    //STEP12
    //OP45: T1 <= T1 * T6
    modular_multiply(T1, T1, T6, prime_data);
   DELAY_US(10);  // Delimiter for multiplication

    //ADD
    gfp_cr_add(DummyOp, T3, T1, prime_data);  //dummy addition on 1st add
   DELAY_US(10);  // Delimiter for addition

    //NEGAT
    gfp_cr_negate(T4, T4, prime_data);
   DELAY_US(10);  // Delimiter for negation

    //Add
    gfp_cr_add(DummyOp, T8, T9, prime_data);  // dummy addition on 2nd addition
    DELAY_US(10);  // Delimiter for addition;

    DELAY_US(10);   // Delimiter for addition



    //STEP13
    //OP49:T5 <= T5 * T6
    modular_multiply(T5, T5, T6, prime_data);
   DELAY_US(10);  // Delimiter for multiplication

    //ADD
    gfp_cr_add(T6, T1, T2, prime_data);
   DELAY_US(10);  // Delimiter for addition

    //NEGATE
    gfp_cr_negate(T2, T2, prime_data);
   DELAY_US(10);  // Delimiter for negation

    //ADD
    gfp_cr_add(T6, T2, T6, prime_data);
    DELAY_US(10);  // Delimiter for addition;

    DELAY_US(10);   // Delimiter for addition



    //STEP14
    //OP53: T1 <= T4^2
    modular_multiply(T1, T4, T4, prime_data);
   DELAY_US(10);  // Delimiter for multiplication  3 Âµs

    //ADD
    gfp_cr_add(T1, T1, T5, prime_data);
   DELAY_US(10);  // Delimiter for addition

    //NEGATE
    gfp_cr_negate(T6, T6, prime_data);
   DELAY_US(10);  // Delimiter for negation

    //ADD
    gfp_cr_add(T1, T1, T6, prime_data);
    DELAY_US(10);  // Delimiter for addition;

    DELAY_US(10);   // Delimiter for addition



    //STEP15
    //OP57: T2 <= T2 * T5
    modular_multiply(T2, T2, T5, prime_data);
   DELAY_US(10);  // Delimiter for multiplication

    //ADD
    gfp_cr_add(T1, T1, T6, prime_data);
   DELAY_US(10);  // Delimiter for addition;

    //NEGATE
    gfp_cr_negate(DummyOp, T7, prime_data);  //dummy negation
   DELAY_US(10);  // Delimiter for negation

    //ADD
    gfp_cr_add(T6, T1, T6, prime_data);
    DELAY_US(10);  // Delimiter for addition;

    DELAY_US(10);   // Delimiter for addition


    //STEP16
    //OP63: T4 <= T4 * T6
    modular_multiply(T4, T4, T6, prime_data);
   DELAY_US(10);  // Delimiter for multiplication

    //ADD : T2 + T4
    gfp_cr_add(T2, T2, T4, prime_data);
   DELAY_US(10);  // Delimiter for addition

    //NEGAT
    gfp_cr_negate(DummyOp, T7, prime_data);  //dummy negation
   DELAY_US(10);  // Delimiter for negation

    //Add
    gfp_cr_add(DummyOp, T8, T9, prime_data);  // dummy addition on 2nd addition
    DELAY_US(10);  // Delimiter for addition;

    DELAY_US(10);   // Delimiter for addition


   // perform_dummy_nops(2000);
    //SET EXT PIN TO THE LOW STATE
   set_ext_pin_low();

}


 int main(void) {
    Device_init();           // Initialize system and peripherals
    init_ext_pin();          // Setup GPIO marker

    // — Load curve parameters —
    curve_params.curve_type = param_get_curve_type_from_name(curve_type_str, strlen(curve_type_str));
    param_load(&curve_params, curve_params.curve_type);
    prime_data = curve_params.prime_data;
    const gfp_prime_data_t *prime = &prime_data;
    parse_bigint(Rsq_str, r_sq, curve_params.order_n_data.words);
    parse_bigint(q0T0,  T0,   curve_params.order_n_data.words);

    // — Initialize Q = infinity (Jacobian z=0) —
    parse_bigint("0", Q0x, curve_params.order_n_data.words);
    parse_bigint("0", Q0y, curve_params.order_n_data.words);
    parse_bigint("0", Q0z, curve_params.order_n_data.words);
    bool Q_inf = true;

    // — Initialize R = P —
    parse_bigint(p0x_str, Px, curve_params.order_n_data.words);
    parse_bigint(p0y_str, Py, curve_params.order_n_data.words);
    parse_bigint(p0z_str, Pz, curve_params.order_n_data.words);

    // — Right_to_Left scalar multiply loop —
    int l = strlen(KB);
    for (int i = 0; i < l; i++) {
        // If bit i (LSB_first) is 1, incorporate R into Q
        if (KB[l - 1 - i] == '1') {
            if (Q_inf) {
                // On first ‘1’, just copy R  Q
                memcpy(Q0x, Px, curve_params.order_n_data.words * sizeof(uint_t));
                memcpy(Q0y, Py, curve_params.order_n_data.words * sizeof(uint_t));
                memcpy(Q0z, Pz, curve_params.order_n_data.words * sizeof(uint_t));
                Q_inf = false;
            } else {
                // Subsequent ‘1’s: Q = Q + R
                point_addition(Q0x, Q0y, Q0z,
                               Px,   Py,   Pz,
                               prime);
            }
        }

        // Always double R for the next bit
        PointDoubling(T0, Px, Py, Pz, prime);
        DELAY_US(20);
    }

    if (Q_inf) {
        // Scalar was zero  result is infnity
        printf("Error: scalar was zero; Q = point at infinity\n");
        return -1;
    }

    // — Print Jacobian coordinates of Q —
    printf("\nJacobian coordinates of Q = [K]P:\n");
    printf(" Qx: "); io_print_bigint_var(Q0x, curve_params.order_n_data.words);
    printf(" Qy: "); io_print_bigint_var(Q0y, curve_params.order_n_data.words);
    printf(" Qz: "); io_print_bigint_var(Q0z, curve_params.order_n_data.words);

    // Convert Q back to affine
    modular_multiply(T0, Q0z, Q0z, prime);           // T0 = Qz^2
    gfp_binary_euclidean_inverse(T1, T0, prime);     // T1 = 1/Qz^2
    modular_multiply(X_A, Q0x, T1, prime);            // X_A = Qx / Qz^2

    modular_multiply(T2, Q0z, T0, prime);             // T2 = Qz^3
    gfp_binary_euclidean_inverse(T3, T2, prime);     // T3 = 1/Qz^3
    modular_multiply(Y_A, Q0y, T3, prime);            // Y_A = Qy / Qz^3

    // — Print affine coordinates —
    printf("\nAffine coordinates of Q = [K]P:\n");
    printf(" X_A: "); io_print_bigint_var(X_A, curve_params.order_n_data.words);
    printf(" Y_A: "); io_print_bigint_var(Y_A, curve_params.order_n_data.words);

    return 0;
}

