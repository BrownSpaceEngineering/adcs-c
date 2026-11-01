#ifndef TEST_H
#define TEST_H

// Runs every test suite. Failures are printed and counted; nothing aborts.
// Returns the number of failed checks.
int test_run_all(void);

// Individual suites
void test_linalg(void);
void test_matrix_product(void);
void test_quaternion(void);
void test_quest(void);
void test_ukf_internals(void);
void test_ukf_behaviour(void);
void test_iteration_1vec(void);
void test_iteration_2vec(void);
void test_bdot(void);
void test_pd(void);
void test_torque_and_currents(void);
void test_kepler(void);
void test_orbital_to_eci(void);
void test_ecef_to_eci(void);
void test_sun_vec(void);
void test_magnetosphere(void);
void test_photodiodes(void);
void test_pointing(void);
void test_filters(void);
void test_body(void);
void test_closest_approach(void);
// Compares against include/matlab_reference.h (tools/matlab/generate_matlab_reference.m)
void test_matlab_parity(void);
#endif
