#ifndef PACE_RANDOM_H
#define PACE_RANDOM_H

// PACE4_SEED supplies the recorded external seed (1..2147483646).
// Each calculation reinitializes its stream; leaving it unset uses system entropy.
bool InitRandom();
unsigned int CurrentRandomSeed();
double MyRandom();
// One-based cumulative distribution; consume one draw for the entire search.
int SampleInitialSpin(const double* cumulative, int count);
#endif
