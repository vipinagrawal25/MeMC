#ifndef RANDOM_GEN_H
#define RANDOM_GEN_H

#include <random>

class RandomGenerator {
public:
    static void init(unsigned int seed);

    // Generate a Gaussian (normal) random number
    static double generateGaussian(double mean, double stddev);

    // Generate a uniform random number in the range [min, max]
    static double generateUniform(double min, double max);

    static int intUniform(int min, int max);

private:
    static std::default_random_engine generator;
};

#endif // RANDOM_GENERATOR_H
       /* // */

/* class RandomGenerator { */
/* public: */
    /* RandomGenerator(unsigned int seed, bool is_random); */

    /* // Generate a Gaussian (normal) random number */
    /* double generateGaussian(double mean, double stddev); */

    /* // Generate a uniform random number in the range [min, max] */
    /* double generateUniform(double min, double max); */

/* private: */
    /* std::default_random_engine generator; */
/* }; */

/* #endif // RANDOM_GENERATOR_H */

