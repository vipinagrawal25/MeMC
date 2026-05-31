#include "random_gen.h"

std::default_random_engine RandomGenerator::generator;

void RandomGenerator::init(unsigned int seed) {
    generator.seed(seed);
}

double RandomGenerator::generateGaussian(double mean, double stddev) {
    std::normal_distribution<double> distribution(mean, stddev);
    return distribution(generator);
}

double RandomGenerator::generateUniform(double min, double max) {
    std::uniform_real_distribution<double> distribution(min, max);
    return distribution(generator);
}

int RandomGenerator::intUniform(int min, int max) {
    std::uniform_int_distribution<int> distribution(min, max);
    return distribution(generator);
}


/* RandomGenerator::RandomGenerator(unsigned int seed, bool is_random) : generator(std::random_device()()) { */
/*     if(!is_random){ */
/*         generator.seed(seed); */
/*     } */
/* } */

/* double RandomGenerator::generateGaussian(double mean, double stddev) { */
/*     std::normal_distribution<double> distribution(mean, stddev); */
/*     return distribution(generator); */
/* } */

/* double RandomGenerator::generateUniform(double min, double max) { */
/*     std::uniform_real_distribution<double> distribution(min, max); */
/*     return distribution(generator); */
/* } */

