#include <vector>
#include <iostream>
#include <algorithm>
#include <random>

using namespace std;

enum Color { red, green, blue };

int main(int argc, char const *argv[])
{

    // Define the values and their respective probabilities
    std::vector<int> values = {1, 2, 3};
    std::vector<double> probabilities = {0.5, 1.0/3.0, 1.0/6.0};

    // Define the discrete distribution
    std::discrete_distribution<int> dist(probabilities.begin(), probabilities.end());

    // Seed the random number generator
    std::random_device rd;
    std::mt19937 gen(rd());

    // Sample from the distribution
    int sample = values[dist(gen)];

    std::cout << "Sampled value: " << sample << std::endl;

    return 0;
}
