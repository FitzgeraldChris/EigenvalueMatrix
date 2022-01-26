// Dr. Block's wrapper class for the Mersenne Twister (mt19937) PRNG
// !!!!!!!!!!!!!!!!!!!!
// !! Requires C++14 !!
// !!!!!!!!!!!!!!!!!!!!

#ifndef MT_HPP
#define MT_HPP

#include <cstdlib>
#include <random>
#include <array>
#include <algorithm>
#include <functional>

namespace block
{
    class MT
    {
        private:
        
        std::random_device rd; // doesn't work with minGW; use <chrono> instead
        std::array<int, std::mt19937::state_size> seed_data;
        std::mt19937 gen;
        
        void init_fancy()
        {
            std::generate_n(seed_data.data(), seed_data.size(), std::ref(rd));
            std::seed_seq seq(std::begin(seed_data), std::end(seed_data));
            gen.seed(seq);
        }
        
        void init_basic(int seed)
        {
            gen.seed(seed);
        }
        
        public:
        
        MT() // use for effective random seeding
        {
            init_fancy();
        }
        
        MT(int seed) // use when you want reproducible results
        {
            if(seed < 0)
                init_fancy();
            else
                init_basic(seed);
        }
        
        // gives a uniform distribution for integers between min and max, inclusive of both
        auto rand_int(unsigned int min, unsigned int max)
        {
            return std::bind(std::uniform_int_distribution<int>(min,max),std::ref(gen));
        }
        
        // gives a uniform distribution for real numbers between min and max;
        // inclusive of min, but not of max
        auto rand_real(double min, double max)
        {
            return std::bind(std::uniform_real_distribution<double>(min,max),std::ref(gen));
        }
        
        // same as previous function with min = 0.0 and max = 1.0
        auto rand_real()
        {
            return rand_real(0.0, 1.0);
        }

    };
}

#endif
