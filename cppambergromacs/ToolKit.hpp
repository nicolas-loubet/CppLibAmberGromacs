#ifndef TOOLKIT_HPP
#define TOOLKIT_HPP

#include <stdexcept>
#include <iostream>
#include <chrono>
#include <thread>
#include <vector>

/**
 * Version: April 2025
 * Author: Nicolás Loubet
 */

namespace ToolKit {
    constexpr float k_B= -1.380649*6.02214076*0.001; // Boltzmann constant in kJ/(mol*K)

    /**
     * STRUCT that contains the pointer to an array of int and the size of that array
     */
    struct ArrInt {
        int* arr;
        int size= 0;
    };

    /**
     * STRUCT that contains the pointer to an array of float and the size of that array
     */
    struct ArrFloat {
        float* arr;
        int size= 0;
    };

    struct FlaggedArrFloat {
        int size;
        float* arr;
        bool flag;
    
        FlaggedArrFloat(int size, float* arr, bool flag) {
            this->size= size;
            this->arr= arr;
            this->flag= flag;
        }
    };


    /**
     * Function that returns the position in a set of bins of a value (for a histogram)
     * @param value The value to analyze
     * @param LIMIT_MIN The value of the first bin
     * @param LIMIT_MAX The value of the last bin
     * @param N_BINS The number of bins
     * @param trim If true, when the position if less than 0 or more than N_BINS, it returns the nearest position. If false, in the case of "out of bounds" it thwors an exception
     * @throws std::out_of_range if trim is false and the position is out of bounds
     * @return position of the bin
     */
    int getBinPosition(float value, const float LIMIT_MIN, const float LIMIT_MAX, const int N_BINS, bool trim= false) {
        int pos= int((value-LIMIT_MIN)*N_BINS/(LIMIT_MAX-LIMIT_MIN));
        if(trim) {
            return std::max(0, std::min(N_BINS-1,pos));
        } else {
            if(pos < 0 || pos >= N_BINS)
                throw std::out_of_range("The calculated position is out of bounds");
            return pos;
        }
    }

    /**
     * Function that measures the time elapsed between the start and the end of an execution
     * @param f The function to execute
     * @return The time elapsed in seconds
     * @tparam Func The type of the function to execute
     */  
    template<typename Func>
    long int takeTime(Func f) {
        auto init_time= std::chrono::high_resolution_clock::now();
        f();
        auto end_time= std::chrono::high_resolution_clock::now();

        long int elapsed= std::chrono::duration_cast<std::chrono::seconds> (end_time - init_time).count();
        std::cout << "Finished execution in " << ((int) elapsed/3600) << "h:" << ((int) (elapsed%3600)/60) << "min:" << elapsed%60 << "s (" << elapsed << "s)" << std::endl;
        return elapsed;
    }

    /**
     * Function that executes a function in parallel for each element of a list
     * @param f The function to execute (first argument must be one of the elements of the list, type T; second argument must be a reference to the result of the function, type Res)
     * @param list The list of elements (for example, the directories to be analyzed)
     * @param args The arguments of the function
     * @tparam Func The type of the function to execute
     * @tparam T The type of the elements of the list
     * @tparam Res The type of the return value of the function
     * @tparam Args The types of the arguments of the function
     */
    template<typename Func, typename T, typename Res, typename... Args>
    void parallel(Func f, std::vector<T>& list, std::vector<Res>& res, Args... args) {
        std::vector<std::thread> threads;
        for(size_t i= 0; i < list.size(); i++)
            threads.emplace_back(  [&,i]() { f(list[i], res[i], args...); }  );
        for(std::thread& thread: threads)
            thread.join();
    }

}

#endif
