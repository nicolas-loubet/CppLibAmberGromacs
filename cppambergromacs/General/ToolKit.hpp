#ifndef TOOLKIT_HPP
#define TOOLKIT_HPP

/**
 * Version: April 2025
 * Author: Nicolás Loubet
 */

#include <stdexcept>
#include <iostream>
#include <chrono>
#include <thread>
#include <vector>
#include <algorithm>
#include <mutex>

static std::mutex mtx;

namespace ToolKit {
	constexpr Real k_B= -1.380649*6.02214076*0.001; // Boltzmann constant in kJ/(mol*K)

	/**
	 * STRUCT that contains the pointer to an array of int and the size of that array
	 */
	struct ArrInt {
		int* arr;
		int size= 0;
	};

	/**
	 * STRUCT that contains the pointer to an array of Real and the size of that array
	 */
	struct ArrFloat {
		Real* arr;
		int size= 0;
	};

	/**
	 * STRUCT that contains the pointer to an array of Real and the size of that array
	 */
	struct FlaggedArrFloat {
		int size;
		Real* arr;
		bool flag;
	
		FlaggedArrFloat(int size, Real* arr, bool flag) {
			this->size= size;
			this->arr= arr;
			this->flag= flag;
		}
	};


	/**
	 * Removes all spaces from the string
	 * @param _s The string to strip
	 * @return The stripped string
	 */
	std::string strip(std::string _s)
	{
		_s.erase(std::remove(_s.begin(), _s.end(), ' '), _s.end());
		return(_s);
	}


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
	int getBinPosition(Real value, const Real LIMIT_MIN, const Real LIMIT_MAX, const int N_BINS, bool trim= false) {
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

	/**
	 * Function that executes a function for each element of a list, one after the other
	 * @param f The function to execute (first argument must be one of the elements of the list, type T; second argument must be a reference to the result of the function, type Res)
	 * @param list The list of elements (for example, the directories to be analyzed)
	 * @param args The arguments of the function
	 * @tparam Func The type of the function to execute
	 * @tparam T The type of the elements of the list
	 * @tparam Res The type of the return value of the function
	 * @tparam Args The types of the arguments of the function
	 */
	template<typename Func, typename T, typename Res, typename... Args>
	void serialExecution(Func f, std::vector<T>& list, std::vector<Res>& res, Args... args) {
		for(size_t i= 0; i < list.size(); i++)
			f(list[i], res[i], args...);
	}

	#ifndef CONSOLE_WIDTH
		#define CONSOLE_WIDTH 150
	#endif

	/**
	 * Prints a charging bar for each parallel execution
	 * @param i The current iteration
	 * @param N The total number of iterations
	 * @param pos The position of the bar (0-20)
	 * @param total_pos The total number of positions
	 * @param symbol The symbol to print
	 */
	void printPercentageParallel(int i, int N, int pos, int total_pos, const char symbol = '=') {
		#ifdef NOHUP
			std::cout << "Progress: "+std::to_string(i)+"/"+std::to_string(N)+" of process "+std::to_string(pos+1)+"/"+std::to_string(total_pos)+"\n" << std::flush;
		#else
			const int bar_width= CONSOLE_WIDTH/total_pos -1;

			static std::vector<std::string> bars(total_pos, std::string(bar_width, ' '));
			int filled= (i * (bar_width-2)) / N;

			bars[pos]= "[" + std::string(filled, symbol) + std::string(bar_width-2 - filled, ' ') + "]";

			std::lock_guard<std::mutex> lock(mtx);
			std::cout << "\r";
			for(int p= 0; p < total_pos; p++)
				std::cout << bars[p];
			std::cout.flush();
		#endif
	}


}

#endif
