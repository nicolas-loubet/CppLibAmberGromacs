#ifndef SORTER_HPP
#define SORTER_HPP

/**
 * Version: April 2025
 * Author: Nicolás Loubet
 */

#include <vector>
#include <algorithm>
#include <cstddef>

namespace Sorter {
	enum class Order { Ascending, Descending };

	/**
	 * Orders a vector of type T in ascending order (least to greatest).
	 * @param list The vector to be sorted.
	 * @tparam T The type of the elements in the vector.
	 */
	template<typename T>
	inline void sort(std::vector<T>& list) {
		std::sort(list.begin(), list.end());
	}

	/**
	 * Orders a vector of type T in ascending/descending order.
	 * @param list The vector to be sorted.
	 * @param order The order to sort the vector. Use Sorter::Order::Ascending or Sorter::Order::Descending.
	 * @tparam T The type of the elements in the vector.
	 */
	template<typename T>
	inline void sort(std::vector<T>& list, Order order) {
		if(order == Order::Descending)
			std::sort(list.begin(), list.end(), std::greater<T>());
		else
			sort(list);
	}

	/**
	 * Orders an array of type T in ascending order (least to greatest).
	 * @param arr The array to be sorted.
	 * @param size The size of the array.
	 * @tparam T The type of the elements in the array.
	 */
	template<typename T>
	inline void sort(T* arr, size_t size) {
		std::sort(arr, arr + size);
	}

	/**
	 * Orders an array of type T in ascending/descending order.
	 * @param arr The array to be sorted.
	 * @param size The size of the array.
	 * @param order The order to sort the vector. Use Sorter::Order::Ascending or Sorter::Order::Descending.
	 * @tparam T The type of the elements in the array.
	 */
	template<typename T>
	inline void sort(T* arr, size_t size, Order order) {
		if(order == Order::Descending)
			std::sort(arr, arr + size, std::greater<T>());
		else
			sort(arr, size);
	}

	/**
	 * Orders two vectors of type T1 and T2 based on the values in the first vector.
	 * @param values The first vector to be sorted.
	 * @param indexes The second vector to be sorted based on the order of the first vector.
	 * @tparam T1 The type of the elements in the first vector.
	 * @tparam T2 The type of the elements in the second vector.
	 */
	template<typename T1, typename T2>
	inline void cosort(std::vector<T1>& values, std::vector<T2>& indexes) {
		std::vector<std::pair<float,int>> combined;
		for(size_t i= 0; i < values.size(); i++)
			combined.emplace_back(values[i],indexes[i]);

		std::sort(combined.begin(), combined.end());

		for(size_t i= 0; i < combined.size(); i++) {
			values[i]= combined[i].first;
			indexes[i]= combined[i].second;
		}
	}

	/**
	 * Orders two vectors of type T1 and T2 based on the values in the first vector.
	 * @param values The first vector to be sorted.
	 * @param indexes The second vector to be sorted based on the order of the first vector.
	 * @param order The order to sort the vector. Use Sorter::Order::Ascending or Sorter::Order::Descending.
	 * @tparam T1 The type of the elements in the first vector.
	 * @tparam T2 The type of the elements in the second vector.
	 */
	template<typename T1, typename T2>
	inline void cosort(std::vector<T1>& values, std::vector<T2>& indexes, Order order) {
		if(values.size() != indexes.size()) {
			throw std::invalid_argument("The sizes of the two vectors must be equal.");
		}

		std::vector<std::pair<float,int>> combined;
		for(size_t i= 0; i < values.size(); i++)
			combined.emplace_back(values[i],indexes[i]);

		if(order == Order::Descending) {
			std::sort(combined.begin(), combined.end(),
					[](const std::pair<T1,T2>& a, const std::pair<T1,T2>& b) {
						return a.first > b.first;
					});
		} else {
			std::sort(combined.begin(), combined.end(),
					[](const std::pair<T1,T2>& a, const std::pair<T1,T2>& b) {
						return a.first < b.first;
					});
		}

		for(size_t i= 0; i < combined.size(); i++) {
			values[i]= combined[i].first;
			indexes[i]= combined[i].second;
		}
	}

}

#endif
