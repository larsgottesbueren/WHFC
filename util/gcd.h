#pragma once

#include <vector>

namespace whfc {
	class GreatestCommonDivisor {
	public:
		template<typename T>
		static T euclid(T a, T b) {
			while (b != 0) {
				T temp = a % b;
				a = b;
				b = temp;
			}
			return a;
		}

		template<typename T>
		static T compute(const std::vector<T>& numbers) {
			T a = numbers[0];
			for (size_t i = 1; i < numbers.size() && a > 1; ++i)
				a = euclid(a, numbers[i]);
			return a;
		}
	};
}