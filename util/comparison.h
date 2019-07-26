#pragma once

namespace Util {
	template<typename T>
	class InvertComparison {
	public:
		T value;

		explicit InvertComparison(const T& value) : value(value) { }

		InvertComparison& operator=(const T& value) {
			this->value = value;
			return *this;
		}

		bool operator<(const InvertComparison<T>& o) const {
			return value > o.value;
		}

		bool operator>(const InvertComparison<T>& o) const {
			return value < o.value;
		}
	};
}