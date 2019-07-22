#pragma once
#include <chrono>
namespace whfc {
    typedef std::chrono::duration<double, std::micro> duration;
    typedef std::chrono::duration<double> second_duration;

    inline std::chrono::time_point<std::chrono::high_resolution_clock> time_now() {
		return std::chrono::high_resolution_clock::now();
	}

    template<typename resolution>
    std::chrono::duration<double, std::milli> inMilliseconds(std::chrono::duration<double, resolution> dur) {
    	return dur;
    }
	template<typename resolution>
	std::chrono::duration<double> inSeconds(std::chrono::duration<double, resolution> dur) {
		return dur;
	}

}//namespace
