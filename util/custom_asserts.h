#pragma once

#include <iostream>
#include <string>
#include <cassert>

inline void ensure(const bool assumption) noexcept {
    if (!assumption) exit(1);
}

//#define Assert(assumption) assert(assumption)
#define AssertMsg(assumption, msg) assert((assumption) || (std::cout << "\n\033[31mASSERTION FAILED: " << msg << "\033[0m\nFile: " << __FILE__ << "\nLine: " << __LINE__ << "\n" << std::flush && false))
#define Assert(assumption) assert((assumption) || (std::cout << "\n\033[31mASSERTION FAILED: " << "\033[0m\nFile: " << __FILE__ << "\nLine: " << __LINE__ << "\n" << std::flush && false))

