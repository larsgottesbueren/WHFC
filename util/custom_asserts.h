#pragma once

#include <iostream>
#include <string>
#include <cassert>

inline void ensure(const bool assumption) noexcept {
    if (!assumption) exit(1);
}

//#define Assert(assumption) assert(assumption)
//These don't seem necessary anymore, since we not get appropriate line numbers and filenames
#define AssertMsg(assumption, msg) assert(( (void)(msg), (assumption) ))
#define Assert(assumption) assert((assumption))

