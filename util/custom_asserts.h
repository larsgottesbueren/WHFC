#pragma once

#include <iostream>
#include <string>
#include <cassert>

inline void ensure(const bool assumption) noexcept {
    if (!assumption) exit(1);
}

#define AssertMsg(assumption, msg) assert(( (void)(msg), (assumption) ))
#define Assert(assumption) assert((assumption))

