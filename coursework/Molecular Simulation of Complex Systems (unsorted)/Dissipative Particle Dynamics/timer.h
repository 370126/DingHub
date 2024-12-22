/************************************************************************************
  Thread-safe Timer class with high resolution
  Copyright (c) 2010- Jinglei Hu < jingleihu@gmail.com >
************************************************************************************/
#ifndef _TIMER_H_
#define _TIMER_H_

#include <stdint.h>

/**
 \brief High Resolution Timer class adapted from Intel TBB tick_count.h which can be found at 
 \brief http://www.threadingbuildingblocks.org/files/documentation/a00402.html
 */
#if _WIN32||_WIN64
#include <stdio.h>
#include <windows.h>
/**
 \brief Microsecond Timer under Windows
 */
class cTimer {
public:
  cTimer() : t0(time()) {
    LARGE_INTEGER Ticks;
    QueryPerformanceFrequency(&Ticks);
    TicksPerSecond = uint64_t(Ticks.QuadPart);
  }
  inline void restart() { t0 = time(); }
  inline double elapsed() { return double(time() - t0) / TicksPerSecond; }
private:
  uint64_t t0;
  uint64_t TicksPerSecond;
  inline uint64_t time() {
    LARGE_INTEGER Ticks;
    QueryPerformanceCounter(&Ticks);
    return uint64_t(Ticks.QuadPart);
  }
};
#elif __linux__
#include <time.h>
/*!
 \brief Nanosecond Timer under Linux
 \brief option -lrt(real time library) is required to link librt.a!
 */
class cTimer {
public:
  cTimer() : t0(time()) { }
  inline void restart() { t0 = time(); }
  inline double elapsed() { return double((time() - t0)) * 1.0E-09; }
private:
  uint64_t t0;
  inline uint64_t time() {
    struct timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    return 1000000000ULL * uint64_t(ts.tv_sec) + uint64_t(ts.tv_nsec);
  }
};
#else
#include <sys/time.h>
/**
 \brief Microsecond Timer under Generic Unix
 */
class cTimer {
public:
  cTimer() : t0(time()) { }
  inline void restart() { t0 = time(); }
  inline double elapsed() { return double(time() - t0) * 1.0E-06; }
private:
  uint64_t t0;
  inline uint64_t time() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return 1000000ULL * uint64_t(tv.tv_sec) + uint64_t(tv.tv_usec);
  }
};
#endif

#endif
