/*****************************************************************************
  Template Class for Random Number Generators
  Copyright (c) 2010- Jinglei Hu < jinglei.hu@mpikg.mpg.de >
*****************************************************************************/
#ifndef _PRNG_H_
#define _PRNG_H_

#include <limits.h>
#include <math.h>
#include <new>
#include <stdint.h>
#include <stdlib.h>

/*!
 \brief <B>Usage Example</B>
 \brief Uniform RNG: cPRNG< LFSR113, UNIFORM >
 \brief Normal  RNG: cPRNG< MT19937, NORMAL::ZIGGURAT >

 \brief <B>Speed</B> based on $10^9$ times of generation (icpc -O2 -msse2)
 \brief  Uniform        t(nsec)
 \brief  LFSR113        --
 \brief  LFSR238        --
 \brief  GFSR4          2.23251
 \brief  KISS           --
 \brief  SUPERKISS      2.38473
 \brief  MWC256         3.97905
 \brief  CMWC4096       2.52524
 \brief  XOR4096        2.81664
 \brief  MOTHER         --
 \brief  MT19937        2.24437
 \brief  SFMT19937      1.93612
 \brief  WELL512        5.23513
 \brief  WELL1024a      4.75843
 \brief  WELL19937c     6.8555
 \brief  WELL44497b     7.65181
 \brief  ISAAC8         4.04082
 \brief  ISAAC16        7.28034

 \brief  Ziggurat        t(nsec)
 \brief  LFSR113        --
 \brief  LFSR238        --
 \brief  GFSR4          4.26816
 \brief  KISS           --
 \brief  SUPERKISS      6.8465
 \brief  MWC256         7.33663
 \brief  CMWC4096       5.4318
 \brief  XOR4096        5.56849
 \brief  MOTHER         --
 \brief  MT19937        6.89164
 \brief  SFMT19937      4.02578
 \brief  WELL512        7.68303
 \brief  WELL1024a      7.38846
 \brief  WELL19937c     10.4152
 \brief  WELL44497b     11.1209
 \brief  ISAAC8         6.5721
 \brief  ISAAC16        10.0114

 \warning <B>TEST IT BEFORE YOU USE IT!</B>
 */

/*!
 \brief Ran2 from Numerical Recipe with period 2^61
 \brief Pass BigCrush test in TestU01
 */
class RAN2 {
public:
  static const int BIT = 64;
  typedef uint64_t INT;
  RAN2() : v(4101842887655102017ULL), w(1) { }
  void init(uint64_t seed) {
    u = seed ^ v;
    gen_int64();
    v = u;
    gen_int64();
    w = v;
    gen_int64();
  }
  inline uint64_t gen_int64() {
    u = u * 2862933555777941757ULL + 7046029254386353087ULL;
    v ^= v >> 17;
    v ^= v << 31;
    v ^= v >> 8;
    w = 4294957665U * (w & 0xffffffff) + (w >> 32);
    uint64_t x = u ^ (u << 21);
    x ^= x >> 35;
    x ^= x << 4;
    return (x + v) ^ w;
  }
  inline uint32_t gen_int32() { return uint32_t(gen_int64() >> 32); }
private:
  uint64_t u, v, w;
};

/*!
 \brief Combined LFSR Tausworthe Generators by Pierre L'Ecuyer
 \brief c.f. http://www.iro.umontreal.ca/~simardr/rng/
 */
class LFSR113 {
public:
  static const int BIT = 32;
  typedef int32_t INT;
  LFSR113() : z1(12345), z2(12345), z3(12345), z4(12345) { } //!< Initial seeds: z1 > 1, z2 > 7, z3 > 15, z4 > 128
  void init(int32_t seed) {
    if (seed <= 0) seed = 1;
    int k = seed / 127773;
    seed = 16807 * (seed - k * 127773) - 2836 * k;
    if (seed < 0) seed += 2147483647;
    if (seed < 2)
      z1 = seed + 2;
    else
      z1 = seed;
    k = seed / 127773;
    seed = 16807 * (seed - k * 127773) - 2836 * k;
    if (seed < 0) seed += 2147483647;
    if (seed < 8)
      z2 = seed + 8;
    else
      z2 = seed;
    k = seed / 127773;
    seed = 16807 * (seed - k * 127773) - 2836 * k;
    if (seed < 0) seed += 2147483647;
    if (seed < 16)
      z3 = seed + 16;
    else
      z3 = seed;
    k = seed / 127773;
    seed = 16807 * (seed - k * 127773) - 2836 * k;
    if (seed < 0) seed += 2147483647;
    if (seed < 128)
      z4 = seed + 128;
    else
      z4 = seed;
  }
  inline uint32_t gen_int32() {
    uint32_t b  = ((z1 << 6) ^ z1) >> 13;
    z1 = ((z1 & 4294967294UL) << 18) ^ b;
    b  = ((z2 << 2) ^ z2) >> 27;
    z2 = ((z2 & 4294967288UL) << 2) ^ b;
    b  = ((z3 << 13) ^ z3) >> 21;
    z3 = ((z3 & 4294967280UL) << 7) ^ b;
    b  = ((z4 << 3) ^ z4) >> 12;
    z4 = ((z4 & 4294967168UL) << 13) ^ b;
    return (z1 ^ z2 ^ z3 ^ z4);
  }
  inline uint64_t gen_int64() { return (uint64_t)gen_int32() | (uint64_t)gen_int32() << 32; }
private:
  uint32_t z1, z2, z3, z4;
};

class LFSR238 {
public:
  static const int BIT = 64;
  typedef int32_t INT;
  /*!
   \brief Initial seeds y1 > 1, y2 > 511, y3 > 4095, y4 > 131971, y5 > 8388607
   */
  LFSR238() : y1(987654321), y2(987654321), y3(987654321), y4(987654321), y5(987654321) { }
  void init(int32_t seed) {
    if (seed <= 0) seed = 1;
    int k = seed / 127773;
    seed = 16807 * (seed - k * 127773) - 2836 * k;
    if (seed < 0) seed += 2147483647;
    if (seed < 2)
      y1 = seed + 2;
    else
      y1 = seed;
    k = seed / 127773;
    seed = 16807 * (seed - k * 127773) - 2836 * k;
    if (seed < 0) seed += 2147483647;
    if (seed < 512)
      y2 = seed + 512;
    else
      y2 = seed;
    k = seed / 127773;
    seed = 16807 * (seed - k * 127773) - 2836 * k;
    if (seed < 0) seed += 2147483647;
    if (seed < 4096)
      y3 = seed + 4096;
    else
      y3 = seed;
    k = seed / 127773;
    seed = 16807 * (seed - k * 127773) - 2836 * k;
    if (seed < 0) seed += 2147483647;
    if (seed < 131072)
      y4 = seed + 131072;
    else
      y4 = seed;
    k = seed / 127773;
    seed = 16807 * (seed - k * 127773) - 2836 * k;
    if (seed < 0) seed += 2147483647;
    if (seed < 8388608)
      y5 = seed + 8388608;
    else
      y5 = seed;
  }
  inline uint64_t gen_int64() {
    uint64_t b = ((y1 << 1) ^ y1) >> 53;
    y1 = ((y1 & 18446744073709551614ULL) << 10) ^ b;
    b = ((y2 << 24) ^ y2) >> 50;
    y2 = ((y2 & 18446744073709551104ULL) << 5) ^ b;
    b = ((y3 << 3) ^ y3) >> 23;
    y3 = ((y3 & 18446744073709547520ULL) << 29) ^ b;
    b = ((y4 << 5) ^ y4) >> 24;
    y4 = ((y4 & 18446744073709420544ULL) << 23) ^ b;
    b = ((y5 << 3) ^ y5) >> 33;
    y5 = ((y5 & 18446744073701163008ULL) << 8) ^ b;
    return (y1 ^ y2 ^ y3 ^ y4 ^ y5);
  }
  inline uint32_t gen_int32() { return (uint32_t)(gen_int64() >> 32); }
private:
  uint64_t y1, y2, y3, y4, y5;
};

/*!
 \brief Generalized Feedback Shift-Register Random Number Generator
 \brief Period : 2^9689 - 1
 \brief c.f. Robert M. Ziff, "Four-tap shift-register-sequence random-number generators", Computers in Physics 1998, 12(4):385–392. 
 */
class GFSR4 {
public:
  static const int BIT = 32;
  typedef uint32_t INT;
  GFSR4() { }
  void init(uint32_t seed) {
    if (seed == 0) seed = 4357; // default seed
    // Avoid low-order bit correlations
    for (int i = 0; i < 16384; ++i) {
      uint32_t t = 0, bit = 0x80000000UL;
      for (int j = 0; j < 32; ++j) {
        seed = (69069 * seed) & 0xffffffffUL;
        if (seed & 0x80000000UL) t |= bit;
        bit >>= 1;
      }
      STATE[i] = t;
    }
    uint32_t msb = 0x80000000UL, mask = 0xffffffffUL;
    for (int i = 0; i < 32; ++i) {
      int k = 7 + i * 3;
      STATE[k] &= mask;
      STATE[k] |= msb;
      mask >>= 1;
      msb >>= 1;
    }
    state_i = 32;
  }
  inline uint32_t gen_int32() {
    state_i = (state_i + 1) & 16383UL;
    return STATE[state_i] = STATE[(state_i + 15913UL) & 16383UL] ^ STATE[(state_i + 14798UL) & 16383UL]
                          ^ STATE[(state_i + 9396UL) & 16383UL] ^ STATE[(state_i + 6695UL) & 16383UL];
  }
  inline uint64_t gen_int64() { return (uint64_t)gen_int32() | (uint64_t)gen_int32() << 32; }
private:
  int32_t  state_i;
  uint32_t STATE[16384];
};

/*!
 \brief George Marsaglia's PRNGs
 \brief KISS passes all of the Dieharder tests and the BigCrunch tests in TestU01.
 \brief c.f. http://www.cs.ucl.ac.uk/staff/d.jones/GoodPracticeRNG.pdf
 */
class KISS {
public:
  static const int BIT = 32;
  typedef uint32_t INT;
  KISS() : x(123456789), y(362436069), z(21288629), w(21288629), c(0) {  }
  void init(uint32_t seed) {
    x = seed & 0xffffffffUL;
    y = (1812433253UL * (x ^ (x >> 30)) + 1) & 0xffffffffUL;
    z = (1812433253UL * (y ^ (y >> 30)) + 2) & 0xffffffffUL;
    w = (1812433253UL * (z ^ (z >> 30)) + 3) & 0xffffffffUL;
    c = w % 698769068 + 1;
  }
  uint32_t gen_int32() {
    y ^= (y << 5);
    y ^= (y >> 7);
    y ^= (y << 22);
    int32_t t = z + w + c;
    z = w;
    c = t < 0;
    w = t & 2147483647;
    x += 1411392427;
    return (x + y + w);
  }
  inline uint64_t gen_int64() { return (uint64_t)gen_int32() | (uint64_t)gen_int32() << 32; }
private:
  uint32_t x, y, z, w, c;
};

/*!
 \brief c.f. http://www.velocityreviews.com/forums/t704080-re-rngs-a-super-kiss.html
 \brief Period 54767*2^1337279
 */
class SUPERKISS {
public:
  static const int BIT = 32;
  typedef int32_t INT;
  SUPERKISS() : m_index(41790), m_carry(362436), xcng(1236789), xs(521288629) { }
  void init(int32_t seed) {
    m_carry = uint32_t(seed & 0xffffffffUL) % 7010176;
    xcng = (1812433253UL * (m_carry ^ (m_carry >> 30)) + 1) & 0xffffffffUL;
    xs = (1812433253UL * (xcng ^ (xcng >> 30)) + 2) & 0xffffffffUL;
    for (int i = 0; i < 41790; ++i)
      m_q[i] = (xcng = 69609 * xcng + 123) + (xs ^= xs << 13, xs ^= (unsigned) xs >> 17, xs ^= xs >> 5);
  }
  inline uint32_t gen_int32() {
    return (m_index < 41790 ? m_q[m_index++] : refill()) + (xcng = 69609 * xcng + 123) +
            (xs ^= xs << 13, xs ^= (unsigned) xs >> 17, xs ^= xs >> 5);
  }
  inline uint64_t gen_int64() { return (uint64_t)gen_int32() | (uint64_t)gen_int32() << 32; }
private:
  uint32_t m_index, m_carry, xcng, xs;
  uint32_t m_q[41790];
  inline uint32_t refill() {
    for (int32_t i = 0; i < 41790; ++i) {
      uint64_t t = 7010176ULL * m_q[i] + m_carry;
      m_carry = t >> 32;
      m_q[i] = ~t;
    }
    m_index = 1;
    return m_q[0];
  }
};

/*!
 \brief Period 2^8222
 \brief c.f. http://forums.wolfram.com/mathgroup/archive/2003/Feb/msg00456.html
 */
class MWC256 {
public:
  static const int BIT = 32;
  typedef uint32_t INT;
  MWC256() { }
  void init(uint32_t seed) {
    m_q[0] = seed;
    for (uint32_t i = 1; i < N; ++i)
      m_q[i] = 1812433253UL * (m_q[i-1] ^ (m_q[i-1] >> 30)) + i;
    m_carry = m_q[N-1] % 61137367UL;
    m_index = N - 1;
  }
  inline uint32_t gen_int32() {
    ++m_index;
    uint64_t temp = 1540315826ULL * m_q[m_index] + m_carry;
    m_carry = temp >> 32;
    uint32_t x = temp + m_carry;
    if (x < m_carry) {
      ++x;
      ++m_carry;
    }
    return m_q[m_index] = x;
  }
  inline uint64_t gen_int64() { return (uint64_t)gen_int32() | (uint64_t)gen_int32() << 32; }
private:
  enum { N = 256 };
  uint8_t  m_index;
  uint32_t m_carry;
  uint32_t m_q[N];
};

/*!
 \brief Period 2^131086
 */
class CMWC4096 {
public:
  static const int BIT = 32;
  typedef uint32_t INT;
  CMWC4096() { }
  void init(uint32_t seed) {
    m_q[0] = seed;  // Save seed for historical purpose
    for (uint32_t i = 1; i < N; ++i) // Set the array using one of Knuth's generators
      m_q[i] = 1812433253UL * (m_q[i - 1] ^ (m_q[i - 1] >> 30)) + i;
    m_carry = m_q[N-1] % 61137367UL;
    m_index = N - 1;
  }
  inline uint32_t gen_int32() {
    m_index = (m_index + 1) & 4095;
    uint64_t temp = 18782ULL * m_q[m_index]+ m_carry;
    m_carry = temp >> 32;
    uint32_t x = temp + m_carry;
    if (x < m_carry) {
      ++x;
      ++m_carry;
    }
    return m_q[m_index] = 0xfffffffeUL - x; 
  }
  inline uint64_t gen_int64() { return (uint64_t)gen_int32() | (uint64_t)gen_int32() << 32; }
private:
  enum { N = 4096 };
  uint32_t m_index;
  uint32_t m_carry;
  uint32_t m_q[N];
};

/*!
 \brief Adapted from Anger Fog's implementation
 */
class MOTHER {
public:
  static const int BIT = 32;
  typedef uint32_t INT;
  void init(uint32_t seed) {
    uint32_t s = seed;
    for (int32_t i = 0; i < 5; ++i) {
      s = s * 29943829UL - 1;
      x[i] = s;
    }
  }
  inline uint32_t gen_int32() {
    uint64_t sum = 2111111111ULL * uint64_t(x[3]) + 1492ULL * uint64_t(x[2])
                 + 1776ULL * uint64_t(x[1]) + 5115ULL * uint64_t(x[0]) + uint64_t(x[4]);
    x[3] = x[2];
    x[2] = x[1];
    x[1] = x[0];
    x[4] = uint32_t(sum >> 32);  // Carry
    x[0] = (uint32_t)sum;        // Low 32 bits of sum
    return x[0];
  }
  inline uint64_t gen_int64() { return (uint64_t)gen_int32() | (uint64_t)gen_int32() << 32; }
private:
  uint32_t x[5];
};

/*!
 \brief Brent's XOR4096
 \brief c.f. ftp://ftp.comlab.ox.ac.uk/pub/Documents/techpapers/Richard.Brent/random/xorgens201.c
 */
class XOR4096 {
public:
  static const int BIT = 32;
  typedef uint32_t INT;
  XOR4096() { }
  void init(uint32_t seed) {
    uint32_t v = seed != 0 ? seed : ~seed;
    uint32_t k;
    for (k = 32; k > 0; --k) {
      v ^= v << 13;
      v ^= v >> 17;
      v ^= v << 5;
    } // v ^= (v ^= (v ^= v << 13) >> 17) << 5;
    for (w = v, k = 0; k < 128; ++k) {
      v ^= v << 13;
      v ^= v >> 17;
      v ^= v << 5;
      w += 0x61c88647UL;
      STATE[k] = v + w;
    } // STATE[k] = (v ^= (v ^= (v ^= v << 13) >> 17) << 5) + (w += 0x61c88647UL);
    uint32_t t;
    for (state_i = 127, k = 512; k > 0; --k) {
      t = STATE[state_i = (state_i + 1) & 127];
      v = STATE[(state_i + 33) & 127];
      t ^= t << 17;
      t ^= t >> 12;
      v ^= v << 13;
      STATE[state_i] = (v ^= t ^ (v >> 15));
    }
  }
  inline uint32_t gen_int32() {
    uint32_t t = STATE[state_i = (state_i + 1) & 127];
    uint32_t v = STATE[(state_i + 33) & 127];
    t ^= t << 17;
    t ^= t >> 12;
    v ^= v << 13;
    STATE[state_i] = (v ^= t ^ (v >> 15));
    return v + (w += 0x61c88647UL);
  }
  inline uint64_t gen_int64() { return (uint64_t)gen_int32() | (uint64_t)gen_int32() << 32; }
private:
  uint32_t state_i, w, STATE[128];
};

/*!
 \brief Mersenne Twister PRNGs
 */
#if INT_MAX == 2147483647
class MT19937 {
public:
  static const int BIT = 32;
  typedef uint32_t INT;
  MT19937() : mti(N + 1) { }
  void init(uint32_t seed) {
    mt[0]= seed & 0xffffffffUL;
    for (mti = 1; mti < N; ++mti) {
      mt[mti] = (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);
      mt[mti] &= 0xffffffffUL;
    }
  }
  /*!
   \brief initialize by an array with array-length
   @param init_key the array for initializing keys
   @param key_length its length
   */
  void init(uint32_t init_key[], uint32_t key_length) {
    uint32_t i = 1, j = 0, k = N > key_length ? N : key_length;
    init(19650218UL);
    for (; k > 0; --k) {
      mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL)) + init_key[j] + j; // non linear
      mt[i] &= 0xffffffffUL;
      ++i;
      ++j;
      if (i >= N) {
        mt[0] = mt[N-1];
        i = 1;
      }
      if (j >= key_length) j = 0;
    }
    for (k = N - 1; k > 0; --k) {
      mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL)) - i; // non linear
      mt[i] &= 0xffffffffUL;
      ++i;
      if (i >= N) {
        mt[0] = mt[N-1];
        i = 1;
      }
    }
    mt[0] = 0x80000000UL; // MSB is 1; assuring non-zero initial array
  }
  inline uint32_t gen_int32() {
    uint32_t y;
    if (mti >= N) { // generate N words at one time
      static const uint32_t mag01[2] = { 0x0UL, 0x9908b0dfUL };
      int32_t kk;
      if (mti == N + 1)   // if init() has not been called,
        init(5489UL); // a default initial seed is used
      for (kk = 0; kk < N - M; ++kk) {
        y = (mt[kk] & 0x80000000UL)|(mt[kk+1] & 0x7fffffffUL);
        mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
      }
      for (; kk < N - 1; ++kk) {
        y = (mt[kk] & 0x80000000UL)|(mt[kk+1] & 0x7fffffffUL);
        mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
      }
      y = (mt[N-1] & 0x80000000UL)|(mt[0] & 0x7fffffffUL);
      mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];
      mti = 0;
    }
    y = mt[mti++];
    // Tempering
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);
    return y;
  } //!< random number in [0,0xffffffff]-interval
  inline uint64_t gen_int64() { return (uint64_t)gen_int32() | (uint64_t)gen_int32() << 32; }
private:
  enum { N = 624, M = 397 };
  uint32_t mti;
  uint32_t mt[N];
};
#else
class MT19937 {
public:
  static const int BIT = 64;
  typedef uint64_t INT;
  MT19937() : mti(N + 1) { }
  void init(uint64_t seed) {
    mt[0] = seed;
    for (mti = 1; mti < N; ++mti)
      mt[mti] = (6364136223846793005ULL * (mt[mti-1] ^ (mt[mti-1] >> 62)) + mti);
  }
  void init(uint64_t init_key[], uint64_t key_length) {
    uint64_t i = 1, j = 0, k = N > key_length ? N : key_length;
    init(19650218ULL);
    for (; k > 0; --k) {
      mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 62)) * 3935559000370003845ULL)) + init_key[j] + j;
      ++i;
      ++j;
      if (i >= N) {
        mt[0] = mt[N-1];
        i = 1;
      }
      if (j >= key_length) j = 0;
    }
    for (k = N - 1; k > 0; --k) {
      mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 62)) * 2862933555777941757ULL)) - i;
      ++i;
      if (i >= N) {
        mt[0] = mt[N-1];
        i = 1;
      }
    }
    mt[0] = 1ULL << 63;
  }
  inline uint64_t gen_int64() {
    uint64_t x;
    if (mti >= N) {
      static uint64_t mag01[2] = { 0ULL, 0xB5026F5AA96619E9ULL };
      // generate N words at one time
      // if init() has not been called,
      // a default initial seed is used
      if (mti == N + 1) init(5489ULL);
      int32_t i;
      for (i = 0; i < N - M; ++i) {
        x = (mt[i] & 0xFFFFFFFF80000000ULL)|(mt[i+1] & 0x7FFFFFFFULL);
        mt[i] = mt[i+M] ^ (x>>1) ^ mag01[(int)(x&1ULL)];
      }
      for (; i < N - 1; ++i) {
        x = (mt[i] & 0xFFFFFFFF80000000ULL)|(mt[i+1] & 0x7FFFFFFFULL);
        mt[i] = mt[i+(M-N)] ^ (x>>1) ^ mag01[(int)(x&1ULL)];
      }
      x = (mt[N-1] & 0xFFFFFFFF80000000ULL)|(mt[0] & 0x7FFFFFFFULL);
      mt[N-1] = mt[M-1] ^ (x>>1) ^ mag01[(int)(x&1ULL)];
      mti = 0;
    }
    x = mt[mti++];
    x ^= (x >> 29) & 0x5555555555555555ULL;
    x ^= (x << 17) & 0x71D67FFFEDA60000ULL;
    x ^= (x << 37) & 0xFFF7EEE000000000ULL;
    x ^= (x >> 43);
    return x;
  } //!< random number on [0, 2^64-1]-interval
  inline uint32_t gen_int32() { return uint32_t(gen_int64() >> 32); }
private:
  enum { N = 312, M = 156 };
  uint32_t mti;
  uint64_t mt[N];
};
#endif

#ifdef __SSE2__
#include <emmintrin.h>
template<uint32_t MEXP, uint32_t POS1, int SL1, int SL2, int SR1, int SR2, uint32_t MSK1, uint32_t MSK2, uint32_t MSK3, uint32_t MSK4, uint32_t PARITY1, uint32_t PARITY2, uint32_t PARITY3, uint32_t PARITY4>
class cSFMT {
public:
  static const int BIT = 64;
  typedef uint32_t INT;
  cSFMT() {
    psfmt32 = &sfmt[0].u[0];
    psfmt64 = (uint64_t*)&sfmt[0].u[0];
    initialized = 0;
    parity[0] = PARITY1;
    parity[1] = PARITY2;
    parity[2] = PARITY3;
    parity[3] = PARITY4;
  }
  void init(uint32_t seed) {
    psfmt32[0] = seed;
    for (uint32_t i = 1; i < N32; i++)
      psfmt32[i] = 1812433253UL * (psfmt32[i - 1] ^ (psfmt32[i-1] >> 30)) + i;
    idx = N32;
    period_certification();
    initialized = 1;
  }
  void init(uint32_t init_key[], int key_length) {
    uint32_t lag, size = N * 4;
    if (size >= 623)
      lag = 11;
    else if (size >= 68)
      lag = 7;
    else if (size >= 39)
      lag = 5;
    else
      lag = 3;
    uint32_t mid = (size - lag) / 2;
    uint32_t count = key_length + 1 > N32 ? key_length + 1 : N32;
    uint32_t r = func1(psfmt32[0] ^ psfmt32[mid] ^ psfmt32[N32-1]);
    psfmt32[mid] += r;
    r += key_length;
    psfmt32[mid+lag] += r;
    psfmt32[0] = r;
    --count;
    int32_t i, j;
    for (i = 1, j = 0; j < count && j < key_length; ++j) {
      r = func1(psfmt32[i] ^ psfmt32[(i+mid)%N32] ^ psfmt32[(i+N32-1)%N32]);
      psfmt32[(i+mid)%N32] += r;
      r += init_key[j] + i;
      psfmt32[(i+mid+lag)%N32] += r;
      psfmt32[i] = r;
      i = (i + 1) % N32;
    }
    for (; j < count; ++j) {
      r = func1(psfmt32[i] ^ psfmt32[(i+mid)%N32] ^ psfmt32[(i+N32-1)%N32]);
      psfmt32[(i+mid)%N32] += r;
      r += i;
      psfmt32[(i+mid+lag)%N32] += r;
      psfmt32[i] = r;
      i = (i + 1) % N32;
    }
    for (j = 0; j < N32; j++) {
      r = func2(psfmt32[i] + psfmt32[(i+mid)%N32] + psfmt32[(i+N32-1)%N32]);
      psfmt32[(i+mid)%N32] ^= r;
      r -= i;
      psfmt32[(i+mid+lag)%N32] ^= r;
      psfmt32[i] = r;
      i = (i + 1) % N32;
    }
    idx = N32;
    period_certification();
    initialized = 1;
  }
  inline uint32_t gen_int32() {
    if (idx >= N32) {
      gen_rand_all();
      idx = 0;
    }
    return psfmt32[idx++];
  }
  inline uint64_t gen_int64() {
    if (idx >= N32) {
      gen_rand_all();
      idx = 0;
    }
    uint64_t r = psfmt64[idx/2];
    idx += 2;
    return r;
  }
private:
  enum { N = MEXP / 128 + 1, N64 = N + N, N32 = N64 + N64 };
  int32_t  idx;
  int32_t  initialized;
  static inline __m128i mm_recursion(__m128i* a, const __m128i* b, __m128i c, __m128i d, __m128i mask) {
    __m128i v, x, y, z;
    x = _mm_load_si128(a);
    y = _mm_srli_epi32(*b, SR1);
    z = _mm_srli_si128(c, SR2);
    v = _mm_slli_epi32(d, SL1);
    z = _mm_xor_si128(z, x);
    z = _mm_xor_si128(z, v);
    x = _mm_slli_si128(x, SL2);
    y = _mm_and_si128(y, mask);
    z = _mm_xor_si128(z, x);
    z = _mm_xor_si128(z, y);
    return z;
  }
  union w128_t {
    __m128i si;
    uint32_t u[4];
  };
  w128_t   sfmt[N];
  uint32_t parity[4];
  uint32_t *psfmt32;
  uint64_t *psfmt64;
  static uint32_t func1(uint32_t x) { return (x ^ (x >> 27)) * (uint32_t)1664525UL; }
  static uint32_t func2(uint32_t x) { return (x ^ (x >> 27)) * (uint32_t)1566083941UL; }
  void period_certification() {
    int32_t inner = 0;
    for (uint32_t i = 0; i < 4; ++i)
      inner ^= psfmt32[i] & parity[i];
    for (uint32_t i = 16; i > 0; i >>= 1)
      inner ^= inner >> i;
    inner &= 1;
    // check OK
    if (inner == 1) return;
    // check NG, and modification
    for (uint32_t i = 0; i < 4; ++i) {
      uint32_t work = 1;
      for (uint32_t j = 0; j < 32; ++j) {
        if ((work & parity[i]) != 0) {
          psfmt32[i] ^= work;
          return;
        }
        work <<= 1;
      }
    }
  }
  void gen_rand_all() {
    __m128i r, r1, r2;
    static const __m128i mask = _mm_set_epi32( MSK4, MSK3, MSK2, MSK1);
    r1 = _mm_load_si128(&sfmt[N - 2].si);
    r2 = _mm_load_si128(&sfmt[N - 1].si);
    uint32_t i;
    for (i = 0; i < N - POS1; ++i) {
      r = mm_recursion(&sfmt[i].si, &sfmt[i + POS1].si, r1, r2, mask);
      _mm_store_si128(&sfmt[i].si, r);
      r1 = r2;
      r2 = r;
    }
    for (; i < N; ++i) {
      r = mm_recursion(&sfmt[i].si, &sfmt[i + POS1 - N].si, r1, r2, mask);
      _mm_store_si128(&sfmt[i].si, r);
      r1 = r2;
      r2 = r;
    }
  }
  void gen_rand_array(w128_t *array, uint32_t size) {
    __m128i r, r1, r2;
    static const __m128i mask = _mm_set_epi32(MSK4, MSK3, MSK2, MSK1);
    r1 = _mm_load_si128(&sfmt[N - 2].si);
    r2 = _mm_load_si128(&sfmt[N - 1].si);
    int32_t i, j;
    for (i = 0; i < N - POS1; ++i) {
      r = mm_recursion(&sfmt[i].si, &sfmt[i + POS1].si, r1, r2, mask);
      _mm_store_si128(&array[i].si, r);
      r1 = r2;
      r2 = r;
    }
    for (; i < N; ++i) {
      r = mm_recursion(&sfmt[i].si, &array[i + POS1 - N].si, r1, r2, mask);
      _mm_store_si128(&array[i].si, r);
      r1 = r2;
      r2 = r;
    }
    for (; i < size - N; ++i) {
      r = mm_recursion(&array[i - N].si, &array[i + POS1 - N].si, r1, r2, mask);
      _mm_store_si128(&array[i].si, r);
      r1 = r2;
      r2 = r;
    }
    for (j = 0; j < 2 * N - size; ++j) {
      r = _mm_load_si128(&array[j + size - N].si);
      _mm_store_si128(&sfmt[j].si, r);
    }
    for (; i < size; ++i) {
      r = mm_recursion(&array[i - N].si, &array[i + POS1 - N].si, r1, r2, mask);
      _mm_store_si128(&array[i].si, r);
      _mm_store_si128(&sfmt[j++].si, r);
      r1 = r2;
      r2 = r;
    }
  }
};

/*!
 \typedef Alias for SIMD-oriented Fast Mersenne Twister Class
 */
typedef cSFMT<607, 2, 15, 3, 13, 3, 0xfdff37ffU, 0xef7f3f7dU, 0xff777b7dU, 0x7ff7fb2fU,
0x00000001U, 0x00000000U, 0x00000000U, 0x5986f054U> SFMT607;
typedef cSFMT<1279, 7, 14, 3, 5, 1, 0xf7fefffdU, 0x7fefcfffU, 0xaff3ef3fU, 0xb5ffff7fU,
0x00000001U, 0x00000000U, 0x00000000U, 0x20000000U> SFMT1279;
typedef cSFMT<2281, 12, 19, 1, 5, 1, 0xbff7ffbfU, 0xfdfffffeU, 0xf7ffef7fU, 0xf2f7cbbfU,
0x00000001U, 0x00000000U, 0x00000000U, 0x41dfa600U> SFMT2281;
typedef cSFMT<4253, 17, 20, 1, 7, 1, 0x9f7bffff, 0x9fffff5f, 0x3efffffb, 0xfffff7bb,
0xa8000001U, 0xaf5390a3U, 0xb740b3f8U, 0x6c11486dU> SFMT4253;
typedef cSFMT<12213, 68, 14, 3, 7, 3, 0xeffff7fbU, 0xffffffefU, 0xdfdfbfffU, 0x7fffdbfdU,
0x00000001U, 0x00000000U, 0xe8148000U, 0xd0c7afa3U> SFMT12213;
typedef cSFMT<19937, 122, 18, 1, 11, 1, 0xdfffffefU, 0xddfecb7fU, 0xbffaffffU, 0xbffffff6U,
0x00000001U, 0x00000000U, 0x00000000U, 0x13c9e684U> SFMT19937;
typedef cSFMT<44497, 330, 5, 3, 9, 3, 0xeffffffbU, 0xdfbebfffU, 0xbfbf7befU, 0x9ffd7bffU,
0x00000001U, 0x00000000U, 0xa3ac4000U, 0xecc1327aU> SFMT44497;
typedef cSFMT<86243, 366, 6, 7,19, 1, 0xfdbffbffU, 0xbff7ff3fU, 0xfd77efffU, 0xbf9ff3ffU, 
0x00000001U, 0x00000000U, 0x00000000U, 0xe9528d85U> SFMT86243;
typedef cSFMT<132049, 110, 19, 1, 21, 1, 0xffffbb5fU, 0xfb6ebf95U, 0xfffefffaU, 0xcff77fffU,
0x00000001U, 0x00000000U, 0xcb520000U, 0xc7e91c7dU> SFMT132049;
typedef cSFMT<216091, 627, 11, 3, 10, 1, 0xbff7bff7U, 0xbfffffffU, 0xbffffa7fU, 0xffddfbfbU, 
0xf8000001U, 0x89e80709U, 0x3bd2b64bU, 0x0c64b1e4U> SFMT216091;
#endif

/*!
 \brief Well Equidistributed Long-period Linear (WELL) Random Number Generators
 \brief Adapted from the C code in http://www.lomont.org/Math/Papers/2008/Lomont_PRNG_2008.pdf
 \brief which is 40% faster than the original implementation
 */
class WELL512 {
public:
  static const int BIT = 32;
  typedef uint32_t INT;
  WELL512() : state_i(0) { }
  void init(uint32_t seed) {
    STATE[0]= seed & 0xffffffffUL;
    for (uint32_t i = 1; i < 16; ++i) {
      STATE[i] = (1812433253UL * (STATE[i-1] ^ (STATE[i-1] >> 30)) + i); 
      STATE[i] &= 0xffffffffUL;
    }
  }
  inline uint32_t gen_int32() {
    uint32_t a = STATE[state_i];
    uint32_t c = STATE[(state_i + 13) & 15];
    uint32_t b = a ^ c ^ (a << 16) ^ (c << 15);
    c = STATE[(state_i + 9) & 15];
    c ^= (c >> 11);
    a = STATE[state_i] = b ^ c;
    uint32_t d = a ^ ((a << 5) & 0xDA442D20UL);
    state_i = (state_i + 15) & 15;
    a = STATE[state_i];
    STATE[state_i] = a ^ b ^ d ^ (a << 2) ^ (b << 18) ^ (c << 28);
    return STATE[state_i];
  }
  inline uint64_t gen_int64() { return (uint64_t)gen_int32() | (uint64_t)gen_int32() << 32; }
private:
  uint32_t state_i;
  uint32_t STATE[16];
};

/*!
 \brief Adapted from C++ code at http://www.mpipks-dresden.mpg.de/~gross/download/myrng/myrng1.2/
 */
class WELL1024a {
public:
  static const int BIT = 32;
  typedef int32_t INT;
  WELL1024a() : state_i(0) { }
  void init(int32_t seed) {
    STATE[0]= seed & 0xffffffffUL;
    for (uint32_t i = 1; i < 32; ++i) {
      STATE[i] = (1812433253UL * (STATE[i-1] ^ (STATE[i-1] >> 30)) + i); 
      STATE[i] &= 0xffffffffUL;
    }
  }
  inline uint32_t gen_int32() {
    uint32_t z0 = STATE[(state_i + 31) & 0x0000001fU];
    uint32_t z1 = (STATE[state_i])^ STATE[(state_i + 3) & 0x0000001fU];
    uint32_t z2 = (STATE[(state_i + 24) & 0x0000001fU] ^ (STATE[(state_i + 24) & 0x0000001fU] << 19))
                ^ (STATE[(state_i + 10) & 0x0000001fU] ^ (STATE[(state_i + 10) & 0x0000001fU] << 14));
    STATE[state_i] = z1 ^ z2;
    STATE[(state_i + 31) & 0x0000001fU] = (z0 ^ (z0 << 11)) ^ (z1 ^ (z1 << 7)) ^ (z2 ^ (z2 << 13));
    state_i = (state_i + 31) & 0x0000001fUL;
    return STATE[state_i];
  }
  inline uint64_t gen_int64() { return (uint64_t)gen_int32() | (uint64_t)gen_int32() << 32; }
private:
  uint32_t state_i;
  uint32_t STATE[32];
};

/*!
 \brief Adapted from Java implementation which can be found at 
 \brief https://issues.apache.org/jira/browse/MATH?report=com.atlassian.jira.plugin.ext.subversion:subversion-project-tab
 */
class WELL19937c {
public:
  static const int BIT = 32;
  typedef int32_t INT;
  WELL19937c() : state_i(0) {
    for (uint32_t i = 0; i < N; ++i) {
      iRm1[i] = (i + N - 1) % N;
      iRm2[i] = (i + N - 2) % N;
      i1[i] = (i + M1) % N;
      i2[i] = (i + M2) % N;
      i3[i] = (i + M3) % N;
    }
  }
  void init(int32_t seed) {
    STATE[0]= seed & 0xffffffffUL;
    for (uint32_t i = 1; i < 32; ++i) {
      STATE[i] = (1812433253UL * (STATE[i-1] ^ (STATE[i-1] >> 30)) + i); 
      STATE[i] &= 0xffffffffUL;
    }
  }
  inline uint32_t gen_int32() {
    int indexRm1 = iRm1[state_i];
    int indexRm2 = iRm2[state_i];
    int v0       = STATE[state_i];
    int vM1      = STATE[i1[state_i]];
    int vM2      = STATE[i2[state_i]];
    int vM3      = STATE[i3[state_i]];
    int z0 = (0x80000000 & STATE[indexRm1]) ^ (0x7FFFFFFF & STATE[indexRm2]);
    int z1 = (v0 ^ (v0 << 25))  ^ (vM1 ^ (vM1 >> 27));
    int z2 = (vM2 >> 9) ^ (vM3 ^ (vM3 >> 1));
    int z3 = z1 ^ z2;
    int z4 = z0 ^ (z1 ^ (z1 << 9)) ^ (z2 ^ (z2 << 21)) ^ (z3 ^ (z3 >> 21));
    STATE[state_i] = z3;
    STATE[indexRm1] = z4;
    STATE[indexRm2] &= 0x80000000;
    state_i = indexRm1;
    // Add Matsumoto-Kurita tempering to get a maximally-equidistributed generator
    z4 = z4 ^ ((z4 <<  7) & 0xe46e1700);
    z4 = z4 ^ ((z4 << 15) & 0x9b868000);
    return z4;
  }
  inline uint64_t gen_int64() { return (uint64_t)gen_int32() | (uint64_t)gen_int32() << 32; }
private:
  enum { K = 19937, W = 32, N = (K + W - 1) / 32, M1 = 70, M2 = 179, M3 = 449 };
  int state_i;
  uint32_t STATE[N];
  int iRm1[N], iRm2[N];
  int i1[N], i2[N], i3[N];
};

/*!
 \brief Adapted from Java implementation which can be found at 
 \brief https://issues.apache.org/jira/browse/MATH?report=com.atlassian.jira.plugin.ext.subversion:subversion-project-tab
 */
class WELL44497b {
public:
  static const int BIT = 32;
  typedef int32_t INT;
  WELL44497b() : state_i(0) {
    for (uint32_t i = 0; i < N; ++i) {
      iRm1[i] = (i + N - 1) % N;
      iRm2[i] = (i + N - 2) % N;
      i1[i] = (i + M1) % N;
      i2[i] = (i + M2) % N;
      i3[i] = (i + M3) % N;    
    }
  }
  void init(int32_t seed) {
    STATE[0]= seed & 0xffffffffUL;
    for (uint32_t i = 1; i < W; ++i) {
      STATE[i] = (1812433253UL * (STATE[i-1] ^ (STATE[i-1] >> 30)) + i); 
      STATE[i] &= 0xffffffffUL;
    }
  }
  inline uint32_t gen_int32() {
    int indexRm1 = iRm1[state_i];
    int indexRm2 = iRm2[state_i];
    int v0       = STATE[state_i];
    int vM1      = STATE[i1[state_i]];
    int vM2      = STATE[i2[state_i]];
    int vM3      = STATE[i3[state_i]];
    int z0       = (0xFFFF8000 & STATE[indexRm1]) ^ (0x00007FFF & STATE[indexRm2]);
    int z1       = (v0 ^ (v0 << 24))  ^ (vM1 ^ (vM1 >> 30));
    int z2       = (vM2 ^ (vM2 << 10)) ^ (vM3 << 26);
    int z3       = z1 ^ z2;
    int z2Prime  = ((z2 << 9) ^ (z2 >> 23)) & 0xfbffffff;
    int z2Second = ((z2 & 0x00020000) != 0) ? (z2Prime ^ 0xb729fcec) : z2Prime;
    int z4       = z0 ^ (z1 ^ (z1 >> 20)) ^ z2Second ^ z3;
    STATE[state_i]   = z3;
    STATE[indexRm1]  = z4;
    STATE[indexRm2] &= 0xFFFF8000;
    state_i        = indexRm1;
    // Add Matsumoto-Kurita tempering to get a maximally-equidistributed generator
    z4 = z4 ^ ((z4 <<  7) & 0x93dd1400);
    z4 = z4 ^ ((z4 << 15) & 0xfa118000);
    return z4;
  }
  inline uint64_t gen_int64() { return (uint64_t)gen_int32() | (uint64_t)gen_int32() << 32; }
private:
  enum { K = 44497, W = 32, N = (K + W - 1) / 32, M1 = 23, M2 = 481, M3 = 229 };
  int state_i;
  int STATE[N];
  int iRm1[N], iRm2[N];
  int i1[N], i2[N], i3[N];
};

/**
 \brief ISAAC (Indirection, Shift, Accumulate, Add, and Count), Cryptographic Pseudo Random Number Generator
 \brief Expected period: 2 ^ ((ALPHA + BIT * (3 + 2^ALPHA) - 1)
 \brief 2^8295 for ALPHA = 8, BIT = 32
 \brief c.f. http://burtleburtle.net/bob/rand/isaacafa.html
 */
template <int ALPHA>
class cISAAC {
public:
  #if INT_MAX == 2147483647
  static const int BIT = 32;
  typedef uint32_t INT;
  #else
  static const int BIT = 64;
  typedef uint64_t INT;
  #endif
  typedef unsigned char BYTE;
  cISAAC() { INIT(false); }
  void init(INT seed) {
    for(uint32_t i = 0; i < N; ++i)
      rsl[i] = 0;
    #if INT_MAX == 2147483647
    z0 = seed & 0xffffffffUL;
    z1 = (1812433253UL * (z0 ^ (z0 >> 30)) + 1) & 0xffffffffUL;
    z2 = (1812433253UL * (z0 ^ (z0 >> 30)) + 2) & 0xffffffffUL;
    #else
    z0 = seed;
    z1 = (6364136223846793005ULL * (z0 ^ (z0 >> 62)) + 1);
    z2 = (6364136223846793005ULL * (z1 ^ (z1 >> 62)) + 2);
    #endif
    INIT(true);
  }
  #if INT_MAX == 2147483647
  inline uint32_t gen_int32() { return (!K-- ? (isaac(), K = N - 1 , rsl[K]) : rsl[K]); }
  inline uint64_t gen_int64() { return (uint64_t)gen_int32() | (uint64_t)gen_int32() << 32; }
  #else
  inline uint64_t gen_int64() { return (!K-- ? (isaac(), K = N - 1 , rsl[K]) : rsl[K]); }
  inline uint32_t gen_int32() { return uint32_t(gen_int64() >> 32); }
  #endif
private:
  enum { N = 1 << ALPHA };
  INT z0, z1, z2, K;
  INT rsl[N], mem[N];
  inline void isaac() {
    INT *mm = mem, *r = rsl;
    INT a = z0, b = z1 + (++z2);
    INT *m = mm, *m2 = m + N / 2, *mend = m2;
    INT x, y;
    for(; m < mend; ) {
      #if INT_MAX == 2147483647
      rngstep((a << 13), a, b, mm, m, m2, r, x, y);
      rngstep((a >> 6) , a, b, mm, m, m2, r, x, y);
      rngstep((a << 2) , a, b, mm, m, m2, r, x, y);
      rngstep((a >> 16), a, b, mm, m, m2, r, x, y);
      #else
      rngstep(~(a^(a << 21)), a, b, mm, m, m2, r, x, y);
      rngstep(  a^(a >> 5)  , a, b, mm, m, m2, r, x, y);
      rngstep(  a^(a << 12) , a, b, mm, m, m2, r, x, y);
      rngstep(  a^(a >> 33) , a, b, mm, m, m2, r, x, y);
      #endif
    }
    m2 = mm;
    for(; m2 < mend; ) {
      #if INT_MAX == 2147483647
      rngstep((a << 13), a, b, mm, m, m2, r, x, y);
      rngstep((a >> 6) , a, b, mm, m, m2, r, x, y);
      rngstep((a << 2) , a, b, mm, m, m2, r, x, y);
      rngstep((a >> 16), a, b, mm, m, m2, r, x, y);
      #else
      rngstep(~(a^(a << 21)), a, b, mm, m, m2, r, x, y);
      rngstep(  a^(a >> 5)  , a, b, mm, m, m2, r, x, y);
      rngstep(  a^(a << 12) , a, b, mm, m, m2, r, x, y);
      rngstep(  a^(a >> 33) , a, b, mm, m, m2, r, x, y);
      #endif
    }
    z1 = b;
    z0 = a;
  }
  inline INT ind(INT* mm, INT x) {
    #if INT_MAX == 2147483647
    return (*(INT*)((BYTE*)(mm) + ((x) & ((N-1)  << 2))));
    #else 
    return (*(INT*)((BYTE*)(mm) + ((x) & ((N-1)  << 3))));
    #endif
  }
  inline void rngstep(INT mix, INT& a, INT& b, INT*& mm, INT*& m, INT*& m2, INT*& r, INT& x, INT& y) {
    x = *m;
    a = (a^(mix)) + *(m2++);
    *(m++) = y = ind(mm, x) + a + b;
    *(r++) = b = ind(mm, y  >>  ALPHA) + x;
  }
  inline void shuffle(INT& a, INT& b, INT& c, INT& d, INT& e, INT& f, INT& g, INT& h) {
    #if INT_MAX == 2147483647
    a ^= b << 11; d += a; b += c;
    b ^= c >> 2;  e += b; c += d;
    c ^= d << 8;  f += c; d += e;
    d ^= e >> 16; g += d; e += f;
    e ^= f << 10; h += e; f += g;
    f ^= g >> 4;  a += f; g += h;
    g ^= h << 8;  b += g; h += a;
    h ^= a >> 9;  c += h; a += b;
    #else
    a -= e; f ^= h >> 9;  h += a;
    b -= f; g ^= a << 9;  a += b;
    c -= g; h ^= b >> 23; b += c;
    d -= h; a ^= c << 15; c += d;
    e -= a; b ^= d >> 14; d += e;
    f -= b; c ^= e << 20; e += f;
    g -= c; d ^= f >> 17; f += g;
    h -= d; e ^= g << 14; g += h;
    #endif
  }
  void INIT(bool has_seed) {
    INT a, b, c, d, e, f, g, h;
    #if INT_MAX == 2147483647
    a = b = c = d = e = f = g = h = 0x9e3779b9UL;
    #else
    a = b = c = d = e = f = g = h = 0x9e3779b97f4a7c13ULL;
    #endif
    INT *m = mem;
    INT *r = rsl;
    if(!has_seed) z0 = z1 = z2 = 0;
    // scramble it
    uint32_t i;
    for (i = 0; i < 4; ++i)
      shuffle(a, b, c, d, e, f, g, h);
    if (has_seed) {
      // initialize using the contents of r[] as the seed
      for(i = 0; i < N; i += 8) {
        a += r[i]; b += r[i+1]; c += r[i+2]; d += r[i+3];
        e += r[i+4]; f += r[i+5]; g += r[i+6]; h += r[i+7];
        shuffle(a, b, c, d, e, f, g, h);
        m[i] = a; m[i+1] = b; m[i+2] = c; m[i+3] = d;
        m[i+4] = e; m[i+5] = f; m[i+6] = g; m[i+7] = h;
      }
      //do a second pass to make all of the seed affect all of m
      for (i = 0; i < N; i += 8) {
        a += m[i]; b += m[i+1]; c += m[i+2]; d += m[i+3];
        e += m[i+4]; f += m[i+5]; g += m[i+6]; h += m[i+7];
        shuffle(a, b, c, d, e, f, g, h);
        m[i] = a; m[i+1] = b; m[i+2] = c; m[i+3] = d;
        m[i+4] = e; m[i+5] = f; m[i+6] = g; m[i+7] = h;
      }
    }
    else {
      // fill in mm[] with messy stuff
      shuffle(a, b, c, d, e, f, g, h);
      m[i] = a; m[i+1] = b; m[i+2] = c; m[i+3] = d;
      m[i+4] = e; m[i+5] = f; m[i+6] = g; m[i+7] = h;
    }
    isaac();       // fill in the first set of results
    K = N;   // prepare to use the first set of results
  }
};

typedef cISAAC<8> ISAAC8;
typedef cISAAC<16> ISAAC16;

enum { UNIFORM = 0 };

struct NORMAL {
  enum { ZIGGURAT = 1, BOXMULLER = 2, INVERSE = 3 };
};

template <class G, int M>
class cPRNG : private G {
public:
  static const int BIT = G::BIT;
  typedef typename G::INT INT;
  cPRNG() { }
  void init(INT seed) { G::init(seed); }
  void init(INT init_key[], INT key_length) { G::init(init_key, key_length); }
  inline uint32_t gen_int32() { return G::gen_int32(); } //!< Return a random 32-bit integer within [0, 2^32-1] 
  inline uint64_t gen_int64() { return G::gen_int64(); } //!< Return a random 64-bit integer within [0, 2^64-1]
  inline double gen_close0_close1() { return this->template close0_close1<BIT>(); }
  inline double gen_close0_open1() { return this->template close0_open1<BIT>(); }
  inline double gen_open0_open1() { return this->template open0_open1<BIT>(); }
  inline double gen_close0_open1_res63() { return this->template close0_open1_res63<BIT>(); } //!< 64-bit precision value [0, 1)
private:
  template <int Bit> double close0_close1();
  template <int Bit> double close0_open1();
  template <int Bit> double open0_open1();
  template <int Bit> double close0_open1_res63();
};

/*!
 \def RND_IMPL
 \brief Specialization of template member functions for special template class
 */
#define RND_IMPL(G)\
template <> template <>\
inline double cPRNG<G, UNIFORM>::close0_close1<32>() {\
  return G::gen_int32() * (1.0 / 4294967295.0);\
}\
template <> template <>\
inline double cPRNG<G, UNIFORM>::close0_open1<32>() {\
  return G::gen_int32() * (1.0 / 4294967296.0);\
}\
template <> template <>\
inline double cPRNG<G, UNIFORM>::open0_open1<32>() {\
  return (G::gen_int32() + 0.5) * (1.0 / 4294967296.0);\
}\
template <> template <>\
inline double cPRNG<G, UNIFORM>::close0_open1_res63<32>() {\
  uint64_t res = (uint64_t)G::gen_int32() | (uint64_t)G::gen_int32() << 32;\
  return res * (1.0 / 18446744073709551616.0);\
}\
template <> template <>\
inline double cPRNG<G, UNIFORM>::close0_close1<64>() {\
  return G::gen_int64() * (1.0 / 18446744073709551615.0);\
}\
template <> template <>\
inline double cPRNG<G, UNIFORM>::close0_open1<64>() {\
  return G::gen_int64() * (1.0 / 18446744073709551616.0);\
}\
template <> template <>\
inline double cPRNG<G, UNIFORM>::open0_open1<64>() {\
  return (G::gen_int64() + 0.5) * (1.0 / 18446744073709551616.0);\
}\
template <> template <>\
inline double cPRNG<G, UNIFORM>::close0_open1_res63<64>() {\
  return G::gen_int64() * (1.0 / 18446744073709551616.0);\
}
RND_IMPL(RAN2)
RND_IMPL(LFSR113)
RND_IMPL(LFSR238)
RND_IMPL(GFSR4)
RND_IMPL(KISS)
RND_IMPL(SUPERKISS)
RND_IMPL(MWC256)
RND_IMPL(CMWC4096)
RND_IMPL(XOR4096)
RND_IMPL(MOTHER)
RND_IMPL(MT19937)
#ifdef __SSE2__
RND_IMPL(SFMT607)
RND_IMPL(SFMT1279)
RND_IMPL(SFMT2281)
RND_IMPL(SFMT4253)
RND_IMPL(SFMT12213)
RND_IMPL(SFMT19937)
RND_IMPL(SFMT44497)
RND_IMPL(SFMT86243)
RND_IMPL(SFMT132049)
RND_IMPL(SFMT216091)
#endif
RND_IMPL(WELL512)
RND_IMPL(WELL1024a)
RND_IMPL(WELL19937c)
RND_IMPL(WELL44497b)
RND_IMPL(ISAAC8)
RND_IMPL(ISAAC16)
#undef RND_IMPL

/*!
 \brief Gaussian/Normal Random Number Generators
 \brief Fast Ziggurat algorithm
 */
template <class G>
class cPRNG<G, NORMAL::ZIGGURAT> : private cPRNG<G, UNIFORM> {
public:
  static const int BIT = G::BIT;
  typedef typename G::INT INT;
  cPRNG() { INIT(); }
  void init(INT seed) { cPRNG<G, UNIFORM>::init(seed); }
  void init(INT init_key[], INT key_length) { cPRNG<G, UNIFORM>::init(init_key, key_length); }
  inline uint32_t gen_int32() { return cPRNG<G, UNIFORM>::gen_int32(); }
  inline uint64_t gen_int64() { return cPRNG<G, UNIFORM>::gen_int64(); }
  inline double gen_close0_close1() { return cPRNG<G, UNIFORM>::gen_close0_close1(); }
  inline double gen_close0_open1() { return cPRNG<G, UNIFORM>::gen_close0_open1(); }
  inline double gen_open0_open1() { return cPRNG<G, UNIFORM>::gen_open0_open1(); }
  inline double gen_close0_open1_res63() { return cPRNG<G, UNIFORM>::gen_close0_open1_res63(); }
  inline double gen_gauss() {
    int32_t hhhzzz = cPRNG<G, UNIFORM>::gen_int32(), iz = hhhzzz & 127;
    return abs(hhhzzz) < kn[iz] ? hhhzzz * wn[iz] : zignor(hhhzzz, iz);
  } //!< gauss random number with zero mean and unit variance
private:
  static uint32_t flag, kn[128];
  static double   wn[128], fn[128];
  inline double zignor(int32_t hhhzzz, int32_t iz) {
    double x, y;
    for( ; ; ) {
      x = hhhzzz * wn[iz];
      if (iz == 0) {
        do {
          x = log(cPRNG<G, UNIFORM>::gen_open0_open1()) * 0.29047645161474317645749646220079;
          y = -log(cPRNG<G, UNIFORM>::gen_open0_open1());
        } while (y + y < x * x);
        return hhhzzz > 0 ? 3.442619855899 - x : x - 3.442619855899;
      }
      if (fn[iz] + cPRNG<G, UNIFORM>::gen_open0_open1() * (fn[iz-1]-fn[iz]) < exp(-0.5 * x * x) ) return x;
      hhhzzz = cPRNG<G, UNIFORM>::gen_int32();
      iz = hhhzzz & 127;
      if (abs(hhhzzz) < kn[iz]) return hhhzzz * wn[iz];
    }
  }
  static void INIT() {
    if (flag == 1) {
      flag = 0;
      double m1 = 2147483648.0, dn = 3.442619855899, tn = 3.442619855899, vn = 9.91256303526217e-3;
      double q = vn / exp(-0.5 * dn * dn);
      kn[0] = uint32_t((dn / q) * m1);
      kn[1] = 0;
      wn[0] = q / m1;
      wn[127] = dn / m1;
      fn[0] = 1.0;
      fn[127] = exp(-0.5 * dn * dn);
      for(uint32_t i = 126; i > 0; --i) {
        dn = sqrt(-2.0 * log(vn / dn + exp(-0.5 * dn * dn)));
        kn[i+1] = uint32_t((dn / tn) * m1);
        tn = dn;
        fn[i] = exp(-0.5 * dn * dn);
        wn[i] = dn / m1;
      }
    }
  }
};
template <class G> uint32_t cPRNG<G, NORMAL::ZIGGURAT>::flag = 1;
template <class G> uint32_t cPRNG<G, NORMAL::ZIGGURAT>::kn[128];
template <class G> double cPRNG<G, NORMAL::ZIGGURAT>::wn[128];
template <class G> double cPRNG<G, NORMAL::ZIGGURAT>::fn[128];

template <class G>
class cPRNG<G, NORMAL::BOXMULLER> : private cPRNG<G, UNIFORM> {
public:
  static const int BIT = G::BIT;
  typedef typename G::INT INT;
  cPRNG() : flag(0) { }
  void init(INT seed) { cPRNG<G, UNIFORM>::init(seed); }
  void init(INT init_key[], INT key_length) { cPRNG<G, UNIFORM>::init(init_key, key_length); }
  inline uint32_t gen_int32() { return cPRNG<G, UNIFORM>::gen_int32(); }
  inline uint64_t gen_int64() { return cPRNG<G, UNIFORM>::gen_int64(); }
  inline double gen_close0_close1() { return cPRNG<G, UNIFORM>::gen_close0_close1(); }
  inline double gen_close0_open1() { return cPRNG<G, UNIFORM>::gen_close0_open1(); }
  inline double gen_open0_open1() { return cPRNG<G, UNIFORM>::gen_open0_open1(); }
  inline double gen_close0_open1_res63() { return cPRNG<G, UNIFORM>::gen_close0_open1_res63(); }
  inline double gen_gauss(double sigma = 1.0, double mean = 0.0) {
    if (flag == 1) {
      flag = 0;
      return x;
    }
    else {
      double r, v1, v2;
      do {
        v1 = cPRNG<G, UNIFORM>::gen_close0_close1() - 0.5;
        v2 = cPRNG<G, UNIFORM>::gen_close0_close1() - 0.5;
        r = v1 * v1 + v2 * v2;
      } while(r >= 0.25 || r == 0.0);
      r = sqrt(- 2.0 * (log(r) + 1.3862943611198906) / r) * sigma;
      flag = 1;
      x = v2 * r + mean;
      return v1 * r + mean;
    }
  } //!< gauss random number with zero mean and unit variance
private:
  int    flag;
  double x;
};

/*!
 \brief Inverse cumulative normal distribution function for generating gauss random number with zero mean and unit variance
 \brief c.f. http://home.online.no/~pjacklam/notes/invnorm/
 */
template <class G>
class cPRNG<G, NORMAL::INVERSE> : private cPRNG<G, UNIFORM> {
public:
  static const int BIT = G::BIT;
  typedef typename G::INT INT;
  cPRNG() { }
  void init(INT seed) { cPRNG<G, UNIFORM>::init(seed); }
  void init(INT init_key[], INT key_length) { cPRNG<G, UNIFORM>::init(init_key, key_length); }
  inline uint32_t gen_int32() { return cPRNG<G, UNIFORM>::gen_int32(); }
  inline uint64_t gen_int64() { return cPRNG<G, UNIFORM>::gen_int64(); }
  inline double gen_close0_close1() { return cPRNG<G, UNIFORM>::gen_close0_close1(); }
  inline double gen_close0_open1() { return cPRNG<G, UNIFORM>::gen_close0_open1(); }
  inline double gen_open0_open1() { return cPRNG<G, UNIFORM>::gen_open0_open1(); }
  inline double gen_close0_open1_res63() { return cPRNG<G, UNIFORM>::gen_close0_open1_res63(); }
  inline double gen_gauss() {
    double p = cPRNG<G, UNIFORM>::gen_open0_open1();
    if (p < 0.02425) {
      double x = -log(p), q = sqrt(x + x);
      x = (2.938163982698783 + q * (4.374664141464968 - q * (2.549732539343734 + q * (2.400758277161838 + q * (0.3223964580411365 + q * 0.007784894002430293))))) / (1.0 + q * (3.754408661907416 + q * (2.445134137142996 + q * (0.3224671290700398 + q * 0.007784695709041462))));
      return x;
    }
    else if (p > 0.97575) {
      double x = -log(1.0 - p), q = sqrt(x + x);
      x = -(2.938163982698783 + q * (4.374664141464968 - q * (2.549732539343734 + q * (2.400758277161838 + q * (0.3223964580411365 + q * 0.007784894002430293))))) / (1.0 + q * (3.754408661907416 + q * (2.445134137142996 + q * (0.3224671290700398 + q * 0.007784695709041462))));
      return x;
    }
    else {
      double q = p - 0.5, r = q * q;
      q *= (2.506628277459239 + r * (-30.66479806614716 + r * (138.3577518672690 + r * (-275.9285104469687 + r * (220.9460984245205 - r * 39.69683028665376))))) / (1.0 + r * (-13.28068155288572 + r * (66.80131188771972 + r * (-155.6989798598866 + r * (161.5858368580409 - r * 54.47609879822406)))));
      return q;
    }
  }
};

#endif
