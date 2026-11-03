// A re-implementation of NumPy's default random stream (SeedSequence + PCG64, XSL-RR output), so that the
// generator reproduces the instance files of the paper from the seed alone, in C++. Checked against numpy 2.0.2:
// SeedSequence(1).generate_state(4, uint64), the doubles and integers of default_rng(1|3|7|11), and every committed
// instance file (experiments/scripts/reproduce.sh instances).
#pragma once
#include <cstdint>
#include <vector>

namespace nprng {

// numpy.random.SeedSequence with a single integer entropy word and the default pool of 4 words.
struct SeedSequence {
  static constexpr uint32_t INIT_A = 0x43b0d7e5u, MULT_A = 0x931e8875u, INIT_B = 0x8b51f9ddu, MULT_B = 0x58f38dedu,
                            MIX_MULT_L = 0xca01f9ddu, MIX_MULT_R = 0x4973f715u, XSHIFT = 16u, POOL = 4u;
  uint32_t pool[POOL];
  static uint32_t hashmix(uint32_t v, uint32_t& h) { v ^= h; h *= MULT_A; v *= h; v ^= v >> XSHIFT; return v; }
  static uint32_t mix(uint32_t x, uint32_t y) { uint32_t r = MIX_MULT_L * x - MIX_MULT_R * y; r ^= r >> XSHIFT; return r; }
  explicit SeedSequence(const std::vector<uint32_t>& entropy) {
    uint32_t h = INIT_A;
    for (uint32_t i = 0; i < POOL; ++i) pool[i] = hashmix(i < entropy.size() ? entropy[i] : 0u, h);
    for (uint32_t s = 0; s < POOL; ++s)
      for (uint32_t d = 0; d < POOL; ++d)
        if (s != d) pool[d] = mix(pool[d], hashmix(pool[s], h));
    for (size_t s = POOL; s < entropy.size(); ++s)
      for (uint32_t d = 0; d < POOL; ++d) pool[d] = mix(pool[d], hashmix(entropy[s], h));
  }
  explicit SeedSequence(uint64_t seed)
      : SeedSequence(seed >> 32 ? std::vector<uint32_t>{(uint32_t)seed, (uint32_t)(seed >> 32)} : std::vector<uint32_t>{(uint32_t)seed}) {}
  // generate_state(n, np.uint64): 2n words of 32 bits, viewed little-endian as n words of 64.
  std::vector<uint64_t> generate_state64(int n) const {
    std::vector<uint32_t> w(2 * n);
    uint32_t h = INIT_B;
    for (int i = 0; i < 2 * n; ++i) {
      uint32_t v = pool[i % POOL]; v ^= h; h *= MULT_B; v *= h; v ^= v >> XSHIFT; w[i] = v;
    }
    std::vector<uint64_t> out(n);
    for (int i = 0; i < n; ++i) out[i] = (uint64_t)w[2 * i] | ((uint64_t)w[2 * i + 1] << 32);
    return out;
  }
};

// numpy.random.PCG64 (the default_rng bit generator): 128-bit LCG state, XSL-RR 64-bit output.
struct PCG64 {
  using u128 = unsigned __int128;
  static constexpr uint64_t MULT_HI = 2549297995355413924ULL, MULT_LO = 4865540595714422341ULL;
  u128 state = 0, inc = 0;
  uint64_t buf32 = 0; bool has32 = false;    // numpy's pcg64_next32: each 64-bit draw feeds two 32-bit draws, low word first
  static u128 c128(uint64_t hi, uint64_t lo) { return ((u128)hi << 64) + lo; }
  void step() { state = state * c128(MULT_HI, MULT_LO) + inc; }
  explicit PCG64(uint64_t seed) {
    std::vector<uint64_t> w = SeedSequence(seed).generate_state64(4);
    u128 initstate = c128(w[0], w[1]), initseq = c128(w[2], w[3]);
    state = 0; inc = (initseq << 1) | 1; step(); state += initstate; step();
  }
  uint64_t next64() {
    step();
    uint64_t xsl = (uint64_t)(state >> 64) ^ (uint64_t)state;
    unsigned rot = (unsigned)(state >> 122);
    return (xsl >> rot) | (xsl << ((-rot) & 63));
  }
  uint32_t next32() {
    if (has32) { has32 = false; return (uint32_t)buf32; }
    uint64_t v = next64(); buf32 = v >> 32; has32 = true; return (uint32_t)v;
  }
  // Generator.random(): 53 random bits scaled to [0, 1).
  double random() { return (double)(next64() >> 11) * (1.0 / 9007199254740992.0); }
  // Generator.uniform(lo, hi)
  double uniform(double lo, double hi) { return lo + (hi - lo) * random(); }
  // Generator.integers(n), 0 <= x < n, int64 dtype: Lemire's bounded method on the buffered 32-bit stream when the
  // range fits in 32 bits (numpy random_bounded_uint64 with use_masked = false), else on 64-bit draws.
  uint64_t integers(uint64_t n) {
    if (n <= 1) return 0;
    uint64_t rng = n - 1, rng_excl = n;
    if (rng <= 0xFFFFFFFFULL) {
      if (rng == 0xFFFFFFFFULL) return next32();
      uint64_t m = (uint64_t)next32() * rng_excl; uint32_t leftover = (uint32_t)m;
      if (leftover < rng_excl) {
        uint32_t threshold = (uint32_t)((0xFFFFFFFFULL - rng) % rng_excl);
        while (leftover < threshold) { m = (uint64_t)next32() * rng_excl; leftover = (uint32_t)m; }
      }
      return m >> 32;
    }
    if (rng == 0xFFFFFFFFFFFFFFFFULL) return next64();
    u128 m = (u128)next64() * rng_excl; uint64_t leftover = (uint64_t)m;
    if (leftover < rng_excl) {
      uint64_t threshold = (0xFFFFFFFFFFFFFFFFULL - rng) % rng_excl;
      while (leftover < threshold) { m = (u128)next64() * rng_excl; leftover = (uint64_t)m; }
    }
    return (uint64_t)(m >> 64);
  }
};

}  // namespace nprng
