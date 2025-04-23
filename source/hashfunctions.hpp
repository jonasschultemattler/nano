// #include <farmhash.h>
// #include <murmurhash.h>

static inline constexpr uint64_t hash_2(const uint64_t x) {
    std::hash<uint64_t> hasher;
    return hasher(x);
}

static inline constexpr uint64_t hash_3(const uint64_t x) {
    __uint128_t result = x;
    result *= 0x9E3779B97F4A7C15ULL;
    return static_cast<uint64_t>(result) ^ static_cast<uint64_t>(result >> 64);
}