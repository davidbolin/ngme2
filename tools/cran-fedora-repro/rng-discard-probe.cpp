#include <random>

int main() {
  std::mt19937 rng(123), reference(123);
  (void)reference();
#ifdef NGME_EXPLICIT_DISCARD
  (void)rng();
#else
  rng();
#endif
  return rng == reference && rng() == reference() ? 0 : 1;
}
