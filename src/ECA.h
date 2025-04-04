/**
  * <ECA.h>
  *   Functions for expressing properties about the elementary cellular
  *   automata.
 **/

#ifndef OMEGA_ECA_H
#define OMEGA_ECA_H

#include <cstdint>
#include <string>
#include <vector>

namespace omega {

void Tabulate(uint32_t);
void Canonical(uint8_t);
void Minimal(uint8_t, uint32_t);
void Run(uint8_t, uint8_t, const std::vector<std::string>&);
bool RightShift(uint32_t, uint32_t, const std::vector<std::string>&);

bool Injective(uint32_t);
bool Surjective(uint32_t);
bool FixedPoint(uint32_t);
bool Cycle(uint32_t, uint32_t);
bool Nilpotent(uint32_t, uint32_t);
bool InDegree(uint32_t, uint32_t);

} // namespace omega

#endif // OMEGA_ECA_H
