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

void Tabulate(std::uint8_t);
void Canonical(std::uint8_t);
void FixedPoint(std::uint8_t, const std::vector<std::string>&);
void Cycle(std::uint8_t, std::uint32_t);
void Minimal(std::uint8_t, std::uint32_t);
void Run(std::uint8_t, std::uint8_t, const std::vector<std::string>&);

bool Injective(std::uint8_t r);
bool Surjective(std::uint8_t r);
bool InDegree(std::uint8_t r, std::uint32_t k);

} // namespace omega

#endif // OMEGA_ECA_H
