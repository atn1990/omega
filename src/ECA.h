/**
  * <ECA.h>
  *   Functions for expressing properties about the elementary cellular
  *   automata.
 **/

#ifndef OMEGA_ECA_H
#define OMEGA_ECA_H

#include "Util.h"

#include <string>
#include <vector>

namespace omega {

enum class ShiftType {
  Left,
  Right,
};

void Tabulate(int_type);
void Canonical(int_type);
void Minimal(int_type, int_type);
void Run(int_type, int_type, const std::vector<std::string>&);

bool Injective(int_type);
bool Surjective(int_type);
bool FixedPoint(int_type);
bool Cycle(int_type, int_type);
bool Nilpotent(int_type, int_type);
bool Shift(int_type, int_type, ShiftType);
bool InDegree(int_type, int_type);

} // namespace omega

#endif // OMEGA_ECA_H
