// SPDX-License-Identifier: BSD-3-Clause
#pragma once
namespace ProcessLib::StructuralElements
{
/// Individual member convention: EA [N], EI [N m^2], corrected GA [N].
struct Beam2DProcessData
{
    double EA;
    double EI;
    double GA;
    double axial_capacity = 0.;
    double moment_capacity = 0.;
};
}
