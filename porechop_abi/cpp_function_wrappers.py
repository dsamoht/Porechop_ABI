"""
Copyright 2017 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Porechop

Porechop originally used C++ functions compiled in cpp_functions.so, compatibility.so,
approx_counter and msa_consensus. They are replaced by the Rust extension module _native (built
from porechop_rs), which gives the same results. This module re-exports its functions.

This file is part of Porechop. Porechop is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Porechop is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Porechop. If
not, see <http://www.gnu.org/licenses/>.
"""

from ._native import align_adapter as native_align_adapter  # noqa: F401
from ._native import align_regions, find_middle_hits, compatibility_matrix  # noqa: F401
from ._native import msa_consensus, approx_counter  # noqa: F401
