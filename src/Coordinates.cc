// This file is part of LatticeTester.
//
// Copyright (C) 2012-2022  The LatticeTester authors, under the supervision
// of Pierre L'Ecuyer at Universit� de Montr�al.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "latticetester/Coordinates.h"

using namespace std;

namespace LatticeTester {

  ostream& operator<< (ostream& os, const Coordinates& coords)
  {
    os << "{";
    Coordinates::const_iterator it = coords.begin();
    if (it != coords.end()) {
      os << *it;
      while (++it != coords.end())
        os << "," << *it;
    }
    os << "}";
    return os;
  }

  istream& operator>> (istream& is, Coordinates& coords)
  {
    coords.clear();

    string digits = "0123456789";
    string sep = " \t,";

    // check if coordinate set is enclosed in braces
    bool with_braces = false;
    if (is.peek() == '{') {
      is.get();
      with_braces = true;
    }

    while (true) {
      if (with_braces && is.peek() == '}') {
        is.get();
        break;
      }
      if (digits.find(is.peek()) != string::npos) {
        // digit found
        Coordinates::value_type val;
        is >> val;
        coords.insert(val);
        continue;
      }
      if (sep.find(is.peek()) != string::npos) {
        // discard separator character
        is.get();
        continue;
      }
      break;
    }
    return is;
  }

} // end namespace
