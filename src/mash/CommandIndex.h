// Copyright © 2015, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen,
// Sergey Koren, and Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

#ifndef INCLUDED_CommandIndex
#define INCLUDED_CommandIndex

#include "Command.h"

namespace mash {

class CommandIndex : public Command
{
public:

    CommandIndex();

    int run() const; // override
};

} // namespace mash

#endif
