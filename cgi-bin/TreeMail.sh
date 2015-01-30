# Small Genome Tools
# Copyright (C) 2015 University of the Witwatersrand, Johannesburg, South Africa
# Author: Dr Trevor G. Bell, TrevorGrahamBell@gmail.com

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

#!/bin/bash

# align and/or create .phy file

mv "$1" infile

if [ -f outfile ]; then
rm outfile
fi

if [ -f outtree ]; then
rm outtree
fi

if [ -f dnadist.cmd ]; then
rm dnadist.cmd
fi

if [ -f neighbor.cmd ]; then
rm neighbor.cmd
fi

echo D > dnadist.cmd
echo L >> dnadist.cmd
echo 2 >> dnadist.cmd
echo Y >> dnadist.cmd

phylip dnadist < dnadist.cmd

rm infile

mv outfile infile

echo L > neighbor.cmd
echo 2 >> neighbor.cmd
echo Y >> neighbor.cmd

phylip neighbor < neighbor.cmd

rm infile
mv outfile result.txt
mv outtree result.tre

if [ -f outfile ]; then
rm outfile
fi

if [ -f outtree ]; then
rm outtree
fi

if [ -f dnadist.cmd ]; then
rm dnadist.cmd
fi

if [ -f neighbor.cmd ]; then
rm neighbor.cmd
fi

