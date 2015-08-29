Monte Carlo Simulation of Krypton
=================================

A toy simulation of krypton, with:
- minimum image convention in a cube
- periodic boundary condition
- Lennard Jones Potential
- NPT or NVT enssemble

some parameters are copied from http://chryswoods.com/book/export/html/100

# Usage

You need ruby2.0 or later to run it.
`./metropolis.rb` in a terminal will run 300 equilibrium moves and 300 moves for analysis, in NPT enssemble.
`./metropolis.rb 20 100 nvt -v` will run 20 equilibrium moves, 100 moves for analysis, in NVT enssemble,
and write a pdb file.

More moves for analysis will write more pdb files.
If you have vmd installed, you may `vmd *.pdb` to view a movie (though the movie has no physical meaning).
