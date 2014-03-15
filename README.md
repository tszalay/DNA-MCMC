DNA-MCMC
========

3D Markov Chain Monte Carlo code written in Python, written for Harvard's AM207 course in a solid 36 hours.

Dependencies:
 - PyOpenGL
 - freeglut

For a description of what the code is doing, please refer to am207_final.pdf.

The functions that perform the calculations referenced in the writeup are in main.py.
Each simulation condition (free, stretched, and in a nanopore) are in the corresponding sim_free, sim_stretched, and sim_pore. The file sim_draw.py runs the 3D viewer for a given array of Xs, and sim_replay.py runs it for the saved runs in /data (press number keys to toggle dataset, arrows+a/z to look, space to advance sample). These require freeglut and PyOpenGL as a (slightly annoying) dependency, sorry. There are also some vestigal bits of code left over from a cython experiment that didn't amount to a significant enough speedup.

Feel free to use or extend this code for any purpose whatsoever!

