# N-body simulation

This software is a straight-forward implementation of numerical solver for the N-body problem, that is the movement of N physical objects (like planets or stars) under the influence of their mutual gravity. The equations of movement are a system of N second order ordinary differential equations where the source term is given by Newton's law of universal gravitation.

Please note that this implementation is - as said - straight forward, meaning that its runtime complexity is O(N^2) (faster implementations with complexity O(N*log(N) are available) and it will only work for bodies whose product of gravity constant G and mass m_j is in the order of at most 10 or 1000 (for real planets like earth it is in the order of 10^13). This is due to the fact that the integrated system is not un-dimensionalized (the same holds true for the distances between the bodies).
