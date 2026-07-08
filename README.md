# MonomialIdeal-BDD: A BDD-engine for computations on monomial ideals

This repository contains a high-performance C++ engine designed for algebraic and combinatorial computations on **square-free monomial ideals**. By leveraging the bijective correspondence between these ideals and **monotone Boolean functions**, the engine utilizes **Binary Decision Diagrams (BDDs)** to overcome classical symbolic bottlenecks.

## Core Concept
The project is built on the fact that every square-free monomial ideal $I$ corresponds to a unique monotone Boolean function $\mathcal{B}_I$. 
*   **Symbolic approach:** Traditional systems (CoCoALib, Macaulay2) manipulate minimal generating sets $G(I)$, which often suffer from "intermediate expression swell" during complex operations like intersections.
*   **BDD approach:** Our engine uses **Reduced Ordered BDDs (ROBDDs)**. This structure provides a **canonical representation** and enables **node sharing**, where identical Boolean sub-functions share physical memory, significantly reducing redundancy.

## Key Features & Operations
The engine translates algebraic operations into native logical manipulations:
*   **Sum of Ideals:** Map to the Boolean **OR** ($\vee$) operator.
*   **Intersection:** Map to the Boolean **AND** ($\wedge$) operator.
*   **Colon Ideal:** Implemented via an iterative approach over generator sets.
*   **Standard Monomials:** Efficiently calculated by counting paths to the terminal "False" node ($\bot$) in the BDD.
*   **Alexander Dual:** Obtained by computing the irreducible DNF, equivalent to the monotone dualization problem.

## Technology Stack
The implementation integrates two specialized C++ libraries:
1.  **CoCoALib:** Used for managing algebraic structures, polynomial rings, and monomial generators.
2.  **TeDDy:** A templated header-only library used as the core BDD engine to represent and manipulate the ideals.
3.  **Macaulay2:** Employed for generating random test cases and as a benchmark for performance comparison.

## Installation & Configuration
The project requires the **C++20 standard** (mandated by the TeDDy library). To integrate it with CoCoALib, you must modify the `autoconf.mk` file:
```make
# Compilation settings
CXXFLAGS= -std=c++20 -Wall -pedantic -fPIC -O2

# Include TeDDy library path
TEDDY=~/path/to/libteddy/
LDLIBS= ... -I$(TEDDY)
```

## Performance Highlights
Experimental results demonstrate significant efficiency gains over traditional symbolic engines:
*   **Intersections:** Consistently **2 to 3 orders of magnitude faster** than CoCoALib and Macaulay2.
*   **Standard Monomials:** Up to **4 orders of magnitude faster**. The BDD counts the basis length without needing to enumerate every monomial.
*   **Alexander Dual:** Shows remarkable stability and predictable performance by preventing the growth of intermediate generator sets.
*   **Staircase Ideals:** While traditional systems suffer from a performance gap between sum and intersection, the BDD engine maintains a stable and consistent profile regardless of the operation type.

## Implementation Examples
### Generating a BDD from an Ideal
```cpp
void generateBDD(const std::vector<RingElem>& generators, bdd_m& m, bdd_m::diagram_t& myS) {
    myS = m.constant(0);
    for (auto j : generators) {
        auto local = m.constant(1);
        for (auto i : IndetsIn(LPP(j))) {
            local = m.apply<AND>(local, m(i));
        }
        myS = m.apply<OR>(myS, local);
    }
}
```

## Credits & Funding
*   **Authors:** Laura Moreno-Resa and Eduardo Sáenz-de-Cabezón (Universidad de La Rioja).
*   **Funding:** This work is partially supported by grant **PID2024-157733NBI00** funded by MCIN/AEI/10.13039/501100011033 and FEDER EU.

## References
For more details, please refer to the original paper: 
*Moreno-Resa, L., & Sáenz-de-Cabezón, E. (2026). A BDD-engine for computations on monomial ideals. In 51st International Symposium on Symbolic and Algebraic Computation (ISSAC ’26).*
