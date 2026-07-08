// Copyright (c) 2025  Laura Moreno, Eduardo Sáenz de Cabezón

#include "CoCoA/library.H"
#include <ctime>
#include <random>

#include <algorithm>
#include <libteddy/core.hpp>
#include <libteddy/reliability.hpp>
#include <cstdint>

using namespace std;
using namespace teddy::ops;
using bdd_m = teddy::bdd_manager;

//----------------------------------------------------------------------
//const string ShortDescription =
//  "Collection of functions for connecting monomial ideals and BDDs.  \n"
//  "Suitable for examples and programs using monomial ideals and BDDs.    \n";

//const string LongDescription =
//  "Collection of functions for connecting monomial ideals and BDDs.  \n"
//  "Suitable for examples and programs using monomial ideals and BDDs.    \n";
//----------------------------------------------------------------------

namespace CoCoA
{
    
    /*PRE: The method expects a vector of elements `generators`, an instance `m` of the `bdd_m` class and a reference `myS` to a `diagram_t` object.
     `generators` should contain the generators from which the Binary Decision Diagram (BDD) will be constructed.
     `m` provides the necessary operations to construct the BDD.
     `myS` will represent the resulting BDD.

     POST: - `myS` contains the resulting BDD from the combination of variables in the generators.
        - The resulting BDD `myS` represents the BDD of the ideal created by the generators`.*/
      
    void generateBDD(const vector<RingElem>& generators, bdd_m& m, bdd_m::diagram_t& myS) {
        myS = m.constant(0);
        for (RingElem j : generators) {
            bdd_m::diagram_t local = m.constant(1);
            for (int i : IndetsIn(LPP(j))) {
                local = m.apply<AND>(local, m(i));
            }
            myS = m.apply<OR>(myS, local);
        }
    }
    
    /* PRE: El anillo de polinomios con variables x0, ..., x_{n-1} y coeficientes racionales (RingQQ) debe poder construirse correctamente.
        - n > 0: el número de variables debe ser positivo.
        - D > 0: el grado máximo de los monomios debe ser positivo.
        - c > 0: el parámetro de densidad debe ser positivo.
        
        POST: Devuelve un ideal I 0-dimensional del anillo Q[x0, ..., x_{n-1}] generado por un subconjunto aleatorio de monomios
        de grado total ≤ D, elegidos según el modelo ER con probabilidad p = c/D.
        Si no se obtiene un ideal cero-dimensional en un intento, el proceso se repite hasta encontrar uno.
        Si no existen monomios seleccionados en algún intento, se descarta y se repite (nunca devuelve I = (0)).
    */
    ideal createIdealER(int n, int D, double c){
        double p = c/double(D);
        //---- Ring and monomials ----
        const PPMonoid PPM1 = NewPPMonoidEv(SymbolRange("x", 0, n-1), lex);
        const DivMaskRule DMR = NewDivMaskEvenPowers();
        SparsePolyRing P = NewPolyRing(RingQQ(), PPM1);

        vector<long> exp(n,0);   // exponents
        vector<RingElem> G;      // generators
        vector<PPMonoidElem> candidates;

        // Generator of random numbers
        mt19937 rng((unsigned)time(nullptr));
        bernoulli_distribution bern(p);

        // Function for generating all the possible monomials
        function<void(int,int)> genExps = [&](int idx, int sum) {
            if(idx == n){
                if(sum >= 1 && sum <= D){
                    candidates.push_back(PPMonoidElem(PPM1, exp));
                }
                return;
            }
            for(int k = 0; k <= D - sum; ++k){
                exp[idx] = k;
                genExps(idx + 1, sum + k);
            }
            exp[idx] = 0;
        };
        genExps(0,0);
        
        //int i = 0;
        while(true)
        {
            PPVector gens(PPM1, DMR);
            for(auto& m : candidates){
                if(bern(rng)) gens.myPushBack(m);
            }
            if(len(gens) == 0){
                cout << "I=(0)" <<endl;
                continue;
            }
            InterreduceSort(gens);
            convert(G, P, gens);
            ideal I(G);
            if(IsZeroDim(I))
            {
                return I;
            }
        }
    }
    
    /*PRE: The method expects a RingElem `mon`, a reference `myS` to a `diagram_t` object representing a Binary Decision Diagram (BDD), a reference `mu` to a PPMonoidElem object, an integer `n` representing the number of variables in the BDD and an instance `m` of the bdd_m class.

    POST: - `mu` contains the PPMonoidElem representation of the input monomial `mon`.
    - The boolean vector `boolmu` is generated, where `boolmu[i]` is true if the variable `i` is present in `mu`.
    - The method returns true if the evaluation of the BDD `myS` with the boolean vector `boolmu` returns true, indicating that `mon` is an element of the BDD.*/
    
    bool IsElem_BDD(RingElem mon, bdd_m::diagram_t& myS, PPMonoidElem& mu, int n, bdd_m& m) {
        mu=PP(BeginIter(mon));
        vector boolmu=vector(n,false);
        for(int i: IndetsIn(mu))
          boolmu[i]=true;
        return m.evaluate(myS,boolmu);
    }
    
    /*PRE: The method expects a RingElem `mon` and an ideal `I`.
    - `mon` represents the monomial to be evaluated.
    - `I` is the ideal against which `mon` is to be tested for membership.
    - The method internally generates a Binary Decision Diagram (BDD) representing the ideal `I` and evaluates `mon` in the BDD.

    POST: - The method returns true if the input monomial `mon` is found to be an element of the ideal `I`, as determined by the BDD evaluation.*/

    bool IsElem_BDD_Ideal(RingElem mon, ideal I) {
        vector<RingElem> g = gens(I);
        const int n = NumIndets(RingOf(I));
        bdd_m m(n,100000);
        bdd_m::diagram_t myS;
        generateBDD(g,m,myS);
        PPMonoidElem mu(PPM(RingOf(I)));
        return IsElem_BDD(mon,myS,mu,n,m);
    }
    
    /* PRE: The method expects a valid ideal `I` in the context of the code.

    POST: The method returns the `qBasis` vector, which represents a minimal set of elements that generates the quotient space. */

    vector<PPMonoidElem> quotientBasis_BDD(ideal I) {
        vector<RingElem> g = gens(I);
        const int n = NumIndets(RingOf(I));
        bdd_m m(n, 100000);
        bdd_m::diagram_t myS;
        generateBDD(g, m, myS);
        
        vector<vector<int>> sa = m.satisfy_all<vector<int>>(0, myS);

        vector<PPMonoidElem> qBasis;
        const PPMonoid p = PPM(RingOf(I));
        for (const vector<int>& assignment : sa) {
            vector<long> exponents;
            for (int i = 0; i < n; ++i) {
                exponents.push_back(assignment[i]);
            }
            PPMonoidElem monomial(p, exponents);
            qBasis.push_back(monomial);
        }
        
        return qBasis;
    }
    
    /*PRE: The method expects a reference `myS` to a `diagram_t` object representing a Binary Decision Diagram (BDD), an integer `n` representing the number of variables in the BDD and an instance `m` of the bdd_m class.
    
    POST: The method returns the `I` ideal equivalent to the BDD passed as a parameter.*/
    ideal reconstruction_ideal(bdd_m& m, bdd_m::diagram_t& myS, int n) {
        vector<vector<int>> s = m.satisfy_all<vector<int>>(1, myS);
        
        const PPMonoid PPM = NewPPMonoidEv(SymbolRange("x", 0, n-1), lex);
        const DivMaskRule DMR = NewDivMaskEvenPowers();
        PPVector PPV(PPM,DMR);
        SparsePolyRing P = NewPolyRing(RingQQ(), PPM);
        vector<vector<long>> exponents;
        for (const vector<int>& assignment : s) {
            vector<long> exponent;
            for (int i = 0; i < n; ++i) {
                exponent.push_back(assignment[i]);
            }
            exponents.push_back(exponent);
        }
        for (const vector<long>& exp : exponents) {
            PPMonoidElem term(PPM, exp);
            PPV.myPushBack(term);
        }
        
        InterreduceSort(PPV);
        vector<RingElem> gen;
        convert(gen, P, PPV);
        
        ideal I = ideal(gen);
        
        return I;
        
    }
    
    /* PRE: The ideals I1 and I2 must belong to the same ring.

    POST: The method returns the ideal `I`, representing the intersection of `I1` and `I2`, after being processed and reduced.*/
    ideal intersection_BDD(ideal I1, ideal I2) {
        const int n = NumIndets(RingOf(I1));
        bdd_m m(n, 100000);
        bdd_m::diagram_t myS1;
        bdd_m::diagram_t myS2;
    
        vector<RingElem> g1 = gens(I1);
        vector<RingElem> g2 = gens(I2);
        
        generateBDD(g1, m, myS1);
        generateBDD(g2, m, myS2);
        
        bdd_m::diagram_t myS_int = m.apply<AND>(myS1, myS2);
        
        ideal I = reconstruction_ideal(m, myS_int, n);
        
        return I;
    }
    
    /* PRE: The ideals I1 and I2 must belong to the same ring.

    POST: The method returns the ideal `I`, representing the union of `I1` and `I2`, after being processed and reduced.*/
    ideal sum_BDD(ideal I1, ideal I2) {
        const int n = NumIndets(RingOf(I1));
        bdd_m m(n, 100000);
        bdd_m::diagram_t myS1;
        bdd_m::diagram_t myS2;
    
        vector<RingElem> g1 = gens(I1);
        vector<RingElem> g2 = gens(I2);
        
        generateBDD(g1, m, myS1);
        generateBDD(g2, m, myS2);
        
        bdd_m::diagram_t myS_union = m.apply<OR>(myS1, myS2);
        
        ideal I = reconstruction_ideal(m, myS_union, n);
        
        return I;
    }
    
    /* PRE: The method expects a valid ideal `I` in the context of the code.

    POST: The method returns the height of the ideal `I`.*/
    int height(ideal I) {
        vector<RingElem> g = gens(I);
        const int n = NumIndets(RingOf(I));
        bdd_m m(n, 100000);
        bdd_m::diagram_t myS;
        generateBDD(g, m, myS);
        vector<vector<int>> sat = m.satisfy_all<vector<int>>(1, myS);
        // Filter minimum supports (by inclusion)
        vector<vector<int>> minimal;
        for (vector<int> &a : sat) {
            bool has_smaller = false;
            for (vector<int> &b : sat) {
                if (a == b) continue;
                bool divides = true;
                for (int i = 0; i < n; ++i) {
                    if (b[i] > a[i]) {
                        divides = false;
                        break;
                    }
                }
                if (divides) {
                    has_smaller = true;
                    break;
                }
            }
            if (!has_smaller) minimal.push_back(a);
        }
        // Calculate the min hitting set
        int tau = n + 1;
        vector<int> subset(n);
        for (int k = 0; k <= n; ++k) {
            bool found = false;
            vector<int> idx(n);
            fill(idx.begin(), idx.begin() + k, 1);
            do {
                vector<int> current(n);
                for (int i = 0; i < n; ++i)
                    current[i] = idx[i];
                bool ok = true;
                for (vector<int> &edge : minimal) {
                    bool hits = false;
                    for (int j = 0; j < n; ++j)
                        if (edge[j] && current[j]) {
                            hits = true;
                            break;
                        }
                    if (!hits) { ok = false; break; }
                }
                if (ok) {
                    tau = k;
                    found = true;
                    break;
                }
            } while (prev_permutation(idx.begin(), idx.end()));
            if(found) break;
        }
        return tau;
    }

    /*PRE: The method expects a reference `myS` to a `diagram_t` object representing a Binary Decision Diagram (BDD), an integer `n` representing the number of variables in the BDD and an instance `m` of the bdd_m class.
    
    POST: The method returns a list with the elements of the minimal generating set of the ideal `I` represented by the BDD.*/
    vector<vector<long>> minGenSet(int n, bdd_m& manager, bdd_m::diagram_t& myS){
        vector<vector<long>> sat = manager.satisfy_all<vector<long>>(1, myS);

        vector<vector<long>> minimal;

        for (vector<long> &a : sat) {
            bool has_smaller = false;
            for (vector<long> &b : sat) {
                if (a == b) continue;

                bool divides = true;
                for (int i = 0; i < n; ++i) {
                    if (b[i] > a[i]) { divides = false; break; }
                }

                if (divides) { has_smaller = true; break; }
            }

            if (!has_smaller) minimal.push_back(a);
        }
        return minimal;
    }
}

