// Copyright 2008-2009 Deutsches Forschungszentrum fuer Kuenstliche Intelligenz
// or its licensors, as applicable.
//
// You may not use this file except under the terms of the accompanying license.
//
// Licensed under the Apache License, Version 2.0 (the "License"); you
// may not use this file except in compliance with the License. You may
// obtain a copy of the License at http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
// Project: ocrofst
// File: beam-search.cc
// Purpose: beam search
// Responsible: mezhirov
// Reviewer:
// Primary Repository:
// Web Sites: www.iupr.org, www.dfki.de, www.ocropus.org

#include "ocr-pfst.h"
#include "fst-heap.h"

using namespace colib;
using namespace ocrofst;

namespace {

    /// A SearchTree contains all vertices that were ever touched during the
    /// search, and can produce a prehistory for every ID.
    struct SearchTree {
        intarray parents;
        intarray inputs;
        intarray outputs;
        intarray v1; // vertices from FST 1
        intarray v2; // vertices from FST 2
        floatarray costs;

        void clear() {
            parents.clear();
            inputs.clear();
            outputs.clear();
            v1.clear();
            v2.clear();
            costs.clear();
        }

        void get(intarray &r_vertices1,
                 intarray &r_vertices2,
                 intarray &r_inputs,
                 intarray &r_outputs,
                 floatarray &r_costs,
                 int id) {
            intarray t_v1; // vertices
            intarray t_v2; // vertices
            intarray t_i; // inputs
            intarray t_o; // outputs
            floatarray t_c; // costs
            int current = id;
            while(current != -1) {
                t_v1.push(v1[current]);
                t_v2.push(v2[current]);
                t_i.push(inputs[current]);
                t_o.push(outputs[current]);
                t_c.push(costs[current]);
                current = parents[current];
            }

            reverse(r_vertices1, t_v1);
            reverse(r_vertices2, t_v2);
            reverse(r_inputs, t_i);
            reverse(r_outputs, t_o);
            reverse(r_costs, t_c);
        }

        int add(int parent, int vertex1, int vertex2,
                   int input, int output, float cost) {
            int n = parents.length();
            //logger.format("stree: [%d]: parent %d, v1 %d, v2 %d, cost %f",
            //               n, parent, vertex1, vertex2, cost);
            parents.push(parent);
            v1.push(vertex1);
            v2.push(vertex2);
            inputs.push(input);
            outputs.push(output);
            costs.push(cost);
            return n;
        }
    };


    ///
    /// A bunch of intervals corresponding to positions of special symbols
    /// in some array. A symbol, say L_PHI, occupies the slice from 
    /// start(L_PHI) to end(L_PHI).
    ///
    /// The input array must be sorted.
    ///
    class SpecialSymbolRanges
    {
        enum { SIZE = 5 };

        int starts[SIZE];
        int ends[SIZE];

        void init() {
            memset(starts, 0, sizeof(starts));
            memset(ends, 255, sizeof(ends));
        }

        void checkSymbols() {
            CHECK(-L_SIGMA < SIZE);
            CHECK(-L_RHO < SIZE);
            CHECK(-L_PHI < SIZE);
            CHECK(-L_EPSILON < SIZE);

            CHECK(L_SIGMA <= 0);
            CHECK(L_RHO <= 0);
            CHECK(L_PHI <= 0);
            CHECK(L_EPSILON <= 0);
        }


    public:
        /// The index of the first non-special (positive) value.
        int specialsEnd;
        
        SpecialSymbolRanges(intarray &a) {
            checkSymbols();
            init();
            int i = 0;
            while(i < a.length() && a[i] <= 0)
            {
                int start = i;
                int symbol = a[start];
                do { i++; } while(i < a.length() && a[i] == symbol);
                CHECK(-symbol < SIZE);
                starts[-symbol] = start;
                ends[-symbol] = i;
            }
            specialsEnd = i;
        }

        int start(int symbol)
        {
            return starts[-symbol];
        }

        int end(int symbol)
        {
            return ends[-symbol];
        }

        int count(int symbol)
        {
            if (end(symbol) < 0)
                return 0;
            return end(symbol) - start(symbol);
        }
    };


    struct BeamSearch {
        OcroFST &fst1;
        OcroFST &fst2;
        SearchTree stree;

        intarray beam; // indices into stree
        floatarray beamcost; // global cost, corresponds to the beam

        PriorityQueue nbest;
        intarray all_inputs;
        intarray all_targets1;
        intarray all_targets2;
        intarray all_outputs;
        floatarray all_costs;
        intarray parent_trails; // indices into the beam
        int beam_width;
        int accepted_from1;
        int accepted_from2;
        float g_accept;   // best cost for accept so far
        int best_so_far;  // ID into stree (-1 for start)
        float best_cost_so_far;

        BeamSearch(OcroFST &fst1, OcroFST &fst2, int beam_width):
                fst1(fst1),
                fst2(fst2),
                nbest(beam_width),
                beam_width(beam_width),
                accepted_from1(-1),
                accepted_from2(-1) {
        }

        void clear() {
            nbest.clear();
            all_targets1.clear();
            all_targets2.clear();
            all_inputs.clear();
            all_outputs.clear();
            all_costs.clear();
            parent_trails.clear();
        }

        // This looks at the transition from state pair
        // (f1,f2) -> (t1,t2), with the given cost.

        void relax(int f1, int f2,   // input state pair
                   int t1, int t2,   // output state pair
                   double cost,      // transition cost
                   int arc_id1,      // (unused)
                   int arc_id2,      // (unused)
                   int input,        // input label
                   int intermediate, // (unused)
                   int output,       // output label
                   double base_cost, // cost of the path so far
                   int trail_index) {
            //logger.format("relaxing %d %d -> %d %d (bcost %f, cost %f)", f1, f2, t1, t2, base_cost, cost);

            if(!nbest.add_replacing_id(t1 * fst2.nStates() + t2,
                                       all_costs.length(),
                                       - base_cost - cost))
                return;

            //logger.format("nbest changed");
            //nbest.log(logger);

            if(input) {
                // The candidate for the next beam is stored in all_XX arrays.
                // (can we store it in the stree instead?)
                all_inputs.push(input);
                all_targets1.push(t1);
                all_targets2.push(t2);
                all_outputs.push(output);
                all_costs.push(cost);
                parent_trails.push(trail_index);
            } else {
                // Beam control hack
                // -----------------
                // if a node is important (changes nbest) AND its input is 0,
                // then it's added to the CURRENT beam.

                //logger.format("pushing control point from trail %d to %d, %d",
                              //trail_index, t1, t2);
                int new_node = stree.add(beam[trail_index], t1, t2, input, output, cost);
                beam.push(new_node);
                beamcost.push(base_cost + cost);

                // This is a stub entry indicating that the node should not
                // be added to the next generation beam.
                all_inputs.push(0);
                all_targets1.push(-1);
                all_targets2.push(-1);
                all_outputs.push(0);
                all_costs.push(0);
                parent_trails.push(-1);
            }
        }

        /// Call relax() for each arc going out of the given node.
        void traverse(int n1, int n2, double cost, int trail_index) {
            //logger.format("traversing %d %d", n1, n2);
            intarray &o1 = fst1.outputs(n1);
            intarray &i1 = fst1.inputs(n1);
            intarray &t1 = fst1.targets(n1);
            floatarray &c1 = fst1.costs(n1);

            intarray &o2 = fst2.outputs(n2);
            intarray &i2 = fst2.inputs(n2);
            intarray &t2 = fst2.targets(n2);
            floatarray &c2 = fst2.costs(n2);

            // for optimization
            int *O1 = o1.data;
            int *O2 = o2.data;
            int *I1 = i1.data;
            int *I2 = i2.data;
            int *T1 = t1.data;
            int *T2 = t2.data;
            float *C1 = c1.data;
            float *C2 = c2.data;

            SpecialSymbolRanges ranges1(o1);
            SpecialSymbolRanges ranges2(i2);


            // Relax outbound arcs in the composition
            int k1, k2;


            for(k1 = ranges1.start(L_SIGMA); k1 < ranges1.end(L_SIGMA); k1++) {
                if (I1[k1] == L_RHO)
                    printf("RHO:SIGMA and SIGMA:RHO arcs are not supported");

                // relaxing (_ : SIGMA) x (SIGMA  : _)
                //
                //     (SIGMA : SIGMA) x (SIGMA : SIGMA)  =  (SIGMA : SIGMA) 
                //     (SIGMA : SIGMA) x (SIGMA : a)      =  (a : a) 
                //     (  b   : SIGMA) x (SIGMA : a)      =  (b : a) 
                //     (  b   : SIGMA) x (SIGMA : SIGMA)  =  (b : b)
                for(int j=ranges2.start(L_SIGMA); j<ranges2.end(L_SIGMA); j++) {
                    int i = I1[k1] != L_SIGMA ? I1[k1] : O2[j];
                    int o = O2[j] != L_SIGMA ? O2[j] : I1[k1];

                    if (O2[j] == L_RHO)
                        printf("RHO:SIGMA and SIGMA:RHO arcs are not supported");

                    relax(n1, n2,         // from pair
                          T1[k1], T2[j],  // to pair
                          C1[k1] + C2[j], // cost
                          k1, j,         // arc ids
                          i, L_SIGMA, o,   // input, intermediate, output
                          cost, trail_index);
                }
                
                // relaxing (_ : SIGMA) x (symbol : _)
                //     (SIGMA : SIGMA) x (a : _)          =  (a : _) 
                //     (  b   : SIGMA) x (a : _)          =  (b : _)
                for(int j=ranges2.specialsEnd;j<o2.length();j++) {
                    // if it's sigma->sigma, then pick up the label,
                    // if it's x->sigma leave it alone
                    int in = (I1[k1]==L_SIGMA||I1[k1]==L_RHO)?I2[j]:I1[k1];
                    relax(n1, n2,         // from pair
                          T1[k1], T2[j],  // to pair
                          C1[k1] + C2[j], // cost
                          k1, j,         // arc ids
                          in, I2[j], O2[j],   // input, intermediate, output
                          cost, trail_index);
                }
            }

            // relaxing (_ : symbol) x (SIGMA : _)
            //
            //     (_ : a) x (SIGMA : SIGMA)  =  (_ : a) 
            //     (_ : a) x (SIGMA :   b  )  =  (_ : b)
            for(k2 = ranges2.start(L_SIGMA); k2 < ranges2.end(L_SIGMA); k2++) {
                for(int j=ranges1.specialsEnd;j<o1.length();j++) {
                    if (O2[k2] == L_RHO)
                        printf("RHO:SIGMA and SIGMA:RHO arcs are not supported");

                    // if it's sigma->sigma, then pick up the label,
                    // if it's sigma->x leave it alone
                    int out = O2[k2]==L_SIGMA ? O1[j] : O2[k2];
                    relax(n1, n2,       // from pair
                          T1[j], T2[k2],   // to pair
                          C1[j] + C2[k2],       // cost
                          j, k2,       // arc ids
                          I1[j], O1[j], out, // input, intermediate, output
                          cost, trail_index);
                }
            }

            if (ranges1.count(L_RHO) && ranges2.count(L_RHO))
                printf("composition _:RHO x RHO:_ is not supported");

            // relaxing RHO x SIGMA
            //     (   a : RHO ) x ( SIGMA : SIGMA ) - error 
            //     (   a : RHO ) x ( SIGMA : b)      = (a : b)
            //     ( RHO : RHO ) x ( SIGMA : SIGMA ) = (RHO : RHO)
            //     ( RHO : RHO ) x ( SIGMA : b)      = (RHO : b)
            //     (   _ : RHO ) x ( RHO : _ ) - error
            for(k1 = ranges1.start(L_RHO); k1 < ranges1.end(L_RHO); k1++) {
                if (I1[k1] == L_SIGMA)
                    printf("RHO:SIGMA and SIGMA:RHO arcs are not supported");

                for (int j=ranges2.start(L_SIGMA); j<ranges2.end(L_SIGMA); j++) {
                    if (I1[k1] != L_RHO && O2[j] == L_SIGMA)
                        printf("composition (a:RHO) x (SIGMA:SIGMA) is not supported");
                    int o = O2[j] == L_SIGMA ? L_RHO : O2[j];
                    relax(n1, n2,         // from pair
                          T1[k1], T2[j],  // to pair
                          C1[k1] + C2[j], // cost
                          k1, j,         // arc ids
                          I1[k1], L_RHO, o,   // input, intermediate, output
                          cost, trail_index);
                }
            }

            // relaxing SIGMA x RHO  (mirrors the previous loop)
            for(k2 = ranges2.start(L_RHO); k2 < ranges2.end(L_RHO); k2++) {
                if (O2[k2] == L_SIGMA)
                    printf("RHO:SIGMA and SIGMA:RHO arcs are not supported");

                for (int j=ranges1.start(L_SIGMA); j<ranges1.end(L_SIGMA); j++) {
                    if (O2[k2] != L_RHO && I1[j] == L_SIGMA)
                        printf("composition (SIGMA:SIGMA) x (RHO:a) is not supported");
                    int in = I1[j]==L_SIGMA ? L_RHO : I1[j];
                    relax(n1, n2,       // from pair
                          T1[j], T2[k2],   // to pair
                          C1[j] + C2[k2],       // cost
                          j, k2,       // arc ids
                          in, L_RHO, O2[k2], // input, intermediate, output
                          cost, trail_index);
                }
            }

            // relaxing fst1 EPSILON moves
            for(k1 = ranges1.start(L_EPSILON); k1 < ranges1.end(L_EPSILON); k1++) {
                relax(n1, n2,       // from pair
                      T1[k1], n2,   // to pair
                      C1[k1],       // cost
                      k1, -1,       // arc ids
                      I1[k1], L_EPSILON, L_EPSILON, // input, intermediate, output
                      cost, trail_index);
            }

            // relaxing fst2 EPSILON moves
            for(k2 = ranges2.start(L_EPSILON); k2 < ranges2.end(L_EPSILON); k2++) {
                relax(n1, n2,       // from pair
                      n1, T2[k2],   // to pair
                      C2[k2],       // cost
                      -1, k2,       // arc ids
                      L_EPSILON, L_EPSILON, O2[k2], // input, intermediate, output
                      cost, trail_index);
            }

            // The big relaxation loop. It is a merge loop for O1[k1] and I2[k2].
            k1 = ranges1.specialsEnd;
            k2 = ranges2.specialsEnd;
            while(k1 < o1.length() || k2 < i2.length()) {
                // The main matching branch (symbol to symbol)
                // (a, b) x (b, c)  =  (a, c)
                if(k1 < o1.length() && k2 < i2.length() && O1[k1] == I2[k2]) {
                    // Now we found a bunch of the same symbol on either side.
                    // We need to match every one on the left with every one on the right.
                    int next_k2;
                    for (next_k2 = k2; next_k2 < i2.length() && O1[k1] == I2[next_k2]; next_k2++) {}
                    // Now the symbol occupies the range I2[k2: next_k2]

                    for(; k1 < o1.length() && O1[k1] == I2[k2]; k1++) {
                        for(int j = k2; j < next_k2; j++) {
                            relax(n1, n2,           // from pair
                                  T1[k1], T2[j],    // to pair
                                  C1[k1] + C2[j],   // cost
                                  k1, j,            // arc ids
                                  I1[k1], O1[k1], O2[j], // input, intermediate, output
                                  cost, trail_index);
                        }
                    }
                    k2 = next_k2;
                }

                // The following part of the big loop deals with RHOs and PHIs.

                // Move only k1 while we can
                // Each iteration we have a letter in fst1 
                //    that matches no letters on the other side
                int max_o1_k1 = k2 < i2.length() ? I2[k2] - 1 : INT_MAX;
                for(; k1 < o1.length() && O1[k1] <= max_o1_k1; k1++) {
                    // Relax the fst1's arc against rho
                    for(int j=ranges2.start(L_RHO); j<ranges2.end(L_RHO); j++) {
                        // (_ : a) x (RHO : SIGMA)  =  (_ : a)
                        // (_ : a) x (RHO :  RHO )  =  (_ : a)
                        // (_ : a) x (RHO :   b  )  =  (_ : b)
                        int out = (O2[j] == L_SIGMA || O2[j] == L_RHO)
                                     ? O1[k1] : O2[j];
                        relax(n1, n2,           // from pair
                              T1[k1], T2[j],   // to pair
                              C1[k1] + C2[j],  // cost
                              k1, j,           // arc ids
                              I1[k1], O1[k1], out,     // input, intermediate, output
                              cost, trail_index);                      
                    }
                    // An EPSILON is produced for every character in fst1
                    // that matched nothing in fst2.
                    // XXX: this is inefficient (produces too many epsilons)
                    for(int j=ranges2.start(L_PHI); j<ranges2.end(L_PHI); j++) {
                        relax(n1, n2,      // from pair
                              n1, T2[j],   // to pair
                              C2[j],       // cost
                              -1, j,       // arc ids
                              L_EPSILON, L_EPSILON, O2[j],     // input, intermediate, output
                              cost, trail_index);                      
                    }
                } // rho/phi loop, left side

                // This loop mirrors the previous one.
                // Every iteration we have a letter in fst2
                //    that didn't match in fst1
                int max_i2_k2 = k1 < i1.length() ? O1[k1] - 1 : INT_MAX;
                for(;k2 < i2.length() && I2[k2] <= max_i2_k2; k2++) { 
                    // Relax the fst2's arc against rho and phi
                    for(int j=ranges1.start(L_RHO); j<ranges1.end(L_RHO); j++) {
                        int in = (I1[j] == L_SIGMA || I1[j] == L_RHO)
                                    ? I2[k2] : I1[j];
                        relax(n1, n2,           // from pair
                              T1[j], T2[k2],   // to pair
                              C1[j] + C2[k2],  // cost
                              j, k2,           // arc ids
                              in, I2[k2], O2[k2],     // input, intermediate, output
                              cost, trail_index);                      
                    }
                    for(int j=ranges1.start(L_PHI); j<ranges1.end(L_PHI); j++) {
                        relax(n1, n2,      // from pair
                              T1[j], T2[k2],   // to pair
                              C1[j],       // cost
                              j, -1,       // arc ids
                              I1[j], L_EPSILON, L_EPSILON,     // input, intermediate, output
                              cost, trail_index);                      
                    }
                } // rho/phi loop, right side
            } // the merge loop
        }

        // The main loop iteration.
        void radiate() {
            clear();

            //logger("beam", beam);
            //logger("beamcost", beamcost);

            int control_beam_start = beam.length();
            for(int i = 0; i < control_beam_start; i++)
                try_accept(i);

            // in this loop, traversal may add "control nodes" to the beam
            for(int i = 0; i < beam.length(); i++) {
                traverse(stree.v1[beam[i]], stree.v2[beam[i]],
                         beamcost[i], i);
            }

            // try accepts from control beam nodes
            // (they're not going to the next beam)
            for(int i = control_beam_start; i < beam.length(); i++)
                try_accept(i);


            intarray new_beam;
            floatarray new_beamcost;
            for(int i = 0; i < nbest.length(); i++) {
                int k = nbest.tag(i);
                if(parent_trails[k] < 0) // skip the control beam nodes
                    continue;
                new_beam.push(stree.add(beam[parent_trails[k]],
                                        all_targets1[k], all_targets2[k],
                                        all_inputs[k], all_outputs[k],
                                        all_costs[k]));
                new_beamcost.push(beamcost[parent_trails[k]] + all_costs[k]);
                //logger.format("to new beam: trail index %d, stree %d, target %d,%d",
                        //k, new_beam[new_beam.length() - 1], all_targets1[k], all_targets2[k]);
            }
            move(beam, new_beam);
            move(beamcost, new_beamcost);
        }

        // Relax the accept arc from the beam node number i.
        void try_accept(int i) {
            float a_cost1 = fst1.getAcceptCost(stree.v1[beam[i]]);
            float a_cost2 = fst2.getAcceptCost(stree.v2[beam[i]]);
            float candidate = beamcost[i] + a_cost1 + a_cost2;
            if(candidate < best_cost_so_far) {
                //logger.format("accept from beam #%d (stree %d), cost %f",
                //              i, beam[i], candidate);
                best_so_far = beam[i];
                best_cost_so_far = candidate;
            }
        }

        void bestpath(intarray &v1, intarray &v2, intarray &inputs,
                      intarray &outputs, floatarray &costs) {
            stree.clear();

            beam.resize(1);
            beamcost.resize(1);
            beam[0] = stree.add(-1, fst1.getStart(), fst2.getStart(), 0, 0, 0);
            beamcost[0] = 0;

            best_so_far = 0;
            best_cost_so_far = fst1.getAcceptCost(fst1.getStart()) +
                               fst2.getAcceptCost(fst1.getStart());

            while(beam.length())
                radiate();

            stree.get(v1, v2, inputs, outputs, costs, best_so_far);
            costs.push(fst1.getAcceptCost(stree.v1[best_so_far]) +
                       fst2.getAcceptCost(stree.v2[best_so_far]));

            //logger("costs", costs);
        }
    };

};

namespace ocrofst{
    void beam_search(intarray &vertices1,
                     intarray &vertices2,
                     intarray &inputs,
                     intarray &outputs,
                     floatarray &costs,
                     OcroFST &fst1,
                     OcroFST &fst2,
                     int beam_width) {
        BeamSearch b(fst1, fst2, beam_width);
        CHECK(L_SIGMA<L_EPSILON);
        CHECK(L_RHO<L_PHI);
        CHECK(L_PHI<L_EPSILON);
        CHECK(L_EPSILON<1);
        fst1.sortByOutput();
        fst2.sortByInput();
        //fprintf(stderr,"starting bestpath\n");
        b.bestpath(vertices1, vertices2, inputs, outputs, costs);
        //fprintf(stderr,"finished bestpath\n");
    }

    double beam_search(ustrg &result, OcroFST &fst1, OcroFST &fst2,
                       int beam_width) {
        intarray v1;
        intarray v2;
        intarray i;
        intarray o;
        floatarray c;
        //fprintf(stderr,"starting beam search\n");
        beam_search(v1, v2, i, o, c, fst1, fst2, beam_width);
        //fprintf(stderr,"finished beam search\n");
        remove_epsilons(result, o);
        return sum(c);
    }
}
