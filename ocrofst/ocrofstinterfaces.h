#ifndef h_ocrofstinterfaces__
#define h_ocrofstinterfaces__

// Copyright 2008 Deutsches Forschungszentrum fuer Kuenstliche Intelligenz
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
// File: ocr-pfst-IGenericFST.h
// Purpose: OpenFST-compatible I/O
// Responsible: mezhirov
// Reviewer: 
// Primary Repository: 
// Web Sites: www.iupr.org, www.dfki.de, www.ocropus.org


#include <stdio.h>
#include <stdlib.h>
#include <typeinfo>
#include <colib/narray.h>
#include <colib/narray-util.h>
#include <colib/smartptr.h>
#include <colib/misc.h>
#include <colib/coords.h>
#include <colib/iustring.h>
#include "iulib/components.h"


namespace ocrofst{
	
    /// \brief A generic interface for language models.

    /// An IGenericFst is a directed graph
    /// with output/cost/id written on arcs,
    /// accept cost written on vertices and
    /// a fixed start vertice.
    struct IGenericFst : virtual IComponent {
        const char *interface() { return "IGenericFst"; }
        /// Clear the language model
        virtual void clear() { throw Unimplemented(); }

        /// Get a single new state
        virtual int newState() { throw Unimplemented(); }

        /// Add a transition between the given states
        virtual void addTransition(int from,int to,int output,float cost,int input) { throw Unimplemented(); }

        /// A variant of addTransition() with equal input and output.
        virtual void addTransition(int from,int to,int symbol,float cost) {
            addTransition(from, to, symbol, cost, symbol);
        }

        /// Set the start state
        virtual void setStart(int node) { throw Unimplemented(); }

        /// Set a state as an accept state
        virtual void setAccept(int node,float cost=0.0) { throw Unimplemented(); }

        /// Obtain codes for "specials" (language model dependent)
        virtual int special(const char *s) { throw Unimplemented(); }

        /// \brief Compute the best path through the language model.
        /// Useful for simple OCR tasks and for debugging.
        virtual void bestpath(ustrg &result) { throw Unimplemented(); }

        /// destroy the language model
        virtual ~IGenericFst() {}

        /// simple interface for line recognizers
        virtual void setString(ustrg &text,floatarray &costs,intarray &ids) {
            int n = text.length();
            intarray states;
            states.clear();
            for(int i=0;i<n+1;i++)
                states.push(newState());
            for(int i=0;i<n;i++)
                addTransition(states[i],states[i+1],text[i].ord(),costs[i],ids[i]);
            setStart(states[0]);
            setAccept(states[n]);
        }

        // reading methods

        /// Get the number of states.
        virtual int nStates() { throw Unimplemented(); }

        /// Get the starting state.
        virtual int getStart() { throw Unimplemented(); }

        /// Get the accept cost of a given vertex (a cost to finish the line and quit).
        virtual float getAcceptCost(int node) { throw Unimplemented(); }

        /// Determine whether the given node is an accepting state.
        virtual bool isAccepting(int node) { return getAcceptCost(node)<1e30; }

        /// Return an array of arcs leading from the given node.
        virtual void arcs(colib::intarray &ids,
                          colib::intarray &targets,
                          colib::intarray &outputs,
                          colib::floatarray &costs,
                          int from) { throw Unimplemented(); } // WARN_DEPRECATED

        /// A variant of addTransition() with equal input and output.
        virtual void getTransitions(intarray &tos,intarray &symbols,floatarray &costs,intarray &inputs,int from) {
            arcs(inputs,tos,symbols,costs,from);
        }

        /// Change a transition score between the given states
        virtual void rescore(int from,int to,int output,float new_cost,int input) { throw Unimplemented(); }

        /// A variant of rescore() with equal input and output.
        virtual void rescore(int from, int to, int symbol, float new_cost) {
            rescore(from, to, symbol, new_cost, symbol);
        }

        /// These methods should load and save in OpenFST format.
        /// (A simple way of doing that is to convert internally to OpenFST,
        /// then call its load/save methods.)
        virtual void load(const char *file) { throw Unimplemented(); }
        virtual void save(const char *file) { throw Unimplemented(); }
    };

}
#endif
