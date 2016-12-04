/*
 * Fragment.cpp
 *
 *  Created on: Jan 7, 2013
 *      Author: msegar
 */

#include "Fragment.h"

using namespace std;

Fragment::Fragment() {
    sAnchoredRead = "";
    iAnchoredStart = -1;
    iFlag = -1;
    sUnanchoredRead = "";
    iUnanchoredStart = -1;
    sParentRead = "";
    iParentStart = -1;
    bAnchorLeft = true;
}

Fragment::Fragment(string anchor, int64_t anStart, int64_t flag, bool left, string unanchor, int64_t unStart){
    sAnchoredRead = anchor;
    iAnchoredStart = anStart;
    iFlag = flag;
    sUnanchoredRead = unanchor;
    iUnanchoredStart = unStart;
    sParentRead = "";
    iParentStart = -1;
    bAnchorLeft = left;
}

Fragment::Fragment(string anchor, int64_t anStart, int64_t flag, bool left, string unanchor, int64_t unStart, string parent, int64_t pStart){
    sAnchoredRead = anchor;
    iAnchoredStart = anStart;
    iFlag = flag;
    sUnanchoredRead = unanchor;
    iUnanchoredStart = unStart;
    sParentRead = parent;
    iParentStart = pStart;
    bAnchorLeft = left;
}


Fragment::~Fragment() {
    // Auto-generated destructor stub
}
