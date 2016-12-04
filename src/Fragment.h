/*
 * Fragment.h
 *
 *  Created on: Jan 7, 2013
 *      Author: msegar
 */

#ifndef FRAGMENT_H_
#define FRAGMENT_H_

#include <string>

using namespace std;

class Fragment {
private:
    string sAnchoredRead;
    int64_t iAnchoredStart;
    int64_t iFlag;
    string sUnanchoredRead;
    int64_t iUnanchoredStart;
    string sParentRead;
    int64_t iParentStart;
    bool bAnchorLeft;
public:
    Fragment();
    Fragment(string anchor, int64_t anStart, int64_t flag, bool left, string unanchor, int64_t unStart);
    Fragment(string anchor, int64_t anStart, int64_t flag, bool left, string unanchor, int64_t unStart, string parent, int64_t pStart);

    virtual ~Fragment();

    // get functions
    bool isAnchorLeft()        { return bAnchorLeft; }
    int64_t getAnchoredStart()     { return iAnchoredStart; }
    int64_t getFlag()            { return iFlag; }
    int64_t getUnanchoredStart() { return iUnanchoredStart; }
    bool getAnchorLeft()        { return bAnchorLeft; }
    int64_t getParentStart()     { return iParentStart; }
    string* getAnchoredRead() { return &sAnchoredRead; }
    string* getParentRead() { return &sParentRead; }
    string* getUnanchoredRead() { return &sUnanchoredRead; }

    // set functions
    void setParentRead(string parentRead) { sParentRead = parentRead; }
    void setParentStart(int parentStart) { iParentStart = parentStart; }
    //void setUnanchoredRead(string unanchoredRead) { sUnanchoredRead = unanchoredRead; }
};

#endif /* FRAGMENT_H_ */
