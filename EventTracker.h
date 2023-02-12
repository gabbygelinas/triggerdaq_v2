//
// EventTracker.h
//
// Class for better handling of UseEventRange and UseTimeRange
// L GOLINO
//

#ifndef EventTracker_H
#define EventTracker_H

#include<deque>
#include<string>
#include<fstream>
#include<iostream>
#include<algorithm>

class EventTracker
{
    private:
    //Bools to check whether user has requested a cut on either of these
    bool fIDCut = false;
    bool fTimeCut = false;
    bool fRejectAll = false;

    int fRunNumber;

    //deque's to reflect the range of ID's or times.
    //Both in format <x1, y1, x2, y2, ..., xn, yn} (we use deque to pop front when we are done with range {x1, y1})
    std::deque<std::pair<int,int>> fEventIDs;   
    std::deque<std::pair<double,double>> fEventTimes;   
    
    public:
    //Function to find whether an event is in range (and therefore reconstructible)
    bool IsEventInRange(int eventID, double eventTime);

    void LoadEventIDs(std::string fileName);
    void SortDeques();

    void AddEventRange(int start_event, int stop_event);
    void AddTimeRange(double start_event, double stop_event);

    bool GetEventCut()      { return fIDCut; }
    bool GetTimeCut()       { return fTimeCut; }
    bool GetRunNumber()     { return fRunNumber; }


    EventTracker(); // ctor
    EventTracker(std::string fileName, int runNumber); // ctor
    ~EventTracker(); // dtor
   
};


#endif

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */