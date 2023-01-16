//
// EventTracker.h
//
// Class for better handling of UseEventRange and UseTimeRange
// L GOLINO
//

#include "EventTracker.h"




EventTracker::EventTracker() // ctor
{
    fIDCut = false;
    fTimeCut = false;
}

EventTracker::EventTracker(std::string fileName, int runNumber) // ctor
{
    fIDCut = true;
    fRunNumber = runNumber;

    LoadEventIDs(fileName);
}

void EventTracker::AddEventRange(int start_event, int stop_event)
{
    fIDCut = true;
    fEventIDs.push_back( std::pair<int, int>(start_event, stop_event) );
}

void EventTracker::AddTimeRange(double star_event, double stop_event)
{
    fTimeCut = true;
    fEventTimes.push_back( std::pair<double, double>(star_event, stop_event) );
}

void EventTracker::LoadEventIDs(std::string fileName)
{
    std::ifstream myFile(fileName);
    std::string line;
    if(myFile.fail()){
        std::cout << "WARNING: Eventlist failed to open. Likely doesn't exist but could be a permission error.\n";
        exit(-1);
    }
    if (myFile.is_open())
    {
        while ( std::getline(myFile, line) )
        {
            //std::cout << line << std::endl;
            int runNumber, firstEvent, lastEvent;
            //The following is apparantly best method: https://quick-bench.com/q/CWCbHcvWTZBXydPA_mju2r75LX0
            if (3 == std::sscanf(line.c_str(), "%d:%d-%d", &runNumber, &firstEvent, &lastEvent))
            {
                if(runNumber == fRunNumber)
                {
                    std::cout << "runNumber=" << runNumber << std::endl;
                    std::cout << "firstEvent=" << firstEvent << std::endl;
                    std::cout << "lastEvent=" << lastEvent << std::endl;
                }
            }
            if(runNumber == fRunNumber)
                fEventIDs.push_back( std::pair<int, int>(firstEvent, lastEvent) );
        }
        myFile.close();
    }
}

void EventTracker::SortDeques()
{
    //Sort based on first (tmin) keeping idx in the same location
    std::sort(std::begin(fEventIDs), std::end(fEventIDs), 
    [&](const std::pair<int,int>& lhs, const std::pair<int,int>& rhs)
    {
        if(lhs.first == rhs.first)
            return lhs.second < rhs.second;
        else
            return lhs.first < rhs.first;
    } );

    //Sort based on first (tmin) keeping idx in the same location
    std::sort(std::begin(fEventTimes), std::end(fEventTimes), 
    [&](const std::pair<double,double>& lhs, const std::pair<double,double>& rhs)
    {
        if(lhs.first == rhs.first)
            return lhs.second < rhs.second;
        else
            return lhs.first < rhs.first;
    } );
}

EventTracker::~EventTracker() // dtor
{
}


bool EventTracker::IsEventInRange(int eventID, double eventTime)
{

    if(fRejectAll)
        return false;

    bool passidcut = !fIDCut;
    if(fIDCut)
    {
        //If something happens and size is 0 when fIDCut is true- start rejecting all
        if(fEventIDs.size() == 0)
        {
            std::cout << "WARNING: Not reconstructing events from now on. Did you name a file that doesn't exist on the --eventlist flag?\n";
            fRejectAll = true; //Reject all from now on. 
            return false;
        }
        //Are we out of range of first pair?
        if(eventID>fEventIDs.at(0).second)
        {
            fEventIDs.pop_front(); //If so we've passed this interval, feel free to remove.
            if(fEventIDs.size() == 0)
            {
                fRejectAll = true; //Reject all from now on. 
                return false;
            }
        }
        if(eventID>=fEventIDs.at(0).first && eventID<=fEventIDs.at(0).second) //If we are in range reconstruct this event.
        {
            passidcut = true;
        }
    }

    bool passtimecut = !fTimeCut;
    if(fTimeCut)
    {
        //Are we out of range of first pair?
        if(eventTime>fEventTimes.at(0).second)
        {
            fEventTimes.pop_front(); //If so we've passed this interval, feel free to remove.
            if(fEventTimes.size() == 0)
            {
                fRejectAll = true; //Reject all from now on. 
                return false;
            }
        }
        if(eventTime>=fEventTimes.at(0).first && eventTime<=fEventTimes.at(0).second) //If we are in range reconstruct this event.
        {
            passtimecut = true;
        }

    }

    //It was in none of the event range or time ranges, so reject 
    if(passidcut && passtimecut) //if you would rather it be a pass if its in EITHER range: change this to a || 
        return true;
    else
        return false;
}

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */