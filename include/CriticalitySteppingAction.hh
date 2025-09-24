#ifndef CriticalitySteppingAction_h
#define CriticalitySteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"

class G4Step;

class CriticalitySteppingAction : public G4UserSteppingAction
{
public:
    CriticalitySteppingAction();
    virtual ~CriticalitySteppingAction();

    virtual void UserSteppingAction(const G4Step*);

private:
    G4int fTrackLengthId; // Member variable to store track length histogram ID
};

#endif
