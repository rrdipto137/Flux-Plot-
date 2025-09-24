#include "CriticalitySteppingAction.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4Neutron.hh"
#include "G4VProcess.hh"
#include "G4AnalysisManager.hh"
#include "G4SystemOfUnits.hh"
#include <cmath>

// Constructor
CriticalitySteppingAction::CriticalitySteppingAction()
: G4UserSteppingAction(), fTrackLengthId(-1)
{
}

// Destructor
CriticalitySteppingAction::~CriticalitySteppingAction()
{
}

// Process each step
void CriticalitySteppingAction::UserSteppingAction(const G4Step* step)
{
    G4Track* track = step->GetTrack();
    if (track->GetDefinition() != G4Neutron::Definition()) return;

    G4StepPoint* postPoint = step->GetPostStepPoint();
    G4String process = postPoint->GetProcessDefinedStep()->GetProcessName();
    G4double time = postPoint->GetGlobalTime();
    G4AnalysisManager* man = G4AnalysisManager::Instance();

    G4String preVol = step->GetPreStepPoint()->GetPhysicalVolume()->GetName();
    G4String postVol = (postPoint->GetPhysicalVolume()) ? postPoint->GetPhysicalVolume()->GetName() : "";

    // Flux tally (track length estimator)
    if (preVol == "U235Sphere") {
        G4double ekin = step->GetPreStepPoint()->GetKineticEnergy();
        if (ekin > 0) {
            G4double logE = std::log10(ekin / eV);
            G4double length = step->GetStepLength() / cm;
            man->FillH1(5, logE, length);
        }

        // Spatial track length for flux map
        if (fTrackLengthId < 0) {
            fTrackLengthId = man->GetH3Id("track_length");
        }
        G4ThreeVector prePos = step->GetPreStepPoint()->GetPosition();
        G4ThreeVector postPos = step->GetPostStepPoint()->GetPosition();
        G4ThreeVector pos = (prePos + postPos) / 2.0 / cm;
        G4double length_cm = step->GetStepLength() / cm;
        man->FillH3(fTrackLengthId, pos.x(), pos.y(), pos.z(), length_cm);
    }

    // Fission sites for Shannon entropy
    if (process == "nFission") {
        G4ThreeVector pos = postPoint->GetPosition() / cm;
        man->FillH3(0, pos.x(), pos.y(), pos.z(), 1.0);
    }

    // Production: neutrons from fission or inelastic
    if (process == "nFission" || process == "hadInelastic") {
        const G4TrackVector* secondaries = step->GetSecondary();
        G4int numNeut = 0;
        for (const auto* sec : *secondaries) {
            if (sec->GetDefinition() == G4Neutron::Definition()) numNeut++;
        }
        if (numNeut > 0) {
            man->FillH1(2, time / ns, static_cast<G4double>(numNeut)); // prod_hist id=2
        }
    }

    // Loss: absorption (if killed)
    if ((process == "nCapture" || process == "nFission" || process == "hadInelastic") && track->GetTrackStatus() == fStopAndKill) {
        man->FillH1(3, time / ns, 1.0); // loss_hist id=3
    }

    // Loss: escape from sphere to world
    if (preVol == "U235Sphere" && postVol == "World") {
        man->FillH1(3, time / ns, 1.0); // loss_hist
        track->SetTrackStatus(fStopAndKill);
    }
}
