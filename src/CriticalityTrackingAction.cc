
//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications, and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// -------------------------------------------------------------
//!
//! \file       CriticalityTrackingAction.cc
//! \author     Adapted for U-235 criticality simulation
//! \date       September 10, 2025
//!
//! \brief      Implementation of the CriticalityTrackingAction class
//!
//! \details    Records birth and death times for population, and removal times
//!
// -------------------------------------------------------------

#include "CriticalityTrackingAction.hh"
#include "G4Track.hh"
#include "G4Neutron.hh"
#include "G4AnalysisManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4VProcess.hh"

// Constructor
CriticalityTrackingAction::CriticalityTrackingAction()
: G4UserTrackingAction()
{
}

// Destructor
CriticalityTrackingAction::~CriticalityTrackingAction()
{
}

// Pre tracking
void CriticalityTrackingAction::PreUserTrackingAction(const G4Track* track)
{
    if (track->GetDefinition() != G4Neutron::Definition()) return;

    G4double time = track->GetGlobalTime();
    G4AnalysisManager* man = G4AnalysisManager::Instance();
    man->FillH1(0, time / ns, 1.0); // birth_hist id=0
}

// Post tracking
void CriticalityTrackingAction::PostUserTrackingAction(const G4Track* track)
{
    if (track->GetDefinition() != G4Neutron::Definition()) return;

    G4double time = track->GetGlobalTime();
    G4AnalysisManager* man = G4AnalysisManager::Instance();
    man->FillH1(1, time / ns, 1.0); // death_hist id=1

    // For removal time
    const G4VProcess* creator = track->GetCreatorProcess();
    G4String creatorName = (creator) ? creator->GetProcessName() : "";
    G4String termName = (track->GetStep()->GetPostStepPoint()->GetProcessDefinedStep()) ? track->GetStep()->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() : "";

    if ((creatorName == "nFission" || creatorName == "hadInelastic") &&
        (termName == "nCapture" || termName == "nFission" || termName == "hadInelastic" || termName == "Transportation")) {
        man->FillH1(4, track->GetLocalTime() / ns, 1.0); // tau_hist id=4
    }
}

