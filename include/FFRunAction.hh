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
//! \file       FFRunAction.hh
//! \author     Adapted for U-235 criticality simulation
//! \date       September 10, 2025
//!
//! \brief      Definition of the FFRunAction class
//!
//! \details    Computes criticality parameters using Geant4 histograms
//!
// -------------------------------------------------------------

#ifndef FFRUNACTION
#define FFRUNACTION

#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;

class FFRunAction : public G4UserRunAction
{
public:
    FFRunAction();
    virtual ~FFRunAction();
    virtual void BeginOfRunAction(const G4Run* run);
    virtual void EndOfRunAction(const G4Run* run);
};

#endif // FFRUNACTION
