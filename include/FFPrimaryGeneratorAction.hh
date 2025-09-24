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
//  =============== Begin Documentation Comments ===============
//!
//! \file       FFPrimaryGeneratorAction.hh
//! \author     B. Wendt (brycen.linn.wendt@cern.ch)
//! \date       June 06, 2014
//!
//! \brief      Definition of the FFPrimaryGeneratorAction class
//!
//! \details    Generates neutrons with Watt spectrum at the center of the U-235 sphere
//!
//  ================ End Documentation Comments ================
//
//  Modified:
//  2025-09-10    Updated for U-235 sphere criticality simulation, added Watt spectrum sampling
//
// -------------------------------------------------------------

#ifndef FFPRIMARYGENERATORACTION
#define FFPRIMARYGENERATORACTION

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4Sphere.hh"
#include "G4VUserPrimaryGeneratorAction.hh"

class FFPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
    // Constructor
    FFPrimaryGeneratorAction();

    // Functions
    virtual void GeneratePrimaries(G4Event* event);
    G4double SampleWattEnergy();

    // Destructor
    virtual ~FFPrimaryGeneratorAction();

private:
    // Fields
#ifndef NDEBUG
    G4long fEventNumber;
#endif // NDEBUG
    G4ParticleGun* const fParticleGun;
    G4Sphere* fSphereSolid;
    G4double fMaxVal;
};

#endif //FFPRIMARYGENERATORACTION
