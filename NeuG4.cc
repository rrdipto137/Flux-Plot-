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
//! \file       NeuG4.cc
//! \author     Adapted for criticality simulation
//! \date       September 10, 2025
//!
//! \brief      Main program for criticality simulation
//!
// -------------------------------------------------------------

#include "G4MTRunManager.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "FTFP_BERT_HP.hh"
#include "FFDetectorConstruction.hh"
#include "FFActionInitialization.hh"
#include "FFPhysicsList.hh"
int main(int argc, char** argv)
{
    // Instantiate the multi-threaded run manager
    G4MTRunManager* runManager = new G4MTRunManager;

    // Set the number of threads (adjust based on your system)
    runManager->SetNumberOfThreads(11);

    // Set mandatory initialization classes
    runManager->SetUserInitialization(new FFDetectorConstruction());
    runManager->SetUserInitialization(new FTFP_BERT_HP());
    runManager->SetUserInitialization(new FFActionInitialization());
runManager->SetUserInitialization(new FFPhysicsList());
    // Initialize G4 kernel
    runManager->Initialize();

    // Get the pointer to the UI manager
    G4UImanager* UImanager = G4UImanager::GetUIpointer();

    // Process macro file or start interactive session
    if (argc != 1) {
        // Batch mode
        G4String command = "/control/execute ";
        G4String fileName = argv[1];
        UImanager->ApplyCommand(command + fileName);
    } else {
        // Interactive mode
        G4VisExecutive* visManager = new G4VisExecutive;
        visManager->Initialize();
        G4UIExecutive* ui = new G4UIExecutive(argc, argv);
        UImanager->ApplyCommand("/control/execute init_vis.mac");
        ui->SessionStart();
        delete ui;
        delete visManager;
    }

    // Clean up
    delete runManager;

    return 0;
}
