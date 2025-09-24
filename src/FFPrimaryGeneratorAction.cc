#include "globals.hh"
#include "G4Event.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Neutron.hh"
#include "G4ParticleGun.hh"
#include "G4Sphere.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "FFPrimaryGeneratorAction.hh"
#include "G4RandomTools.hh"

// Constructor
FFPrimaryGeneratorAction::FFPrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
#ifndef NDEBUG
  fEventNumber(0),
#endif // NDEBUG
  fParticleGun(new G4ParticleGun(1)),
  fSphereSolid(nullptr),
  fMaxVal(0.0)
{
    fParticleGun->SetParticleDefinition(G4Neutron::Definition());
    // Precompute max for Watt spectrum
    G4double a = 0.988 * MeV;
    G4double b = 2.249 / MeV;
    G4double E_max = 20.0 * MeV;
    G4double de = 0.01 * MeV;
    for (G4double e = 0; e < E_max; e += de) {
        G4double val = std::exp(-e / a) * std::sinh(std::sqrt(b * e));
        if (val > fMaxVal) fMaxVal = val;
    }
}

// Destructor
FFPrimaryGeneratorAction::~FFPrimaryGeneratorAction()
{
    delete fParticleGun;
}

// Sample Watt spectrum
G4double FFPrimaryGeneratorAction::SampleWattEnergy()
{
    G4double a = 0.988 * MeV;
    G4double b = 2.249 / MeV;
    G4double E_max = 20.0 * MeV;
    G4double E;
    do {
        E = -a * std::log(G4UniformRand()) * (1.0 - std::exp(-E_max / a));
        if (E > E_max) continue;
    } while (G4UniformRand() > std::exp(-E / a) * std::sinh(std::sqrt(b * E)) / fMaxVal);
    return E;
}

// Generate primaries
void FFPrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
#ifndef NDEBUG
    G4cout << "Shooting event " << ++fEventNumber << G4endl;
#endif // NDEBUG

    // Get the U-235 sphere geometry
    if (fSphereSolid == nullptr) {
        G4LogicalVolume* temp = G4LogicalVolumeStore::GetInstance()->GetVolume("U235Sphere");
        if (temp != nullptr) {
            fSphereSolid = dynamic_cast<G4Sphere*>(temp->GetSolid());
        }
        if (fSphereSolid == nullptr) {
            G4Exception("FFPrimaryGeneratorAction::GeneratePrimaries",
                        "U235 sphere solid volume not found",
                        FatalException,
                        "This run will be aborted");
        }
    }

    // Position at center
    G4ThreeVector position(0.0, 0.0, 0.0);

#ifndef NDEBUG
    G4cout << "Emission Location: " << position/cm << " cm" << G4endl;
#endif // NDEBUG

    // Sample isotropic emission direction
    G4ThreeVector direction = G4RandomDirection();

#ifndef NDEBUG
    G4cout << "Emission Direction: " << direction << G4endl;
#endif // NDEBUG

    // Sample energy from Watt spectrum
    G4double energy = SampleWattEnergy();

    // Set particle properties and generate
    fParticleGun->SetParticleEnergy(energy);
    fParticleGun->SetParticlePosition(position);
    fParticleGun->SetParticleMomentumDirection(direction);
    fParticleGun->GeneratePrimaryVertex(event);
}

