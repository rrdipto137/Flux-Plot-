f#include "globals.hh"
#include "G4Box.hh"
#include "G4Element.hh"
#include "G4Isotope.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4Sphere.hh"
#include "G4SystemOfUnits.hh"
#include "G4VPhysicalVolume.hh"
#include "FFDetectorConstruction.hh"

static const G4double inch = 2.54 * cm;

// Constructor
FFDetectorConstruction::FFDetectorConstruction()
: G4VUserDetectorConstruction()
{
    DefineMaterials();
}

// Destructor
FFDetectorConstruction::~FFDetectorConstruction()
{
    // Nothing to clean up
}

// Define materials
void FFDetectorConstruction::DefineMaterials(void)
{
    G4NistManager* nist = G4NistManager::Instance();
    
    // Define vacuum for the world
    fVacuum = new G4Material("Vacuum", 1.e-25 * g/cm3, 1);
    fVacuum->AddElement(nist->FindOrBuildElement("N"), 1.0);
    
   // Define HEU for Godiva-like
    //G4Isotope* iU234 = new G4Isotope("U234", 92, 234, 234.040952 * g/mole);
    G4Isotope* iU235 = new G4Isotope("U235", 92, 235, 235.043930 * g/mole);
    G4Isotope* iU238 = new G4Isotope("U238", 92, 238, 238.050788 * g/mole);
    G4Element* elU = new G4Element("EnrichedUranium", "EU", 2);
    //elU->AddIsotope(iU234, 1.0252 * perCent);
    elU->AddIsotope(iU235, 90 * perCent);
    elU->AddIsotope(iU238, 10 * perCent);
    fU235 = new G4Material("HEU", 18.7398 * g/cm3, 1);
    fU235->AddElement(elU, 1.0);
    
    
    
}

// Construct geometry
G4VPhysicalVolume* FFDetectorConstruction::Construct()
{
    G4ThreeVector position;
#ifdef NDEBUG
    G4bool const overlapChecking = false;
#else
    G4bool const overlapChecking = true;
#endif // NDEBUG
    
    // Create the world (vacuum)
    const G4double worldSize = 100.0 * cm; // Large enough to avoid boundary effects
    G4Box* solidWorld = new G4Box("World", worldSize, worldSize, worldSize);
    G4LogicalVolume* logicalWorld = new G4LogicalVolume(solidWorld, fVacuum, "World");
    position.set(0.0, 0.0, 0.0);
    G4VPhysicalVolume* physicalWorld = new G4PVPlacement(nullptr, position, logicalWorld,
                                                        "World", nullptr, false, 0, overlapChecking);
    
    // Create U-235 sphere
    const G4double sphereRadius =3* cm; // Godiva radius
    G4Sphere* solidSphere = new G4Sphere("U235Sphere", 0.0, sphereRadius,
                                        0.0 * deg, 360.0 * deg, 0.0 * deg, 180.0 * deg);
    G4LogicalVolume* logicalSphere = new G4LogicalVolume(solidSphere, fU235, "U235Sphere");
    position.set(0.0, 0.0, 0.0);
    new G4PVPlacement(nullptr, position, logicalSphere, "U235Sphere",
                      logicalWorld, false, 0, overlapChecking);
    
    return physicalWorld;
}

// Remove PlaceFuelPlate as it's no longer needed
void FFDetectorConstruction::PlaceFuelPlate(double, double, G4LogicalVolume*, G4LogicalVolume*)
{
    // Empty: not used in this simulation
}
