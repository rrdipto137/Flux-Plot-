#include "FFRunAction.hh"
#include "G4Run.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"
#include "G4Threading.hh"
#include "G4AnalysisManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Sphere.hh"
#include <iomanip>
#include <cmath>

// Constructor
FFRunAction::FFRunAction()
: G4UserRunAction()
{
}

// Destructor
FFRunAction::~FFRunAction()
{
}

// Begin of run action
void FFRunAction::BeginOfRunAction(const G4Run* /*run*/)
{
    G4AnalysisManager* man = G4AnalysisManager::Instance();

    // Get the sphere radius dynamically
    G4LogicalVolumeStore* store = G4LogicalVolumeStore::GetInstance();
    G4LogicalVolume* log_sphere = store->GetVolume("U235Sphere", false);
    G4double radius = 8.7; // Default radius in cm if not found
    if (log_sphere) {
        G4VSolid* solid = log_sphere->GetSolid();
        G4Sphere* solid_sphere = dynamic_cast<G4Sphere*>(solid);
        if (solid_sphere) {
            radius = solid_sphere->GetOuterRadius() / cm;
        } else {
            G4cout << "Warning: U235Sphere solid is not a G4Sphere, using default radius." << G4endl;
        }
    } else {
        G4cout << "Warning: U235Sphere logical volume not found, using default radius." << G4endl;
    }

    man->SetDefaultFileType("root");
    man->OpenFile("criticality");

    man->CreateH1("birth", "Neutron birth times (ns)", 300, 0., 300.);
    man->CreateH1("death", "Neutron death times (ns)", 300, 0., 300.);
    man->CreateH1("prod", "Production times (ns)", 300, 0., 300.);
    man->CreateH1("loss", "Loss times (ns)", 300, 0., 300.);
    man->CreateH1("tau", "Removal times (ns)", 100, 0., 50.);
    man->CreateH1("flux_logE", "Track length vs log10(E/eV)", 500, -3., 7.3);
    man->CreateH3("fission_sites", "Fission sites (cm)", 20, -radius, radius, 20, -radius, radius, 20, -radius, radius);
    man->CreateH3("track_length", "Track length (cm)", 20, -radius, radius, 20, -radius, radius, 20, -radius, radius);
}

// End of run action
void FFRunAction::EndOfRunAction(const G4Run* run)
{
    G4AnalysisManager* man = G4AnalysisManager::Instance();

    if (G4Threading::IsMasterThread()) {
        // Get the sphere radius and volume dynamically
        G4LogicalVolumeStore* store = G4LogicalVolumeStore::GetInstance();
        G4LogicalVolume* log_sphere = store->GetVolume("U235Sphere", false);
        G4double radius = 8.7; // Default in cm
        if (log_sphere) {
            G4VSolid* solid = log_sphere->GetSolid();
            G4Sphere* solid_sphere = dynamic_cast<G4Sphere*>(solid);
            if (solid_sphere) {
                radius = solid_sphere->GetOuterRadius() / cm;
            } else {
                G4cout << "Warning: U235Sphere solid is not a G4Sphere, using default radius." << G4endl;
            }
        } else {
            G4cout << "Warning: U235Sphere logical volume not found, using default radius." << G4endl;
        }
        G4double volume = (4.0 / 3.0) * M_PI * radius * radius * radius; // cm³

        // Retrieve histograms by name
        G4int birthId = man->GetH1Id("birth");
        G4int deathId = man->GetH1Id("death");
        G4int prodId = man->GetH1Id("prod");
        G4int lossId = man->GetH1Id("loss");
        G4int tauId = man->GetH1Id("tau");
        G4int fluxId = man->GetH1Id("flux_logE");
        G4int fissId = man->GetH3Id("fission_sites");
        G4int trackId = man->GetH3Id("track_length");

        auto* birth_h = man->GetH1(birthId);
        auto* death_h = man->GetH1(deathId);
        G4double cum_b = 0.0;
        G4double cum_d = 0.0;
        if (birth_h && death_h) {
            for (G4int i = 0; i < 300; ++i) {
                cum_b += birth_h->bin_height(i);
                cum_d += death_h->bin_height(i);
            }
        } else {
            G4cout << "Warning: Birth or death histogram not found." << G4endl;
        }

        auto* tau_h = man->GetH1(tauId);
        G4double tau = 0.0;
        G4int entries = 0;
        if (tau_h) {
            for (G4int i = 0; i < 100; ++i) {
                tau += tau_h->bin_height(i) * (i * 0.5);
                entries += static_cast<G4int>(tau_h->bin_height(i));
            }
            if (entries > 0) {
                tau /= entries;
            } else {
                G4cout << "Warning: No entries in tau histogram." << G4endl;
            }
        } else {
            G4cout << "Warning: Tau histogram not found." << G4endl;
        }

        auto* prod_h = man->GetH1(prodId);
        auto* loss_h = man->GetH1(lossId);
        G4double sum_p = 0.0;
        G4double sum_l = 0.0;
        if (prod_h && loss_h) {
            for (G4int i = 50; i <= 250; ++i) {
                sum_p += prod_h->bin_height(i);
                sum_l += loss_h->bin_height(i);
            }
        } else {
            G4cout << "Warning: Prod or loss histogram not found." << G4endl;
        }
        G4double k_p = (sum_l > 0) ? sum_p / sum_l : 0.0;

        G4double alpha = (cum_b > cum_d) ? (cum_d - cum_b) / (run->GetNumberOfEvent() * 300.0) : 0.0;
        G4double beta = -alpha * tau;
        G4double keff = (1 - beta != 0) ? k_p / (1 - beta) : 0.0;

        G4cout << "Run completed: " << run->GetNumberOfEvent() << " events" << G4endl;
        G4cout << "Neutron decay constant alpha: " << alpha / ns << " ns^-1" << G4endl;
        G4cout << "Neutron removal time tau: " << tau / ns << " ns" << G4endl;
        G4cout << "Prompt k_p: " << k_p << G4endl;
        G4cout << "Effective beta: " << beta << G4endl;
        G4cout << "Total effective multiplication factor keff: " << std::fixed << std::setprecision(5) << keff << G4endl;

        // Calculate Shannon entropy
        auto* fiss_h = man->GetH3(fissId);
        G4double total_fiss = 0.0;
        G4double entropy = 0.0;
        if (fiss_h) {
            for (G4int ix = 0; ix < 20; ++ix) {
                for (G4int iy = 0; iy < 20; ++iy) {
                    for (G4int iz = 0; iz < 20; ++iz) {
                        total_fiss += fiss_h->bin_height(ix, iy, iz);
                    }
                }
            }
            if (total_fiss > 0) {
                for (G4int ix = 0; ix < 20; ++ix) {
                    for (G4int iy = 0; iy < 20; ++iy) {
                        for (G4int iz = 0; iz < 20; ++iz) {
                            G4double count = fiss_h->bin_height(ix, iy, iz);
                            if (count > 0) {
                                G4double p = count / total_fiss;
                                entropy -= p * std::log2(p);
                            }
                        }
                    }
                }
            } else {
                G4cout << "Warning: No fission events recorded." << G4endl;
            }
        } else {
            G4cout << "Warning: Fission sites histogram not found." << G4endl;
        }
        G4cout << "Shannon entropy of fission sites: " << entropy << G4endl;

        // Flux spectrum output
        auto* flux_h = man->GetH1(fluxId);
        if (flux_h) {
            G4int num_events = run->GetNumberOfEvent();
            if (num_events == 0) {
                G4cout << "Warning: No events in the run." << G4endl;
                man->Write();
                man->CloseFile();
                return;
            }
            G4double min_logE = -3.0;
            G4double max_logE = 7.3;
            G4int nbins = 500;
            G4double dlogE = (max_logE - min_logE) / nbins;

            // Normalization for absolute flux assuming 1 W power
            G4double assembly_power = 1.0; // W
            G4double ev_per_joule = 1.0 / 1.60217662e-19;
            G4double power_ev_s = assembly_power * ev_per_joule;
            G4double energy_per_fission_ev = 200e6;
            G4double nu = 2.43; // Average neutrons per fission for U-235
            G4double fission_rate = power_ev_s / energy_per_fission_ev;
            G4double neutron_source_rate = fission_rate * nu;

            G4cout << "Assumed neutron source rate for 1W power: " << neutron_source_rate << " n/s" << G4endl;

            G4cout << "Energy flux spectrum (group fluence per source [n/cm^2 per source] vs E_low [eV]):" << G4endl;
            for (G4int i = 0; i < nbins; ++i) {
                G4double logE_low = min_logE + i * dlogE;
                G4double E_low = std::pow(10.0, logE_low);
                G4double bin_height = flux_h->bin_height(i);
                G4double group_fluence_per_source = (volume > 0.0) ? bin_height / volume / num_events : 0.0;
                G4cout << E_low << " " << group_fluence_per_source << G4endl;
            }

            G4cout << "Absolute energy flux spectrum for 1W power (group flux [n/cm^2/s] vs E_low [eV]):" << G4endl;
            for (G4int i = 0; i < nbins; ++i) {
                G4double logE_low = min_logE + i * dlogE;
                G4double E_low = std::pow(10.0, logE_low);
                G4double bin_height = flux_h->bin_height(i);
                G4double group_fluence_per_source = (volume > 0.0) ? bin_height / volume / num_events : 0.0;
                G4double group_flux_abs = group_fluence_per_source * neutron_source_rate;
                G4cout << E_low << " " << group_flux_abs << G4endl;
            }
        } else {
            G4cout << "Warning: Flux histogram not found." << G4endl;
        }
    }
    man->Write();
    man->CloseFile();
}
