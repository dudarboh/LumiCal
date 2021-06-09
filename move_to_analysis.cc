//////////////////////////////////////////////////////
        //Simulation of efficiency of the calorimeter
        G4double S0_cal = 0.819;
        G4double p1_cal = 2.166;
        G4double p0 = 0.999 / 2.;
            // If hit is in calorimeter - reject hit with a probability based on efficiency curve
            // measured for the paper
            if(G4UniformRand() > (1. + std::erf((energy_in_mips - S0_cal) / p1_cal)) * p0) continue;
