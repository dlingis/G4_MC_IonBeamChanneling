# G4_MC_IonBeamChanneling
Geant4 channeling package based on channeling package developed by Enrico Bagli. 
This is just my solution to the problem I had with the incorrect energy loss evaluation with the newest channeling 
package (which comes with geant10.5). The modification is based on not biasing the cross section of the reaction, 
but by wrapping the physics process based on the Electron and Nuclei density at the current step of the particle, 
as was done in older channeling package (where only planar channeling was considered). The code works, I have been 
able to recreate the experimental results of protons (see doi:10.1103/physrev.161.330). I haven't been able to compare 
alpha particle channeling, and current interest is heavy ion (P, B) channeling in silicon. 

The code is based on Hadr06 example, merged with G4_MC_Channeling (by Enrico Bagli, https://github.com/ebagli/G4_MC_CHANNELING) and, 
once again, merged with older channeling package. There are also some additions by me, merged from my older source files.  

I do not hold any rights, as it is a open source code. However, if you use the GEANT4 code, or channeling package, don't forget 
to cite the developers of the code. 


