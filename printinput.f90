subroutine printinput


print*, 'Non-equilibrium BD simulator'
print*, 'NEQBD v3.11'
print*, 'MT 2013/2014'
print*, ' Scope: '
print*, ' This program is a Brownian Dynamics simulator for particles interacting via oscillatory potentials '
print*, ' It is designed to simulate systems of two different types of particles'
print*, 
print*, 
print*, '********************* Format for input.txt **************************"
print*,
print*,

print*, '# dimensions'
print*, 'Number of dimensions, program well tested only for 2 and some features for 3'
print*, '# Number of particles'
print*, ' Number of particles'
print*, '# Type of potential'
print*, 'Options:'
print*, 'LJ0: Repulsive LJ12, e*(D12/r)^12, no parameters'
print*, 'LJN: LJ with even coefficients: e12*(D12/r)^A - e*alfa*(D12/r)^B, params: A, B, alfa'
print*, 'YUK: attractive/repulsive Yukawa between A and B + LJ:  e*(D12/r)^12 + z1*z2*eexp*exp(-r/elen)/r, params: eexp, elen'
print*, 'ATR: Attractive Yukawa between all particles'
print*, 'TAB: Tabulated'
print*, '# Parameters'
print*, '[parameters here]'
print*, '#jump = (D*dt)^0.5'
print*, '[jump in BD dynamics in sigma units, equal to (D*dt)^0.5]'
print*, '#temperature'
print*, '[Temperature, use 1.0 in general]'
print*, '# vdW diameters: D1ini, D1fin, D2ini, D2fin'
print*, '[vdW diameters in sigma units]'
print*, '# vdW energies: e1ini, e1fin, e2ini, e2fin'
print*, '[energies in vdW in kBT units, note that LJ uses "e" and not "4*e" as prefactor]'
print*, '# readinput, 1 = read MCin.xyz and use last coordinates'
print*, 'Options:'
print*, '0 : start from a square distribution'
print*, '1 : read from file BDin.xyz'
print*, '2 : read from file BDin.xyz and reanalize'
print*, '3 : read from file and perform one cycle, ignore number of cycles below'
print*, '-1: start from a square with A and B arranged in alternated order'  
print*, '# switchtype'
print*, 'Type of switching, options:'
print*, '0 = step jump'
print*, '1 = sine'
print*, '2 = linear in radius'
print*, '3 = linear in area'
print*, '4 = linear in volume'
print*, '5 = linear in pH'
print*, '6,7 = non-linear in pH, see paper'
print*, 'other value = no jump, D,e and z are fixed according to the initial step and oscillation type 1'
print*, '# interparticle distance'
print*, 'Average interparticle distance, density = 1/(i.d)^n'
print*, '#fa'
print*, 'Fraction of A-type particles'
print*, 'Switching half period in BD steps'
print*, 'Number of steps in a half-period'
print*, '# number of BD steps'
print*, 'Number of BD steps to simulate'
print*, 'If readinput = 1 is used, the simulation starts in the last step of BDin.xyz'
print*, 'If readinput = 3 is used, only 1 cycle is simulated'
print*, '# save coordinates every'
print*, 'Number of steps to save in BD.xyz'
print*, '# save coordinate offset'
print*, 'Offset to save in BD.xyz, typically is period/2'
print*, '# save GRs'
print*, 'Calculate and save radial distrubution functions, 0 = no, 1 = yes'
print*, '# save BO'
print*, 'Calculate and save bond-order parameters, 0 = no, 1 = yes'
print*, '# detect crystal'
print*, 'Calculate and save bond-order parameters, 0 = no, 1 = yes'
print*, '# save dissipation'
print*, 'Calculate dissipation, 0=no, 1=yes'
print*, '# save every period'
print*, 'saves every this number of oscillation period'
print*, '# points per period'
print*, '20
# cutoff
8.0
# calculate pressure every
20000
# Barostat?, constant, target pressure
0 1000.0 10.0
