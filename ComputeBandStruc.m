function [] = ComputeBandStruc()
%COMPUTEBANDSTRUC Compute 2-d band structure of optical lattice 
%   Using the potential output by the general lattice for N beams code, I
%   will use this program to compute the electronic band structure of the
%   system. Note, this needs to be in the same directory (or path) as the
%   general lattice complex function given that the latter needs to be
%   called to generate the potential (specifically the underlying
%   sinusoidal representation)

%% generate the lattice potential
%using the 6 cosine components of the lattice potential generated as output
%from Peter's code.
A = [1,1,0.6,0.5];
ph_deg = [0, 0, 90, -70];
th_deg = [0,90,180,270];
pol_deg = [0,0,0,0];
plots = 0; %boolean to turn plots on or off
[waveAmplitudes,deltaKsUnitless,deltaPhis]=GeneralLatticeComplexRepWithComponents(A,ph_deg,th_deg,pol_deg,plots);



end

