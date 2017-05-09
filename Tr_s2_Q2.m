function [ output_args ] = Tr_s2_Q2(dr,dt,simstruct)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%use a structure to pass all simulation parameters
ss = simstruct;
E=ss.E;Q=ss.Q;nu=ss.nu;Gamma=ss.Gamma;eta=ss.eta
%returns the compoenents in front of v1, v2, and v3 basis for G_ijkl
%generalized response function, including factors of Gamma, eta, etc
[I1 I2 I3]=Tr_s2_factors(dr,dt,E,nu,Gamma,eta);

%now need contractions of Q with various v basis.
%v1 = x_i x_j x_k x_l / x^4, v2 = delta_ij x_k x_l / x^2 + permutations
%v3 = delta_ij delta_kl + permutations

end

