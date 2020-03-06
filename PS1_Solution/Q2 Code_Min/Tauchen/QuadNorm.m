% Quadnorm.m
% This function computes the discrete approximation to a standard AR(1)
% Normal process using Tauchen and Hussey's (1991) quadrature procedure
% as implemented by Burnside.
% Required inputs are the dimension of the state-space (nz), the mean (zbar), 
% standard deviation (std) and the degree of autocorrelation (rho)

function [Qz, z] = Quadnorm(nz, zbar, std, rho);

ghquad;

% wmat is the vector of weights
z = ghq(0.5*(nz - 1)*nz + 1 : 0.5*nz*(nz + 1), 2);
wmat = ghq(0.5*(nz - 1)*nz + 1 : 0.5*nz*(nz + 1), 3);

% Now transform the abscissas given the law of motion
z = zbar + z*std;

% generate the transition matrix
Qz = Trans(z, wmat, zbar, rho, std)';

clear wmat ghq