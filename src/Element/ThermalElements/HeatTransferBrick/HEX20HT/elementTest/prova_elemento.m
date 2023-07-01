close all
clearvars
clc


%% Deifne Material
% material properties
E = 200e9;              % 70e9 % 200e9 % Young's modulus
rho = 7850;             % 2700 % 7850 % density
k = 1;                  % conductivity
alpha = 10;              % heat convective coefficient
c = 1;                  % heat capacity
nu = 0.3;               % nu
kappa = 1e8;            % material damping modulus 1e8
Emiss = 0;              % Emissivety constant
alpha_T = 11.7e-6;      % thermal expansion coefficient
loadfactor = 50;        % mechanical load factor 


myMaterial  = Material();
set(myMaterial,'THERMAL_CONDUCTIVITY', k, 'SPECIFIC_HEAT', c, 'DENSITY', rho, 'CONVECTION_COEFFICIENT', alpha, 'EMISSIVITY', Emiss);

%% Create element

Ngauss = 2;

HTElement = HEX8ElementHT(myMaterial,Ngauss);

HTElement.nodes = 0.2*[-1 -1 -1
                 +1 -1 -1
                 +1 +1 -1
                 -1 +1 -1
                 -1 -1 +1
                 +1 -1 +1
                 +1 +1 +1
                 -1 +1 +1];
             
HTElement.nodeIDs = 1:1:8;


%% conductivity matrix

K = conductivity_matrix(HTElement);

%% Thermal capacity matrix

M = thermal_capacity_matrix(HTElement);

%% Neumann Element

NeumEl = HEX8NeuElementHT(myMaterial,Ngauss);

NeumEl.nodes = 0.2*[-1 -1 -1
                 -1 +1 -1
                 -1 +1 +1
                 -1 -1 +1];
             
NeumEl.nodeIDs = 1:1:4;


q = randn(8,1);
NeumEl.thermal_flux_matrix()











