%% Load data and set initial values
clear
clc
addpath(fullfile("..","functions"))
% load data containing time-dependent size and position(x, y)
nanoparticle_data = readmatrix('Monomeric_2.csv');
pixel_len = 0.53; % pixel length in nm
rho_s = 63.2825; % Pt monomer number density in atom/nm^-3
dt = 0.5; % Time interval in second
%% Convert Pixel data to radius and monomer number data
timedata = (1:size(nanoparticle_data, 1))'.*dt;
radidata = pixel_len.*sqrt(nanoparticle_data(:,1:3:end) ./ pi);
sigmas = (3 / (4 * pi * rho_s))^(1/3);
% Convert to monomer number: assuming the spherical geomerty
radidata(radidata == 0) = nan;
Preprocessing_in_situ(timedata, radidata,rho_s,dt,[1,40],"expdata_PtCOD.mat", ...
    window = 10)