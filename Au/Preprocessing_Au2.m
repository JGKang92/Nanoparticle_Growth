%% Load data and set initial values
clear
clc
addpath(fullfile("..","functions"))
% load data containing time-dependent size and position(x, y)
monomeric_data = load('Au_monomeric_pixel.csv');
dissolution_data = load('Au_dissolution_pixel.csv');
nanoparticle_data = [monomeric_data(:,1:3:end), dissolution_data(:,1:3:end)];
pixel_len = 0.538; % pixel length in nm
rho_s = 59.0088; % Au monomer number density in atom/nm^-3
dt = 0.5; % Time interval in second
%% Convert Pixel data to radius and monomer number data
sigmas = (3/(4*pi*rho_s))^(1/3);
lent = size(nanoparticle_data, 1);
timedata = (1:lent)'.*dt;
radidata = pixel_len.*sqrt(nanoparticle_data ./ pi);
Preprocessing_in_situ(timedata, radidata, rho_s, dt, [30, 30], "expdata_Au.mat", ...
    trim = 200);