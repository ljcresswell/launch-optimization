clear
clc
close all
addpath(genpath(pwd));
params = readParamsFile('UserInputValues.txt');
N_nodes = 10;
t_max = 10000;
[u0] = generateInitialGuess(params, N_nodes, t_max);