function [laplacian_pyramid_energy] = CalcEnergy(laplacian_pyramid)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fucnction name:    calcEnergy
% Fucnction input:   'laplacian_pyramid' - laplacian pyramid of an image
% Fucnction output:  'laplacian_pyramid_energy' - pyramids of energy. each
% level holds the corresponding local energy for the laplacian pyramid
% Fucnction description: for each level in the lap pyramid this func calcs
% the local energy as described in the pdf. finaly returns energy pyramid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[rows, cols, levels] = size(laplacian_pyramid);
energy_pyramid = zeros(rows, cols, levels);

for i=1:levels
    curr_sigma = 2^(i+1);
    POW_2_lap_level = laplacian_pyramid(:,:,i).^2;
    curr_local_energy = imgaussfilt(POW_2_lap_level, curr_sigma);
    energy_pyramid(:,:,i) = curr_local_energy;
end
laplacian_pyramid_energy = energy_pyramid;
end

