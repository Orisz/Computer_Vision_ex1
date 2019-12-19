function [gain_map] = CalcGain(in_image_energy_pyramid,example_image_energy_pyramid)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fucnction name:    CalcGain
% Fucnction input:   'in_image_energy_pyramid' - energy pyramid of the in
% image
%                    'example_image_energy_pyramid' - enrgy oyramid of the
%                    example image
% Fucnction output:  'gain_map' - pyramids of the gains
% Fucnction description: calcs the gain pyramid as dedcribed in the pdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eps = 1e-4;
min = 0.9;
max = 2.8;

gain_map = sqrt(example_image_energy_pyramid ./(in_image_energy_pyramid + eps));
gain_map(gain_map < min) = min;
gain_map(gain_map > max) = max;
end

