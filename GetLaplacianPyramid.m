function [pyramid_decomposition] = GetLaplacianPyramid(input_image,num_of_levels)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fucnction name:    GetLaplacianPyramid
% Fucnction input:   'input_image' - this is the image we nees
%                     to decompose at grey level
%                     'num_of_levels' - the desired level of the pyramid
% Fucnction output:  'pyramid_decomposition' - the image decomposed to its
%                     laplacian pyramid
% Fucnction description: this function takes in image from the user and
% decompose it to its laplacian pyramid of level as given from the user
% (see 2nd argument). please note that all levels are from the same size
% (as the input image size) as required in the instructions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[rows , cols] = size(input_image);%get image size

%init gauss & laplac pyramids
gauss_pyramid =  zeros(rows, cols, num_of_levels);
gauss_pyramid(:, :, 1) = input_image;%at gauss the first level is the original image
lap_pyrmaid = zeros(rows, cols, num_of_levels);

%build the pyramids
for i = 2:num_of_levels
   curr_sigma =  2 ^ i;
   prev_gauss_level = gauss_pyramid(:, :, i - 1);
   %please note that 'imgaussfilt' makes sure that the filt size is atleast
   %5 times bigger then 'curr_sigma'
   curr_gauss_level = imgaussfilt(prev_gauss_level, curr_sigma);
   gauss_pyramid(:, :, i) = curr_gauss_level;
   lap_pyrmaid(:, :, i - 1) = prev_gauss_level - curr_gauss_level;
end
%take care of the last level of the lap pyramid
lap_pyrmaid(:, :, num_of_levels) = gauss_pyramid(:, :, num_of_levels);

pyramid_decomposition = lap_pyrmaid;
end

