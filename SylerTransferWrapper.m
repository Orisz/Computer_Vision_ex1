function [] = SylerTransferWrapper(input_im, input_bg, input_mask, example_im, levels)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fucnction name:    SylerTransferWrapper
% Fucnction input:   'input_im' - input image we want to transfer syle to
%                    'input_bg' - the bg for the input image
%                    'input_mask' - the mask for the input image
%                    'example_im' - the image we want it 'style'
%                    'levels' - number of levels for the pyramids
% Fucnction output:  ' ' - non
% Fucnction description: wrapper functin to avoid code duplication
%                       transfers the style from example image to input
%                       image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%laplacian pyramids for input image and example image. we calc for each
%image 3 pyramids, one for each chanel:
%input image:
in_image_pyramid_R      = GetLaplacianPyramid(input_im(:,:,1),levels);
in_image_pyramid_G      = GetLaplacianPyramid(input_im(:,:,2),levels);
in_image_pyramid_B      = GetLaplacianPyramid(input_im(:,:,3),levels);
%example image:
example_image_pyramid_R = GetLaplacianPyramid(example_im(:,:,1),levels);
example_image_pyramid_G = GetLaplacianPyramid(example_im(:,:,2),levels);
example_image_pyramid_B = GetLaplacianPyramid(example_im(:,:,3),levels);

%calculate image energy for each channel
%input image:
in_image_energy_pyramid_R      = CalcEnergy(in_image_pyramid_R);
in_image_energy_pyramid_G      = CalcEnergy(in_image_pyramid_G);
in_image_energy_pyramid_B      = CalcEnergy(in_image_pyramid_B);
%example image:
example_image_energy_pyramid_R = CalcEnergy(example_image_pyramid_R);
example_image_energy_pyramid_G = CalcEnergy(example_image_pyramid_G);
example_image_energy_pyramid_B = CalcEnergy(example_image_pyramid_B);

%calculate gain map for each lap level
gain_map_pyramid_R = CalcGain(in_image_energy_pyramid_R,example_image_energy_pyramid_R);
gain_map_pyramid_G = CalcGain(in_image_energy_pyramid_G,example_image_energy_pyramid_G);
gain_map_pyramid_B = CalcGain(in_image_energy_pyramid_B,example_image_energy_pyramid_B);

% 'section e'
%condtruct output image pyramid:
%R
output_image_pyramid_R = gain_map_pyramid_R .* in_image_pyramid_R;
output_image_pyramid_R(:, :, levels) = example_image_pyramid_R(:, :, levels);
%G
output_image_pyramid_G = gain_map_pyramid_G .* in_image_pyramid_G;
output_image_pyramid_G(:, :, levels) = example_image_pyramid_G(:, :, levels);
%B
output_image_pyramid_B = gain_map_pyramid_B .* in_image_pyramid_B;
output_image_pyramid_B(:, :, levels) = example_image_pyramid_B(:, :, levels);

% 'section f'
%reconstruct each chanel of the out put image:
output_image_R = ImReconWithLaplacPyramid(output_image_pyramid_R);
output_image_G = ImReconWithLaplacPyramid(output_image_pyramid_G);
output_image_B = ImReconWithLaplacPyramid(output_image_pyramid_B);

output_image = zeros(size(input_im));
output_image(:, :, 1) = output_image_R;
output_image(:, :, 2) = output_image_G;
output_image(:, :, 3) = output_image_B;

output_image_new_bg = ChangeImBg(output_image, input_mask, input_bg);

figure();
subplot(1,3,1);
imshow(input_im);
title('Input Image');
subplot(1,3,2);
imshow(example_im);
title('Example Image');
subplot(1,3,3);
imshow(output_image_new_bg);
title('Output Image');

end

