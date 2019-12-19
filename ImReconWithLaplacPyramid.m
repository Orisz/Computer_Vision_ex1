function [reconstructed_image] = ImReconWithLaplacPyramid(pyramid_decomposition)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fucnction name:    ImReconWithLaplacPyramid
% Fucnction input:   'pyramid_decomposition' - this is the images laplacian
% pyramid decomposition.
% Fucnction output:  'reconstructed_image' - the image reconstruction using
%                       the laplacian pyrmid
% Fucnction description: this function takes in images laplacian
% pyramid decomposition and attends to reconstruct the original image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%as shown in class, when the levels are of the same size we need just to
%sum all the levels to get the images reconstruction
reconstructed_image = sum(pyramid_decomposition, 3);

end

