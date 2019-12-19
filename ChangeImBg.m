function [im_with_new_bg] = ChangeImBg(in_image_rgb, in_mask, in_bg)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fucnction name:   ChangeImBg
% Fucnction input:  'in_image_rgb' - input image in rgb
%                   'in_mask' - binary mask '0' at the backgraound
%                   'in_bg' - the new bg for the image
% Fucnction output:  'im_with_new_bg' - the image with the new bg
% Fucnction description: takes in an image and replace its bg using the
% users mask and put instead the users bg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Im_no_bg = in_image_rgb .* in_mask;
not_mask = 1 - in_mask;
adjusted_bg = in_bg .* not_mask;
fin_image = Im_no_bg + adjusted_bg;
im_with_new_bg = fin_image;
end

