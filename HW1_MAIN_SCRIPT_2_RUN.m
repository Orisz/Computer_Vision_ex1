%% this is the main script for task #1 : Hough Transform.
% each section will excute a different section of this exrecise.
% submiters:
%   Ori Sztyglic 
%   Yossi Magriso 

%% Transform
% a.

close all

f=zeros(101,101);
f(1,1)=1;
[H, theta, rho] = hough(f);

figure;
imagesc(theta, rho, H);
colmap = colormap('gray');
colmap = colmap(end:-1:1,:);
colormap(colmap);
colorbar;
xlabel('\theta');
ylabel('\rho');
title('Hough Transform')

f(101,1)=1;
[H, theta, rho] = hough(f);

h = figure;
imagesc(theta, rho, H);
colmap = colormap('gray');
colmap = colmap(end:-1:1,:);
colormap(colmap);
colorbar;
xlabel('\theta');
ylabel('\rho');
title('Hough Transform')

f(1,101)=1;
[H, theta, rho] = hough(f);

figure;
imagesc(theta, rho, H);
colmap = colormap('gray');
colmap = colmap(end:-1:1,:);
colormap(colmap);
colorbar;
xlabel('\theta');
ylabel('\rho');
title('Hough Transform')

f(101,101)=1;
[H, theta, rho] = hough(f);

figure;
imagesc(theta, rho, H);
colmap = colormap('gray');
colmap = colmap(end:-1:1,:);
colormap(colmap);
colorbar;
xlabel('\theta');
ylabel('\rho');
title('Hough Transform')

% Not used:
f(51,51)=1;

%% b.

% Creating a square:
sqr_mask = zeros(101,101);
sqr_size = 40;

sqr_mask((51- sqr_size/2):(51 + sqr_size/2),(51- sqr_size/2):(51 + sqr_size/2)) = 1;
sqr_mask = imrotate(sqr_mask,30);

figure(100);
imagesc(sqr_mask);
colmap = colormap('gray');
colmap = colmap(end:-1:1,:);
colormap(colmap);
colorbar;

% c.

% Detecting the square's edges:
sqr_edges = edge(sqr_mask, 'Canny');

figure;
imagesc(sqr_edges);
colmap = colormap('gray');
colmap = colmap(end:-1:1,:);
colormap(colmap);
colorbar;

% d.

Theta_resolution  = 5; % [deg]
[Hough_of_sqr_edges, theta, rho] = hough(sqr_edges,'Theta',-90:Theta_resolution:(90-1));

figure;
imagesc(theta, rho, Hough_of_sqr_edges);
colmap = colormap('gray');
colmap = colmap(end:-1:1,:);
colormap(colmap);
colorbar;
xlabel('\theta');
ylabel('\rho');

%% Voting
% f.

Max_num_of_peaks = 4;
peaks = houghpeaks(Hough_of_sqr_edges,Max_num_of_peaks);

figure;
imagesc(theta, rho, Hough_of_sqr_edges);
colmap = colormap('gray');
colmap = colmap(end:-1:1,:);
colormap(colmap);
colorbar;
xlabel('\theta');
ylabel('\rho');
hold on;
plot(theta(peaks(:,2)),rho(peaks(:,1)),'o','color','red','MarkerSize',20);

%% Line linking:
% g.

lines = houghlines(sqr_edges,theta,rho,peaks,'MinLength',35);

figure(100);
hold on;

for k = 1:length(lines)
    xy = [lines(k).point1; lines(k).point2];
    plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
end

%% h.

% Finding lines in a picture of a building:

building_image = imread('building.jpg');

figure(101);
imshow(building_image);

building_edges = edge(building_image,'Canny',0.1,sqrt(2));

figure;
imshow(building_edges);
colmap = colormap('gray');
colmap = colmap(end:-1:1,:);
colormap(colmap);

% Show 9 different Canny parameters edge detection:
fig_ind = 1;
for sigma_factor = -1:1:1
    for Th_factor = -1:1:1
        Threshold = 0.1 * 2^Th_factor;
        Sigma = sqrt(2) * 2^sigma_factor;
        building_edges = edge(building_image,'Canny',Threshold,Sigma);
        figure;
        imshow(building_edges);
        colmap = colormap('gray');
        colmap = colmap(end:-1:1,:);
        colormap(colmap);
        title(['Threshold = ',num2str(Threshold),' , \sigma = ',num2str(Sigma)]);
        fig_ind = fig_ind+1;
    end
end


Theta_resolution = 1;
[Hough_of_building_edges, theta, rho] = hough(building_edges,'Theta',-90:Theta_resolution:(90-1));

Max_num_of_peaks = 100;
peaks = houghpeaks(Hough_of_building_edges,Max_num_of_peaks,'Threshold',70,'NHoodSize',[31,31]);

figure;
imagesc(theta, rho, Hough_of_building_edges);
colmap = colormap('gray');
colmap = colmap(end:-1:1,:);
colormap(colmap);
colorbar;
xlabel('\theta');
ylabel('\rho');
hold on;
plot(theta(peaks(:,2)),rho(peaks(:,1)),'o','color','red','MarkerSize',20);
title('Hough Transform')

lines = houghlines(building_edges,theta,rho,peaks,'MinLength',30);

figure(101);
hold on;

for k = 1:length(lines)
    xy = [lines(k).point1; lines(k).point2];
    plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
end


%% this is the main script for task #2 : laplacian Pyramids.
% each section will excute a different section of this exrecise.
% submiters:
%   Ori Sztyglic 201110830
%   Yossi Magriso 035847839
%% Section a
clc;clear all;close all;
% see attached fuction "GetLaplacianPyramid"
%% section b
% see attached function "ImReconWithLaplacPyramid"
%test
levels = 6;
Im_rgb = im2double(imread('data/Inputs/imgs/0004_6.png'));
Im = rgb2gray(Im_rgb);
pyramid_decomp = GetLaplacianPyramid(Im,levels);
Im_recon = ImReconWithLaplacPyramid(pyramid_decomp);
err = immse(Im, Im_recon);
fprintf('"sanaty check" => reconstruction error = %d\n',err);
figure(1);
subplot(1,2,1);
imshow(Im);
title('original');
subplot(1,2,2);
imshow(Im_recon);
title('reconstructed');
%% section c 

%import image
in_image_rgb = im2double(imread('data/Inputs/imgs/0004_6.png'));

%import mask
in_mask = im2double(imread('data/Inputs/masks/0004_6.png'));

%import bg
in_bg = im2double(imread('data/Examples/bgs/6.jpg'));

in_image_rgb_new_bg = ChangeImBg(in_image_rgb, in_mask, in_bg);

%plot results
figure(2);
subplot(1,2,1);
imshow(in_image_rgb);
title('Input Image, Original Background');
subplot(1,2,2);
imshow(in_image_rgb_new_bg);
title('Input Image New Background');

%% section d
example_image_rgb = im2double(imread('data/Examples/imgs/6.png'));

levels = 6;

%laplacian pyramids for input image and example image. we calc for each
%image 3 pyramids, one for each chanel:
%input image:
in_image_pyramid_R      = GetLaplacianPyramid(in_image_rgb(:,:,1),levels);
in_image_pyramid_G      = GetLaplacianPyramid(in_image_rgb(:,:,2),levels);
in_image_pyramid_B      = GetLaplacianPyramid(in_image_rgb(:,:,3),levels);
%example image:
example_image_pyramid_R = GetLaplacianPyramid(example_image_rgb(:,:,1),levels);
example_image_pyramid_G = GetLaplacianPyramid(example_image_rgb(:,:,2),levels);
example_image_pyramid_B = GetLaplacianPyramid(example_image_rgb(:,:,3),levels);

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

%% section e
%construct output image pyramid:
%R
output_image_pyramid_R = gain_map_pyramid_R .* in_image_pyramid_R;
output_image_pyramid_R(:, :, levels) = example_image_pyramid_R(:, :, levels);
%G
output_image_pyramid_G = gain_map_pyramid_G .* in_image_pyramid_G;
output_image_pyramid_G(:, :, levels) = example_image_pyramid_G(:, :, levels);
%B
output_image_pyramid_B = gain_map_pyramid_B .* in_image_pyramid_B;
output_image_pyramid_B(:, :, levels) = example_image_pyramid_B(:, :, levels);

%% section f
%reconstruct each chanel of the out put image:
output_image_R = ImReconWithLaplacPyramid(output_image_pyramid_R);
output_image_G = ImReconWithLaplacPyramid(output_image_pyramid_G);
output_image_B = ImReconWithLaplacPyramid(output_image_pyramid_B);

output_image = zeros(size(in_image_rgb));
output_image(:, :, 1) = output_image_R;
output_image(:, :, 2) = output_image_G;
output_image(:, :, 3) = output_image_B;

output_image_new_bg = ChangeImBg(output_image, in_mask, in_bg);

figure(3);
subplot(1,3,1);
imshow(in_image_rgb);
title('Input Image');
subplot(1,3,2);
imshow(example_image_rgb);
title('Example Image');
subplot(1,3,3);
imshow(output_image_new_bg);
title('Output Image');

%% section g

%additional styles for image 0004_6:
example_im_16 = im2double(imread('data/Examples/imgs/16.png'));
example_im_21 = im2double(imread('data/Examples/imgs/21.png'));
in_bg_16 = im2double(imread('data/Examples/bgs/16.jpg'));
in_bg_21 = im2double(imread('data/Examples/bgs/21.jpg'));
% the second input image & mask :
in_image_rgb_dude = im2double(imread('data/Inputs/imgs/0006_001.png'));
in_image_mask_dude = im2double(imread('data/Inputs/masks/0006_001.png'));

%the second in image examples images:
example_im_0 = im2double(imread('data/Examples/imgs/0.png'));
example_im_9 = im2double(imread('data/Examples/imgs/9.png'));
example_im_10 = im2double(imread('data/Examples/imgs/10.png'));
in_bg_0 = im2double(imread('data/Examples/bgs/0.jpg'));
in_bg_9 = im2double(imread('data/Examples/bgs/9.jpg'));
in_bg_10 = im2double(imread('data/Examples/bgs/10.jpg'));
%preform style transfer as described in the pdf:
SylerTransferWrapper(in_image_rgb, in_bg_16, in_mask, example_im_16, levels);
SylerTransferWrapper(in_image_rgb, in_bg_21, in_mask, example_im_21, levels);

%dude image:
SylerTransferWrapper(in_image_rgb_dude, in_bg_0, in_image_mask_dude, example_im_0, levels);
SylerTransferWrapper(in_image_rgb_dude, in_bg_9, in_image_mask_dude, example_im_9, levels);
SylerTransferWrapper(in_image_rgb_dude, in_bg_10, in_image_mask_dude, example_im_10, levels);


%% section h 


in_image_ori = im2double(imread('ori2.png'));
in_image_ori = imresize(in_image_ori,[1320 1000]);
SylerTransferWrapper_NO_MASK(in_image_ori, example_im_9, levels);


