function [k,Nx] = image2permeability(filepath,k_min,k_max,Ny)

%% - Image to read in

img = imread(filepath);

aspect_ratio = length(img(:,1,1))/length(img(1,:,1));

%% - Inputs
Nx = round(Ny/aspect_ratio);

%% - Resize image
if size(img,3) == 3
    img = rgb2gray(img);
end

img = im2double(img);

% Resize to match the simulation grid
img = imresize(img,[Ny Nx]);

% Display image
figure
imshow(img)
title('Input image')

%% - Permeability range
% Map grayscale to log10 permeability
logk = log10(k_min) + img * (log10(k_max) - log10(k_min));

% Convert back to permeability
k = 10.^logk;

%% - Display permeability field

figure
imagesc(log10(k))
axis image
colorbar
title('log_{10}(k)')

%% - Reshape
k = reshape(k',Nx*Ny,1);

end
