function output = demosaicImage(im, method)
% DEMOSAICIMAGE computes the color image from mosaiced input
%   OUTPUT = DEMOSAICIMAGE(IM, METHOD) computes a demosaiced OUTPUT from
%   the input IM. The choice of the interpolation METHOD can be 
%   'baseline', 'nn', 'linear', 'adagrad'. 

switch lower(method)
    case 'baseline'
        output = demosaicBaseline(im);
    case 'nn'
        output = demosaicNN(im);         % Implement this
    case 'linear'
        output = demosaicLinear(im);     % Implement this
    case 'adagrad'
        output = demosaicAdagrad(im);    % Implement this
end

%--------------------------------------------------------------------------
%                          Baseline demosaicing algorithm. 
%                          The algorithm replaces missing values with the
%                          mean of each color channel.
%--------------------------------------------------------------------------
function mosim = demosaicBaseline(im)
mosim = repmat(im, [1 1 3]); % Create an image by stacking the input
[imageHeight, imageWidth] = size(im);

% Red channel (odd rows and columns);
redValues = im(1:2:imageHeight, 1:2:imageWidth);
meanValue = mean(mean(redValues));
mosim(:,:,1) = meanValue;
mosim(1:2:imageHeight, 1:2:imageWidth,1) = im(1:2:imageHeight, 1:2:imageWidth);

% Blue channel (even rows and colums);
blueValues = im(2:2:imageHeight, 2:2:imageWidth);
meanValue = mean(mean(blueValues));
mosim(:,:,3) = meanValue;
mosim(2:2:imageHeight, 2:2:imageWidth,3) = im(2:2:imageHeight, 2:2:imageWidth);

% Green channel (remaining places)
% We will first create a mask for the green pixels (+1 green, -1 not green)
mask = ones(imageHeight, imageWidth);
mask(1:2:imageHeight, 1:2:imageWidth) = -1;
mask(2:2:imageHeight, 2:2:imageWidth) = -1;
greenValues = mosim(mask > 0);
meanValue = mean(greenValues);
% For the green pixels we copy the value
greenChannel = im;
greenChannel(mask < 0) = meanValue;
mosim(:,:,2) = greenChannel;

%--------------------------------------------------------------------------
%                           Nearest neighbour algorithm
%--------------------------------------------------------------------------
function mosim = demosaicNN(im)
mosim = demosaicBaseline(im);
%
% Implement this 
%
[imageHeight, imageWidth,~] = size(mosim);

H_even = imageHeight - mod(imageHeight,2);
W_even = imageWidth - mod(imageWidth,2);

% green 2 
    % red covered from right
mosim(1:2:imageHeight, 1:2:W_even, 2) = mosim(1:2:imageHeight, 2:2:W_even, 2);
    % last col of red covered from left
mosim(1:2:imageHeight, imageWidth, 2) = mosim(1:2:imageHeight, W_even, 2);
    % blue covered from left
mosim(2:2:H_even, 2:2:W_even, 2) = mosim(2:2:H_even, 1:2:W_even, 2);

% R B / G
% g_temp = mosim(:,:,2);
% g_temp = linear_interp_zero(g_temp);
% mosim(:,:,1) = mosim(:,:,1) ./ g_temp;
% mosim(:,:,3) = mosim(:,:,3) ./ g_temp;

% red
mosim(:,2:2:W_even,1) = mosim(:,1:2:W_even,1);
mosim(2:2:H_even,:,1)=mosim(1:2:H_even,:,1);

% blue
mosim(:,1:2:W_even,3)=mosim(:,2:2:W_even,3);
mosim(1:2:H_even,:,3)=mosim(2:2:H_even,:,3);
mosim(:,imageWidth,3)=mosim(:,W_even,3);
mosim(imageHeight,:,3)=mosim(H_even,:,3);
% 
% mosim(:,:,1) = mosim(:,:,1) .* g_temp;
% mosim(:,:,3) = mosim(:,:,3) .* g_temp;

%--------------------------------------------------------------------------
%                           Linear interpolation
%--------------------------------------------------------------------------
function mosim = demosaicLinear(im)
[imageHeight, imageWidth,~] = size(im);

mosim(1:2:imageHeight, 1:2:imageWidth,1) = im(1:2:imageHeight, 1:2:imageWidth);
mosim(2:2:imageHeight, 2:2:imageWidth,3) = im(2:2:imageHeight, 2:2:imageWidth);
% We will first create a mask for the green pixels (+1 green, -1 not green)
mask = ones(imageHeight, imageWidth);
mask(1:2:imageHeight, 1:2:imageWidth) = -1;
mask(2:2:imageHeight, 2:2:imageWidth) = -1;
% For the green pixels we copy the value
greenChannel = im;
greenChannel(mask < 0) = 0;
mosim(:,:,2) = greenChannel;
%
% Implement this 
%

H_even = imageHeight - mod(imageHeight,2);
W_even = imageWidth - mod(imageWidth,2);

temp = zeros(imageHeight, imageWidth);
temp1 = zeros(imageHeight, imageWidth);
temp2 = zeros(imageHeight, imageWidth);
ones_mask = ones(imageHeight, imageWidth);
mask = zeros(imageHeight, imageWidth);

%green
%?
temp(:,1:imageWidth-1) = temp(:,1:imageWidth-1)+mosim(:,2:imageWidth,2);
mask(:,1:imageWidth-1) = mask(:,1:imageWidth-1)+ones_mask(:,2:imageWidth);

%?
temp(:,2:imageWidth) = temp(:,2:imageWidth)+mosim(:,1:imageWidth-1,2);
mask(:,2:imageWidth) = mask(:,2:imageWidth)+ones_mask(:,1:imageWidth-1);

%?
temp(1:imageHeight-1,:) = temp(1:imageHeight-1,:)+mosim(2:imageHeight,:,2);
mask(1:imageHeight-1,:) = mask(1:imageHeight-1,:)+ones_mask(2:imageHeight,:);

%?
temp(2:imageHeight,:) = temp(2:imageHeight,:)+mosim(1:imageHeight-1,:,2);
mask(2:imageHeight,:) = mask(2:imageHeight,:)+ones_mask(1:imageHeight-1,:);

%divide by 4
mosim(:,:,2) = temp./mask + mosim(:,:,2); %TODO: not 4  ???1 matrix ??division mask


%R  B / Green
g_temp = mosim(:,:,2);
g_temp = linear_interp_zero(g_temp);
mosim(:,:,1) = mosim(:,:,1) ./ g_temp;
mosim(:,:,3) = mosim(:,:,3) ./ g_temp;

%red
ones_mask_red = zeros(imageHeight, imageWidth);
ones_mask_red(1:2:imageHeight,1:2:imageWidth) = 1;
mask2 = zeros(imageHeight, imageWidth);
mask2(:,1:imageWidth-1) = mask2(:,1:imageWidth-1)+ones_mask_red(:,2:imageWidth);
mask2(:,2:imageWidth) = mask2(:,2:imageWidth)+ones_mask_red(:,1:imageWidth-1);
mask2(1:imageHeight-1,:) = mask2(1:imageHeight-1,:)+ones_mask_red(2:imageHeight,:);
mask2(2:imageHeight,:) = mask2(2:imageHeight,:)+ones_mask_red(1:imageHeight-1,:);
mask2(mask2==0) = 1;

%?
temp1(:,1:imageWidth-1) = temp1(:,1:imageWidth-1)+mosim(:,2:imageWidth,1);
%?
temp1(:,2:imageWidth) = temp1(:,2:imageWidth)+mosim(:,1:imageWidth-1,1);
%?
temp1(1:imageHeight-1,:) = temp1(1:imageHeight-1,:)+mosim(2:imageHeight,:,1);
%?
temp1(2:imageHeight,:) = temp1(2:imageHeight,:)+mosim(1:imageHeight-1,:,1);

temp1 = temp1./mask2; %TODO: not 2

%?
temp2(:,1:imageWidth-1) = temp2(:,1:imageWidth-1)+temp1(:,2:imageWidth);
%?
temp2(:,2:imageWidth) = temp2(:,2:imageWidth)+temp1(:,1:imageWidth-1);
%?
temp2(1:imageHeight-1,:) = temp2(1:imageHeight-1,:)+temp1(2:imageHeight,:);
%?
temp2(2:imageHeight,:) = temp2(2:imageHeight,:)+temp1(1:imageHeight-1,:);

temp2 = temp2./mask;  %TODO: not 4

mosim(:,:,1) = mosim(:,:,1) + temp1;
mosim(2:2:imageHeight,2:2:imageWidth,1) = temp2(2:2:imageHeight,2:2:imageWidth);

%blue
temp1 = zeros(imageHeight, imageWidth);
temp2 = zeros(imageHeight, imageWidth);

ones_mask_blue = zeros(imageHeight, imageWidth);
ones_mask_blue(2:2:imageHeight,2:2:imageWidth) = 1;
mask3 = zeros(imageHeight, imageWidth);
mask3(:,1:imageWidth-1) = mask3(:,1:imageWidth-1)+ones_mask_blue(:,2:imageWidth);
mask3(:,2:imageWidth) = mask3(:,2:imageWidth)+ones_mask_blue(:,1:imageWidth-1);
mask3(1:imageHeight-1,:) = mask3(1:imageHeight-1,:)+ones_mask_blue(2:imageHeight,:);
mask3(2:imageHeight,:) = mask3(2:imageHeight,:)+ones_mask_blue(1:imageHeight-1,:);
mask3(mask3==0) = 1;

%?
temp1(:,1:imageWidth-1) = temp1(:,1:imageWidth-1)+mosim(:,2:imageWidth,3);
%?
temp1(:,2:imageWidth) = temp1(:,2:imageWidth)+mosim(:,1:imageWidth-1,3);
%?
temp1(1:imageHeight-1,:) = temp1(1:imageHeight-1,:)+mosim(2:imageHeight,:,3);
%?
temp1(2:imageHeight,:) = temp1(2:imageHeight,:)+mosim(1:imageHeight-1,:,3);

temp1 = temp1./mask3; %TODO: not 2

%?
temp2(:,1:imageWidth-1) = temp2(:,1:imageWidth-1)+temp1(:,2:imageWidth);
%?
temp2(:,2:imageWidth) = temp2(:,2:imageWidth)+temp1(:,1:imageWidth-1);
%?
temp2(1:imageHeight-1,:) = temp2(1:imageHeight-1,:)+temp1(2:imageHeight,:);
%?
temp2(2:imageHeight,:) = temp2(2:imageHeight,:)+temp1(1:imageHeight-1,:);

temp2 = temp2./mask; %TODO: not 4

mosim(:,:,3) = mosim(:,:,3) + temp1;
mosim(1:2:imageHeight,1:2:imageWidth,3) = temp2(1:2:imageHeight,1:2:imageWidth);

mosim(:,:,1) = mosim(:,:,1) .* mosim(:,:,2);
mosim(:,:,3) = mosim(:,:,3) .* mosim(:,:,2);


%--------------------------------------------------------------------------
%                           Adaptive gradient
%--------------------------------------------------------------------------
function mosim = demosaicAdagrad(im)
%
% Implement this 
%
[imageHeight, imageWidth,~] = size(im);

mosim(1:2:imageHeight, 1:2:imageWidth,1) = im(1:2:imageHeight, 1:2:imageWidth);
mosim(2:2:imageHeight, 2:2:imageWidth,3) = im(2:2:imageHeight, 2:2:imageWidth);
% We will first create a mask for the green pixels (+1 green, -1 not green)
mask = ones(imageHeight, imageWidth);
mask(1:2:imageHeight, 1:2:imageWidth) = -1;
mask(2:2:imageHeight, 2:2:imageWidth) = -1;
% For the green pixels we copy the value
greenChannel = im;
greenChannel(mask < 0) = 0;
mosim(:,:,2) = greenChannel;



% Green ?????  ???   ???
for i=1:imageHeight
    for j=1:imageWidth
        if(mod(i+j,2)==0)
            if(i>1&&j>1&&i<imageHeight&&j<imageWidth)
                vert = mosim(i+1,j,2) - mosim(i-1,j,2);
                hori = mosim(i,j+1,2) - mosim(i,j-1,2);
                if(abs(vert)>abs(hori))
                    mosim(i,j,2) = mosim(i,j-1,2) + hori/2;
                else
                    mosim(i,j,2) = mosim(i-1,j,2) + vert/2;
                end
            elseif((i==1 || i==imageHeight) && (j~=1&&j~=imageWidth))
                hori = mosim(i,j+1,2) - mosim(i,j-1,2);
                mosim(i,j,2) = mosim(i,j-1,2) + hori/2;
            elseif((j==1 || j==imageWidth) && (i~=1&&i~=imageHeight))
                vert = mosim(i+1,j,2) - mosim(i-1,j,2);
                mosim(i,j,2) = mosim(i-1,j,2) + vert/2;
            else
                ii=abs(i-2)+1;
                jj=abs(j-2)+1;
                mosim(i,j,2) = (mosim(i,jj,2)+mosim(ii,j,2))/2;
            end
        end
    end
end


% R B / Green 
g_temp = mosim(:,:,2);
g_temp = linear_interp_zero(g_temp);
mosim(:,:,1) = mosim(:,:,1) ./ g_temp;
mosim(:,:,3) = mosim(:,:,3) ./ g_temp;

% Red  ????????adaptive
ones_mask_red = zeros(imageHeight, imageWidth);
ones_mask_red(1:2:imageHeight,1:2:imageWidth) = 1;
mask2 = zeros(imageHeight, imageWidth);
mask2(:,1:imageWidth-1) = mask2(:,1:imageWidth-1)+ones_mask_red(:,2:imageWidth);
mask2(:,2:imageWidth) = mask2(:,2:imageWidth)+ones_mask_red(:,1:imageWidth-1);
mask2(1:imageHeight-1,:) = mask2(1:imageHeight-1,:)+ones_mask_red(2:imageHeight,:);
mask2(2:imageHeight,:) = mask2(2:imageHeight,:)+ones_mask_red(1:imageHeight-1,:);
mask2(mask2==0) = 1;

temp1 = zeros(imageHeight, imageWidth);

%?
temp1(:,1:imageWidth-1) = temp1(:,1:imageWidth-1)+mosim(:,2:imageWidth,1);
%?
temp1(:,2:imageWidth) = temp1(:,2:imageWidth)+mosim(:,1:imageWidth-1,1);
%?
temp1(1:imageHeight-1,:) = temp1(1:imageHeight-1,:)+mosim(2:imageHeight,:,1);
%?
temp1(2:imageHeight,:) = temp1(2:imageHeight,:)+mosim(1:imageHeight-1,:,1);

temp1 = temp1./mask2; %TODO: not 2
mosim(:,:,1) = mosim(:,:,1) + temp1;

for i=2:2:imageHeight
    for j=2:2:imageWidth
        if(i>1&&j>1&&i<imageHeight&&j<imageWidth)
            vert = mosim(i+1,j,1) - mosim(i-1,j,1);
            hori = mosim(i,j+1,1) - mosim(i,j-1,1);
            if(abs(vert)>abs(hori))
                mosim(i,j,1) = mosim(i,j-1,1) + hori/2;
            else
                mosim(i,j,1) = mosim(i-1,j,1) + vert/2;
            end
        elseif((i==1 || i==imageHeight) && (j~=1&&j~=imageWidth))
            hori = mosim(i,j+1,1) - mosim(i,j-1,1);
            mosim(i,j,1) = mosim(i,j-1,1) + hori/2;
        elseif((j==1 || j==imageWidth) && (i~=1&&i~=imageHeight))
            vert = mosim(i+1,j,1) - mosim(i-1,j,1);
            mosim(i,j,1) = mosim(i-1,j,1) + vert/2;
        else
            ii=abs(i-2)+1;
            jj=abs(j-2)+1;
            mosim(i,j,1) = (mosim(i,jj,1)+mosim(ii,j,1))/2;
        end
    end
end
% Blue same with Red
temp1 = zeros(imageHeight, imageWidth);

ones_mask_blue = zeros(imageHeight, imageWidth);
ones_mask_blue(2:2:imageHeight,2:2:imageWidth) = 1;
mask3 = zeros(imageHeight, imageWidth);
mask3(:,1:imageWidth-1) = mask3(:,1:imageWidth-1)+ones_mask_blue(:,2:imageWidth);
mask3(:,2:imageWidth) = mask3(:,2:imageWidth)+ones_mask_blue(:,1:imageWidth-1);
mask3(1:imageHeight-1,:) = mask3(1:imageHeight-1,:)+ones_mask_blue(2:imageHeight,:);
mask3(2:imageHeight,:) = mask3(2:imageHeight,:)+ones_mask_blue(1:imageHeight-1,:);
mask3(mask3==0) = 1;

%?
temp1(:,1:imageWidth-1) = temp1(:,1:imageWidth-1)+mosim(:,2:imageWidth,3);
%?
temp1(:,2:imageWidth) = temp1(:,2:imageWidth)+mosim(:,1:imageWidth-1,3);
%?
temp1(1:imageHeight-1,:) = temp1(1:imageHeight-1,:)+mosim(2:imageHeight,:,3);
%?
temp1(2:imageHeight,:) = temp1(2:imageHeight,:)+mosim(1:imageHeight-1,:,3);

temp1 = temp1./mask3; %TODO: not 2

mosim(:,:,3) = mosim(:,:,3) + temp1;

for i=1:2:imageHeight
    for j=1:2:imageWidth
        if(i>1&&j>1&&i<imageHeight&&j<imageWidth)
            vert = mosim(i+1,j,3) - mosim(i-1,j,3);
            hori = mosim(i,j+1,3) - mosim(i,j-1,3);
            if(abs(vert)>abs(hori))
                mosim(i,j,3) = mosim(i,j-1,3) + hori/2;
            else
                mosim(i,j,3) = mosim(i-1,j,3) + vert/2;
            end
        elseif((i==1 || i==imageHeight) && (j~=1&&j~=imageWidth))
            hori = mosim(i,j+1,3) - mosim(i,j-1,3);
            mosim(i,j,3) = mosim(i,j-1,3) + hori/2;
        elseif((j==1 || j==imageWidth) && (i~=1&&i~=imageHeight))
            vert = mosim(i+1,j,3) - mosim(i-1,j,3);
            mosim(i,j,3) = mosim(i-1,j,3) + vert/2;
        else
            ii=abs(i-2)+1;
            jj=abs(j-2)+1;
            mosim(i,j,3) = (mosim(i,jj,3)+mosim(ii,j,3))/2;
        end
    end
end

mosim(:,:,1) = mosim(:,:,1) .* mosim(:,:,2);
mosim(:,:,3) = mosim(:,:,3) .* mosim(:,:,2);








function g_temp = linear_interp_zero(g_temp)
[imageHeight,imageWidth] = size(g_temp);
for i=1:imageHeight
    for j=1:imageWidth
        if g_temp(i,j)==0
            a=0;b=0;c=0;d=0;count=0;
            if i-1>0 
                a=g_temp(i-1,j);
                count=count+1;
            end
            if i+1<=imageHeight
                b=g_temp(i+1,j);
                count=count+1;
            end
            if j-1>0
                c=g_temp(i,j-1);
                count=count+1;
            end
            if j+1<=imageWidth
                d=g_temp(i,j+1);
                count=count+1;
            end
            g_temp(i,j)=(a+b+c+d)/count;
        end
    end
end

