function [imShift, predShift] = alignChannels_fst(im, maxShift)


global scale_ratio;
scale_ratio = 2;  

blockno=2;
predShift = zeros(2,2);
predShift0 = zeros(2,2);
imShift = im;

[maxv1, predShift1] = alignChannels_fst_subfunc(impyramid(im,'reduce'), ceil(maxShift/scale_ratio), predShift0,blockno); %(1:scale_ratio:siz(1),1:scale_ratio:siz(2),:)

tmp_max=zeros(1,2);
for i = 1:blockno
    [maxv, predShift2] = alignChannels_fst_subfunc(im, [scale_ratio+1, scale_ratio+1], predShift1(:,:,i).*scale_ratio,1);
    for j = 1:2
        if maxv(j,1) > tmp_max(j)
            tmp_max(j) = maxv(j,1);
            predShift(j,:) = predShift2(j,:,1);
        end
    end
end
for i=1:2
    imShift(:,:,i+1) = circshift(imShift(:,:,i+1), predShift(i,:));
end

% % ALIGNCHANNELS align channels in an image.
% %   [IMSHIFT, PREDSHIFT] = ALIGNCHANNELS(IM, MAXSHIFT) aligns the channels in an
% %   NxMx3 image IM. The first channel is fixed and the remaining channels
% %   are aligned to it within the maximum displacement range of MAXSHIFT (in
% %   both directions). The code returns the aligned image IMSHIFT after
% %   performing this alignment. The optimal shifts are returned as in
% %   PREDSHIFT a 2x2 array. PREDSHIFT(1,:) is the shifts  in I (the first) 
% %   and J (the second) dimension of the second channel, and PREDSHIFT(2,:)
% %   are the same for the third channel.
%  % Faster alignment
% global scale_ratio;
% scale_ratio = 3;  
% siz = size(im);
% predShift0 = zeros(2,2);
% [imtemp, predShift0] = alignChannels_subfunc(im(1:scale_ratio:siz(1),1:scale_ratio:siz(2),:), maxShift/scale_ratio, maxShift/scale_ratio, predShift0);
% [imShift, predShift] = alignChannels_subfunc(im, maxShift, [scale_ratio, scale_ratio], predShift0.*scale_ratio);
%          
% 
% %------------------------%
% 
% 
% % % Sanity check
% % assert(size(im,3) == 3);
% % assert(all(maxShift > 0));
% % 
% % % Dummy implementation (replace this with your own)
% % 
% % predShift = zeros(2, 2);
% % imShift = im;
% % 
% % siz_im = size(im);
% % im_edge = im;
% % for i=1:3
% %     im_edge(:,:,i) = (edge(im(:,:,i), 'sobel')+im(:,:,i)*0.5)/1.5;
% % end
% % im0 = im_edge(maxShift(1)+1:siz_im(1)-maxShift(1), maxShift(2)+1:siz_im(2)-maxShift(2),:);            
% % siz_tpl = size(im0);
% % 
% % mid1 = round((siz_im(1)+siz_tpl(1))/2);
% % mid2 = round((siz_im(2)+siz_tpl(2))/2);
% % for i=2:3
% %     c = normxcorr2(im0(:,:,i),im_edge(:,:,1));
% %     [a,b] = find(c==max(max(c(mid1-maxShift:mid1+maxShift,mid2-maxShift:mid2+maxShift))));
% %     a = a - mid1;
% %     b = b - mid2;
% %     predShift(i-1,:) = [a,b];
% %     imShift(:,:,i) = circshift(im(:,:,i), predShift(i-1,:));
% % end
% %         