function [maxv, predShiftsub] = alignChannels_subfunc(im, maxShift, opoint, block_num)

% Sanity check
assert(size(im,3) == 3);
assert(all(maxShift > 0));

% Dummy implementation (replace this with your own)
if block_num == 1
    predShiftsub = zeros(2, 2);
else
    predShiftsub = zeros(2, 2, block_num);
end
maxv = zeros(2,block_num);
%imShift = im;
for i=1:3
    im(:,:,i) = (edge(im(:,:,i),'canny'))+(edge(im(:,:,i),'prewitt'))*0.6+(edge(im(:,:,i),'log'))*0.6;
end

maxValue = zeros(2*maxShift(1)+1,2*maxShift(2)+1,2);
for layer = 2:3
    for i = opoint(layer-1, 1)-maxShift(1):opoint(layer-1, 1)+maxShift(1)
        for j = opoint(layer-1, 2)-maxShift(2):opoint(layer-1, 2)+maxShift(2)
            if i<-15 || j<-15 || i>15 || j>15
                continue;
            end
            tmp_im = circshift(im(:,:,layer), [i,j]);
            dot_prod = dot(reshape(im(:,:,1).',1,[]),reshape(tmp_im.',1,[]));
            maxValue(i-opoint(layer-1, 1)+maxShift(1)+1, j-opoint(layer-1, 2)+maxShift(2)+1,layer-1) = dot_prod;      
        end
    end
    
    for k = 1:block_num
        tmp = maxValue(:,:,layer-1);
        [a,b] = find(tmp==max(max(tmp)));
        if(block_num==1)
            predShiftsub(layer-1,:)=[a+opoint(layer-1, 1)-maxShift(1)-1,b+opoint(layer-1, 2)-maxShift(2)-1];
        else
            predShiftsub(layer-1,:,k)=reshape([a+opoint(layer-1, 1)-maxShift(1)-1,b+opoint(layer-1, 2)-maxShift(2)-1],1,2,1);
        end
        maxv(layer-1, k) = maxValue(a,b,layer-1);
        maxValue(a,b,layer-1) = -10000;
    end
    %imShift(:,:,layer) = circshift(imShift(:,:,layer), predShift(layer-1,:));
end



% 
% % Sanity check
% assert(size(im,3) == 3);
% assert(all(maxShift > 0));
% 
% % Dummy implementation (replace this with your own)
% 
% predShift = zeros(2, 2);
% imShift = im;
% 
% siz_im = size(im);
% im_edge = im;
% for i=1:2
%     im_edge(:,:,i) = (edge(im(:,:,i), 'sobel'))*1;%+im(:,:,i)*0.5;
% end
% im_edge(:,:,3) = (edge(im(:,:,3), 'sobel'))*1;
% 
% im0 = im_edge;%(cut_width(1)+1:siz_im(1)-cut_width(1), cut_width(2)+1:siz_im(2)-cut_width(2),:);            
% siz_tpl = size(im0);
% 
% cor_mid1 = floor((siz_im(1)+siz_tpl(1))/2);
% cor_mid2 = floor((siz_im(2)+siz_tpl(2))/2);
% 
% for i=2:3
%     maxVal = -10000;
%     %c=normxcorr2(im0(:,:,i), im_edge(:,:,1));
%     for k=-maxShift(1):maxShift(1)
%         for l=-maxShift(2):maxShift(2)
%             ncc = my_normxcorr2(im0(:,:,i),im_edge(:,:,1),o_point(i-1,:),[k,l]);
%             %c1 = c(cor_mid1+o_point(i-1,1)+k,cor_mid2+o_point(i-1,2)+l);
%             if ncc > maxVal
%                 maxVal = ncc;
%                 predShift(i-1,:) = o_point(i-1,:)+[k,l];
%             end
%         end
%     end
%     imShift(:,:,i) = circshift(im(:,:,i), predShift(i-1,:));
% end
        