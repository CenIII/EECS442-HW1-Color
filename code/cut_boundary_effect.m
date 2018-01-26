function [im_cut] = cut_boundary_effect(im)

im_cut=rgb2gray(im);
[height width] = size(im_cut);


top = im_cut(1:floor(0.1*height),:);
bottom = im_cut(end-floor(0.1*height)+1:end,:);
left = im_cut(:,1:floor(0.1*width));
right = im_cut(:,(end-floor(0.1*width)+1):end);

top_cut = find_cut_line(top,1);
bottom_cut = find_cut_line(bottom,2);
left_cut = find_cut_line(left,3);
right_cut = find_cut_line(right,4);

im_cut = im(top_cut:end-bottom_cut+1, left_cut:end-right_cut+1, :);





function [cut_depth] = find_cut_line(bd,idx)

bd = edge(bd,'canny');
siz=size(bd);

Threshold = 0.6;
switch (idx) 
    case 1
        match_line = ones(1,siz(2));
        maxValue = 0;
        tmp=zeros(1,siz(1));
        for i=1:siz(1)
            tmp(i)=sum(bd(i,:).*match_line);
            if tmp(i) > maxValue
                maxValue = tmp(i);
                %cut_depth = i;
            end
        end
        cut_depth = max(find(tmp>(maxValue*Threshold)));
    case 2
        match_line = ones(1,siz(2));
        maxValue = 0;
        tmp=zeros(1,siz(1));
        for i=1:siz(1)
            tmp(i)=sum(bd(end-i+1,:).*match_line);
            if tmp(i) > maxValue
                maxValue = tmp(i);
                %cut_depth = i;
            end
        end
        cut_depth = max(find(tmp>(maxValue*Threshold)));
    case 3
        match_line = ones(1,siz(1));
        maxValue = 0;
        tmp=zeros(1,siz(2));
        for i=1:siz(2)
            tmp(i)=sum(bd(:,i).*match_line');
            if tmp(i) > maxValue
                maxValue = tmp(i);
                %cut_depth = i;
            end
        end
        cut_depth = max(find(tmp>(maxValue*Threshold)));
    case 4
        match_line = ones(1,siz(1));
        maxValue = 0;
        tmp=zeros(1,siz(2));
        for i=1:siz(2)
            tmp(i)=sum(bd(:,end-i+1).*match_line');
            if tmp(i) > maxValue
                maxValue = tmp(i);
                %cut_depth = i;
            end
        end
        cut_depth = max(find(tmp>(maxValue*Threshold)));

end

        
