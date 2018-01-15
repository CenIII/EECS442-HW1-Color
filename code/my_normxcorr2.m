function [gama] = my_normxcorr2(template, image, o_point, direction)

siz_im = size(image);
siz_tpl = size(template);

u=floor(siz_im(1)/2)+o_point(1)+direction(1);
v=floor(siz_im(2)/2)+o_point(2)+direction(2);

temp_mean = mean(mean(template));
image_mean = mean(mean(image));

template = template - temp_mean;
%image = image - image_mean;

f1 = floor(siz_tpl(1)/2);
f2 = floor(siz_tpl(2)/2);
a1 = u+1-f1;
a2 = u+f1+mod(siz_tpl(1),2);
b1 = v+1-f2;
b2 = v+f2+mod(siz_tpl(2),2);

if a1<1 || a2 > siz_im(1) || b1 < 1 || b2 > siz_im(2)
    gama = -10000;
else
    mat = image(a1:a2,b1:b2);
    mat = mat - sum(sum(mat))/(siz_tpl(1)*siz_tpl(2));
    gama = sum(sum(mat.*template))/((sum(sum(mat.^2))*sum(sum(template.^2)))^(1/2));
end