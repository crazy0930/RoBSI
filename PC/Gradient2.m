function [gradient_cell,angle_cell] =Gradient2(im,layers)
%implemented by Yongxiang Yao
%% Image format conversion
P=layers;
%% cell≥ı ºªØ
gradient_cell=cell(1,layers);
angle_cell=cell(1,layers);

%% Define filter operator
h1=[-1,0,1;-2,0,2;-1,0,1];            % Sobel_x 

 %% Image gradient calculation suppressed by LBP filter
for j=1:1:P
    %First order Sobel filter
    gradient_x_11=imfilter(im{j},h1,'replicate');
    gradient_y_11=imfilter(im{j},h1','replicate');
    gradient_22=sqrt(gradient_x_11.^2+gradient_y_11.^2);
    %Second order Sobel filter
    gradient_x_22=imfilter(gradient_22,h1,'replicate');
    gradient_y_22=imfilter(gradient_22,h1','replicate');
    gradient_33=sqrt(gradient_x_22.^2+gradient_y_22.^2);
    angle_22=atan2(gradient_y_22,gradient_x_22);
    angle_22=angle_22*180/pi;%-180~180
    angle_22(angle_22<0)=angle_22(angle_22<0)+360;% 0-360
    gradient_cell{j}=gradient_33;
    angle_cell{j}=angle_22;
end
