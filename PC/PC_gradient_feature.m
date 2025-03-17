function [Bolb_space,Corner_space,gradient_cell,angle_cell] =PC_gradient_feature(Scalespace, Scale_Invariance,n,o,nOctaves)
%implemented by Yongxiang Yao
%% 默认参数设置
if nargin < 5
    nOctaves           = 3;          %  构建的影像金字塔层数.  
end

%% 判断是否构建影像金字塔
if (strcmp(Scale_Invariance  ,'YES'))
    Layers=1;
else
    Layers=nOctaves;
end
%% cell初始化
[M,N]=size(Scalespace{1});
gradient_cell=cell(1,Layers);
angle_cell=cell(1,Layers);
Bolb_space=cell(1,Layers);
Corner_space=cell(1,Layers);
for j=1:Layers
    gradient_cell{j}=zeros(M,N);
    angle_cell{j}=zeros(M,N);
    Bolb_space{j}=zeros(M,N);
    Corner_space{j}=zeros(M,N);
end

 %% 通过相位一致性计算获取频率域内的影像梯度和主方向
 for i=1:Layers
    int_image=im2uint8(Scalespace{i});
%     int_image=Scalespace{i};
    [Max_image,Min_image,pc1,or1] = phasecong_PCGF(int_image,n,o);
    or1 = or1/pi*180;
%     or1(or1<0) = or1(or1<0)+360;
    Bolb_space{i}=Max_image;
    Corner_space{i}=Min_image;
    gradient_cell{i}=single(pc1);
    angle_cell{i}=single(or1);
%         figure, imshow(angle_cell{i},[]);
 end
