
function [Nonelinear_Scalespace]=image_space(im,Scale_Invariance,nOctaves,scale_value,sigma_1)
%% Ĭ�ϲ�������
if nargin < 3
    nOctaves               = 3;          %  ������Ӱ�����������.  
end
if nargin < 4
    scale_value       = 1.6;       %  Ӱ��߶�����ϵ��   
end
if nargin < 5
    sigma_1           = 1.6;         %  ��һ���ͼ��ĳ߶ȣ�Ĭ����1.6.
end
%%  
% im  = double(im );
% im  = (im  - min(im ( : ) ) ) / ( max(im ( : ) ) - min(im ( : ) ) );
% im = single(im2uint8(im));
%% ��Ӱ��ת��Ϊ�Ҷ�ͼ
[~,~,num1]=size(im);
if(num1>=3)
    dst=rgb2gray(im(:,:,1:3));
else
    dst=im;
end
% ��Ӱ��ת��Ϊ������Ӱ����ֵ��[0~1]֮�� 
image=im2double(dst);
[M,N]=size(image);
%% �ж��Ƿ���Ҫ����������
if (strcmp(Scale_Invariance  ,'YES'))
    Layers=1;
else
    Layers=nOctaves;
end
%% ��ʼ��������Ӱ��cell�ռ�
Nonelinear_Scalespace=cell(1,Layers);
for i=1:1:Layers
    Nonelinear_Scalespace{i}=zeros(M,N);
end
%���ȶ�����ͼ����и�˹ƽ��

windows_size=2*round(2*sigma_1)+1;
W=fspecial('gaussian',[windows_size windows_size],sigma_1);      % Fspecial�������ڴ���Ԥ������˲�����
image=imfilter(image,W,'replicate');                                              %base_image�ĳ߶���sigma_1  % ����������������άͼ������˲���
Nonelinear_Scalespace{1}=image;                                                 %base_image��Ϊ�߶ȿռ�ĵ�һ��ͼ��

%% ���������Գ߶ȿռ�
for i=2:1:Layers
    %֮ǰ��ķ�������ɢ��ĵ�ͼ��,�����ݶ�֮ǰ����ƽ����Ŀ����Ϊ����������
    prev_image=Nonelinear_Scalespace{i-1};
    prev_image=imresize(prev_image,1/scale_value,'bilinear');
    Nonelinear_Scalespace{i}=(prev_image);
end
end

