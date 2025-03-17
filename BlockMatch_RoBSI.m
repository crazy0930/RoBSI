clear;close all;
warning('off');
addpath(genpath('display'));
addpath(genpath('HOWP_Func'));
addpath(genpath('Dense_Features_Extraction'));
addpath(genpath('MOTF_Func'));
addpath(genpath('FSC'));
addpath(genpath('PC'));
addpath(genpath('ANMS'));
addpath(genpath('OSS'));
%% 1 Import and display reference and image to be registered
file_image = 'testdata';
[filename,pathname]=uigetfile({'*.*','All Files(*.*)'},'Select Image',file_image);
image_1=imread(strcat(pathname,filename));
if size(image_1)==3
    image_1 = rgb2gray(image_1);
end
[filename,pathname]=uigetfile({'*.*','All Files(*.*)'},'Select Image',file_image);
image_2=imread(strcat(pathname,filename));
if size(image_2)==3
    image_2 = rgb2gray(image_2);
end
t1 = clock;
%% 2 
% key parameters
image_block_size = 512;      % 匹配小块的尺寸
image_match_detect = 30;     % 影像匹配的忽略边缘？
block_num_x = 4;             % 横向分区数
block_num_y = 4;             % 纵向分区数 
block_size_x = floor(size(image_1,2)/block_num_x);      % 横向区宽
block_size_y = floor(size(image_1,1)/block_num_y);      % 纵向区宽
image_block_stride_x = 512;         % 区内横向取影像小块步长
image_block_stride_y = 512;         % 区内纵向取影像小块步长
image_block_num_x = floor((block_size_x-image_block_size)/image_block_stride_x)+1;    % 区内横向小块数
image_block_num_y = floor((block_size_y-image_block_size)/image_block_stride_y)+1;    % 区内纵向小块数
block_min_matchpoints_num = 100;  % 区成功匹配点数
block_min_points_num = 80;
match_block_points_1 = [];
match_block_points_2 = [];
match_points_1 = [];
match_points_2 = [];

%%  MOTIF初始参数
disthre = 15;                          %   he threshod of match errors the deflaut is 5. for
Detectors = 'B_Harris';         %   选择特征检测器方法: 'GBMS', 'ANMS'， 'B_Harris',  'Harris', 'FAST', 'KAZE'
Descriptors = 'MOTF';          %    (1)'MFOG';  (2)'MSFG';  (3)'CSOG'; （4）‘MOTF’   (5)'MSOG';  (6) 'MOGG'   (7) 'SOND' %%  对比算法  'CFOG'    'FHOG'   'DLSS'    'FSURF'  'HOPC'    
kpts_nums = 100;                %   特征点数目;
Tw=96;                              %   模板窗口大小;
Sw=5;                                  %   搜索窗口
T =[1,0,0;0,1,0;0,0,1];

for i = 1:block_num_y
    for j = 1:block_num_x
        is_match = 0;      %判断当前区是否已经匹配成功
        [i,j]
        numbers = randperm(image_block_num_y*image_block_num_x); % 生成1到n的随机序列
        block_id = 1;
        while(is_match == 0)
            if block_id == image_block_num_y*image_block_num_x
                is_match = 1;
            end
            m = floor(numbers(block_id)/image_block_num_x)+1;
            n = mod(numbers(block_id), image_block_num_x);
            if n==0
                n = image_block_num_x;
                m = m-1;
            end
            init_x_index = (j-1)*block_size_x+(n-1)*image_block_stride_x+1;
            init_y_index = (i-1)*block_size_y+(m-1)*image_block_stride_y+1;
            image_1_block = image_1(init_y_index:init_y_index+image_block_size,init_x_index:init_x_index+image_block_size);
            image_2_block = image_2(init_y_index:init_y_index+image_block_size,init_x_index:init_x_index+image_block_size);
            %%  template matching using MOTF
            if all(image_1_block(:)==0)||all(image_2_block(:)==0)
                block_match_points_1 = [];
                block_match_points_2 = [];
            else
                [block_match_points_1,block_match_points_2,CMR,MAPE] = MOTF_match(image_1_block,image_2_block,T,Detectors,Descriptors, disthre,kpts_nums,Tw,Sw);
            end
            block_id = block_id+1;
            % [block_match_points_1, block_match_points_2] = AHM_match(image_1_block, image_2_block);
            if size(block_match_points_1,1)>block_min_points_num
                block_match_points_1(:,1) = block_match_points_1(:,1)+init_x_index;
                block_match_points_1(:,2) = block_match_points_1(:,2)+init_y_index;
                block_match_points_2(:,1) = block_match_points_2(:,1)+init_x_index;
                block_match_points_2(:,2) = block_match_points_2(:,2)+init_y_index;
                match_block_points_1 = [match_block_points_1;block_match_points_1];
                match_block_points_2 = [match_block_points_2;block_match_points_2];
                if size(match_block_points_1,1)>block_min_matchpoints_num
                    match_points_1 = [match_points_1;match_block_points_1];
                    match_points_2 = [match_points_2;match_block_points_2];
                    is_match = 1;
                    match_block_points_1 = [];
                    match_block_points_2 = [];
                end
            end
        end
    end
end
[corrRefPt,corrSenPt] = ErrorDect(double(match_points_1),double(match_points_2),0,1.0);
t2 = clock;
cp_showMatch(image_1, image_2, corrRefPt(:,[1,2]), corrSenPt(:,[1,2]),[],'After Subpixel Fineing');
num2str(etime(t2,t1))

im_Ref = image_1;
im_Sen = image_2;
%融合
Siz1 = size(im_Ref);
Disp.Disp = 1;
Disp.DCel = fix( max( Siz1 ) / 15 );
Disp.CropP = 0;

Model = { 'Projective', 'Polynomial', 'LPiecewise', 'MQuadrics', 'PointWise' };
Opt.GeoMdl = Model{ 5 };
Opt.n = 2;
Opt.d = 0;
Opt.k = 5;
Opt.pw = 2;
if size(im_Ref, 3 ) > 1
    I1 = rgb2gray(im_Ref( :, :, 1:3 ) );
else
    I1 =im_Ref;
end
if size(im_Sen, 3 ) > 1
    I2 = rgb2gray(im_Sen( :, :, 1:3 ) );
else
    I2 = im_Sen;
end

%% 计算融合图
[G,R] = ImgWarp(I1, I2, corrRefPt, corrSenPt, Opt, Disp);
% figure
% imshow(G,[]);
% figure
% imshow(R,[]);

grid_num=15;
grid_size=floor(min(size(G,1),size(R,2))/grid_num);
[~,~,f_3]=mosaic_map(G,R,grid_size);
figure
imshow(f_3,[]);









