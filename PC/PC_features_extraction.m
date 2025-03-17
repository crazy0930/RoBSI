function  [key_Point_end,gradient_1,angle_1]  =  PC_features_extraction(nonelinear_space,s,o,sigma_1,ratio,window,Scale_Invariance,A_patch,k)

%% 4 鲁棒相位一致性梯度图 
% tic;
[Bolb_space_1,Corner_space_1,gradient_1,angle_1]=PC_gradient_feature(nonelinear_space,Scale_Invariance,s,o);% 4 相位一致性梯度图 
% disp(['构建相位一致性梯度图:  ',num2str(toc),'S']);
% figure,imshow(Bolb_space_1{1},[]);
% figure,imshow(Corner_space_1{1},[]);
% figure,imshow(gradient_1{1},[]);
% figure,imshow(angle_1{1},[]);
%% 5  feature point extraction
points_layer = 1200;       % 默认值是1500
position_1=AF_Initial_Feature_Points_extreme2(Bolb_space_1,Corner_space_1,sigma_1,ratio,points_layer,gradient_1,angle_1,Scale_Invariance,A_patch);
% position_1=AF_Initial_Feature_Points_extreme(Bolb_space_1,Corner_space_1,sigma_1,ratio,points_layer,gradient_1,angle_1,Scale_Invariance);

%% 聚合特征优化策略筛选
% 1、非极大值抑制
keypoints = AF_selectMax(position_1, window);                                  % 非极大值抑制 
% 2、边界点消除
% nOctaves =3;
% if (strcmp(Scale_Invariance  ,'YES'))
%     Layers=1;
% else
%     Layers=nOctaves;
% end
% KeyPts_cell=[];
% for i=1:1:Layers
%     temp=Bolb_space_1{i};%每层的尺度
%     KeyPts= AF_Boundary_Points_Filtering(temp,keypoints.kpts, A_patch); 
%     KeyPts_cell=[KeyPts_cell;KeyPts];
% end
% 3、重复点滤除
KeyPts_cell=keypoints.kpts;
uni1=KeyPts_cell(:,[1,2]);
[~,i,~]=unique(uni1,'rows','first');
KeyPts_cell=KeyPts_cell(sort(i)',:);
%  4、根据得分值排序策略
key_Point_end=sortrows(KeyPts_cell,6,'descend');
KeyNum=round(k*size(KeyPts_cell,1));
if(size(key_Point_end,1)>KeyNum)
    key_Point_end=key_Point_end(1:KeyNum,:);
else
    key_Point_end=key_Point_end(:,:);
end
end