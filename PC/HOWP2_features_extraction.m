function  [key_Point_end,gradient_1,angle_1]  =  HOWP2_features_extraction(nonelinear_space,layers,s,o,sigma_1,ratio,window,Scale_Invariance,A_patch,k)

%% 4 ³����λһ�����ݶ�ͼ 
% tic;
[Bolb_space_1,Corner_space_1,~,~]=PC_gradient_feature(nonelinear_space,Scale_Invariance,s,o,layers);
% disp(['������λһ�����ݶ�ͼ:  ',num2str(toc),'S']);

%%  ����Ӱ���ݶ�����
[gradient_1,angle_1]  =Gradient2(nonelinear_space,layers);

%% 5  feature point extraction
points_layer = 2000;       % Ĭ��ֵ��1500
position_1=AF_Initial_Feature_Points_extreme3(Bolb_space_1,Corner_space_1,sigma_1,ratio,points_layer,gradient_1,angle_1,Scale_Invariance,A_patch,layers);
%% �ۺ������Ż�����ɸѡ
% 1���Ǽ���ֵ����
keypoints = AF_selectMax(position_1, window);                                  % �Ǽ���ֵ���� 
% 3���ظ����˳�
KeyPts_cell=keypoints.kpts;
uni1=KeyPts_cell(:,[1,2]);
[~,i,~]=unique(uni1,'rows','first');
KeyPts_cell=KeyPts_cell(sort(i)',:);
%  4�����ݵ÷�ֵ�������
key_Point_end=sortrows(KeyPts_cell,6,'descend');
KeyNum=round(k*size(KeyPts_cell,1));
if(size(key_Point_end,1)>KeyNum)
    key_Point_end=key_Point_end(1:KeyNum,:);
else
    key_Point_end=key_Point_end(:,:);
end
end