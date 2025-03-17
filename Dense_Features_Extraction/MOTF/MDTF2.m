function [MDTF_,dMDTG]=MDTF2(nonelinear_space,Multi_direction,angle)
%% ��ά������ɢ������ǿ����
%% ��ʼ����洢�ռ�
[M,N,P]=size(nonelinear_space);
MDTF_cell=cell(1,Multi_direction);
MDTG=zeros(M,N,P);
MDTF_=zeros(M,N,P);
for i=1:Multi_direction
    MDTF_cell{i}=zeros(M,N);
end
%% ����ͼ��Ĳ���ݶȺͲ�ֽ�
h=[-1,0,1;-2,0,2;-1,0,1];%����˲�ģ��
for o =1:Multi_direction
    %һ�ײ���ݶȺ�һ�ײ�ֽǶ�
    gradient_x_1=imfilter(nonelinear_space,h,'replicate');
    gradient_y_1=imfilter(nonelinear_space,h','replicate');
    %% ������ڱ�Ե��ǿ��ɢ�Ķ����ݶ�����
    Options=struct( 'eigenmode',2,'C', 1e-10, 'm',10,'alpha',0.01,'lambda_e',0.01,'lambda_h',0.3);
    [mu1,mu2,v1x,v1y,v2x,v2y]=EigenVectors2D(gradient_x_1.^2,gradient_x_1.*gradient_y_1,gradient_y_1.^2);
    gradA=gradient_x_1.^2+gradient_y_1.^2;      %% Gradient magnitude squared  �ݶȷ���ƽ��
    [Dxx,Dxy,~]=ConstructDiffusionTensor2D(mu1,mu2,v1x,v1y,v2x,v2y,gradA,Options);  % ��Ե��ǿ��ɢ EED 
    
    MDTF_cell{o}=abs(cos((o-1)*angle).*Dxx+sin((o-1)*angle).*Dxy);
    MDTG(:,:,o) =MDTF_cell{o};
%     figure, imshow(MDTG(:,:,o),[]);
end
sigma =0.8;
f = fspecial( 'gaussian', max( 1, fix( 6 * sigma + 1 ) ), sigma );
f = cat( 2, 1 * f, 2 * f, 1 * f );
dMDTG1 = convn(MDTG, f, 'same' );
sum1 = sum(dMDTG1, 3 );
dMDTG = dMDTG1 ./ ( sum1 + 0.000000001 );
end







































