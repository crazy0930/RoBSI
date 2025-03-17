function [dCSOG] = sobel_9dir(img,n,simag)
% ����Sobel����
Gx = [-1 0 1; -2 0 2; -1 0 1];
Gy = Gx';
% ����ˮƽ����ֱ�ݶ�
Ix = conv2(double(img), Gx, 'same');
Iy = conv2(double(img), Gy, 'same');
% �����ݶȷ���
Gdir = zeros(size(img, 1), size(img, 2), n);
for i = 1:n
    theta = (i-1)*(360/n);
    Gdir(:,:,i) = abs((Ix.*cosd(theta) + Iy.*sind(theta)));
end
f = fspecial( 'gaussian', max( 1, fix( 6 * simag + 1 ) ), simag );
f = cat( 3, 1 * f, 3 * f, 1 * f );
dCSOG1 = convn( Gdir, f, 'same' );
dCSOG = dCSOG1( :, :, 2:end  - 1 );
sum1 = sum( dCSOG, 3 );
dCSOG = dCSOG ./ ( sum1 + 0.000000001 );
end