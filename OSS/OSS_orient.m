function kps2 = OSS_orient( M, kps, sr, nb, sr_mask )
T = 0.8;
if sr == 4
    n_grid = [ 4, 0;3, 1;3, 2;2, 3;1, 3;
        0, 4; - 1, 3; - 2, 3; - 3, 2; - 3, 1;
        - 4, 0; - 3,  - 1; - 3,  - 2; - 2,  - 3; - 1,  - 3;
        0,  - 4;1,  - 3;2,  - 3;3,  - 2;3,  - 1;4, 0 ];
elseif sr == 5
    n_grid = [ 5, 0;4, 1;4, 2;4, 3;3, 3;3, 4;2, 4;1, 4;
        0, 5; - 1, 4; - 2, 4; - 3, 4; - 3, 3; - 4, 3; - 4, 2; - 4, 1;
        - 5, 0; - 4,  - 1; - 4,  - 2; - 4,  - 3; - 3,  - 3; - 3,  - 4; - 2,  - 4; - 1,  - 4;
        0,  - 5;1,  - 4;2,  - 4;3,  - 4;3,  - 3;4,  - 3;4,  - 2;4,  - 1;5, 0 ];
else
    n_grid = [ 3, 0;2, 1;2, 2;1, 2;
        0, 3; - 1, 2; - 2, 2; - 2, 1;
        - 3, 0; - 2,  - 1; - 2,  - 2; - 1,  - 2;
        0,  - 3;1,  - 2;2,  - 2;2,  - 1;3, 0 ];
end
%��Ե��Ƕ�
n_angle = atan2( n_grid( :, 2 ), n_grid( :, 1 ) ) * 180 / pi;

n_angle( n_angle < 0 ) = n_angle( n_angle < 0 ) + 360;
n_angle( end  ) = 360;

kpn = 0;
kps2 = zeros( size( kps, 1 ), 4 );

n = nb;
for i = 1:size( kps, 1 )%�����㼯��
    c0 = round( kps( i, 1 ) );%ͼ���е�λ��
    r0 = round( kps( i, 2 ) );
    %����������λ�ü����һ���ֲ����򣬲�����ת��������������
    idx = sr_mask( sub2ind( size( sr_mask ), n_grid( :, 2 ) + sr + 1, n_grid( :, 1 ) + sr + 1 ) );
    idx2 = sub2ind( size( M ), repmat( r0, length( idx ), 1 ), repmat( c0, length( idx ), 1 ), idx );
    n_bins = M( idx2 );%��6�и���������������һ��Ԥ�ȼ���õ�ֱ��ͼ M �л�ȡһ��ֱ��ͼֵ��
    bins = interp1( n_angle, n_bins, 0:360 / nb:360 - 360 / nb, 'pchip' );%������ֱ��ͼֵ��ֵ��һ�� 36�� bin ��ֱ��ͼ��
    hist = ( max( bins ) - bins ) ./ ( max( bins ) - min( bins ) );
    %�����һ����һ����ֱ��ͼ������������ǰ������ľֲ�ͼ������
    [ val, idx ] = sort( hist, 'descend' );
    val_th = val( 1 ) * T;
    count = 0;
    for j = 1:length( idx )
        kk = idx( j );
        if kk == 1
            k1 = nb;k2 = 2;
        elseif kk == nb
            k1 = nb - 1;k2 = 1;
        else
            k1 = kk - 1;k2 = kk + 1;
        end
        if ( hist( kk ) > hist( k1 ) && hist( kk ) > hist( k2 ) && hist( kk ) > val_th )
            bin = kk - 1 + 0.5 * ( hist( k1 ) - hist( k2 ) ) / ( hist( k1 ) + hist( k2 ) - 2 * hist( kk ) );
            if ( bin < 0 )
                bin = nb + bin;
            elseif ( bin >= nb )
                bin = bin - nb;
            end
            kpn = kpn + 1;
            kps2( kpn, 1 ) = kps( i, 1 );
            kps2( kpn, 2 ) = kps( i, 2 );
            kps2( kpn, 3 ) = kps( i, 3 );
            kps2( kpn, 4 ) = ( 360 / nb ) * bin;
            count = count + 1;
        end
        if count == 2
            break ;
        end
    end
end
end
% �������һ����ֱ��ͼ����һЩ�������������ֱ��ͼ������
% ������ֵ�ȵȡ��������һ���������򽫵�ǰ�����㼰���Ӧ��ֱ��
% ͼ bin ���浽 kps2 �����У������͵õ���һ�������������������