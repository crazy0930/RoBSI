function fea = OSS_descriptor( nM, nkps, sr, nb, o, PS, br, ba )
sr_mask = strel( 'disk', sr, 0 ).Neighborhood;
sr_mask( sr + 1, sr + 1 ) = false;
sr_idx = find( sr_mask );
sr_mask = zeros( size( sr_mask ) );
sr_mask( sr_idx ) = 1:length( sr_idx );

fkps = [  ];
fdes = [  ];
L = max( nkps( :, 4 ) );
parfor i = 1:L
    kps = nkps( nkps( :, 4 ) == i, 1:3 );
    
    if size( kps, 1 ) < 1
        continue ;
    end
    kps( :, 1:2 ) = kps( :, 1:2 ) ./ repmat( kps( :, 3 ), 1, 2 );
    M = nM{ i };
    kps = OSS_orient( M, kps, sr, nb, sr_mask );
    ssf = OSS_desG( M, kps', PS, o, br, ba, sr, sr_mask );%ÃèÊö×ÓÏòÁ¿
    ssf.kps( :, 1:2 ) = ssf.kps( :, 1:2 ) .* repmat( ssf.kps( :, 3 ), 1, 2 );
    fkps = [ fkps;ssf.kps ];
    fdes = [ fdes;ssf.des ];
end
fea.kps = fkps;
fea.des = fdes;
end
