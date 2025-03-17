function [G,R] = ImgWarp( I1, I2, P1, P2, Opt, Opt2 )

CropP = Opt2.CropP;
Disp = Opt2.Disp;
nm = 2000;

I1 = single( I1 );Siz1 = size( I1 );
I2 = single( I2 );Siz2 = size( I2 );

TP1 = Fit2DModel( P2, P1, Opt );

xb2 = [ 1, Siz2( 2 ), Siz2( 2 ), 1 ];
yb2 = [ 1, 1, Siz2( 1 ), Siz2( 1 ) ];
TP = Fit2DModel( P1, P2, Opt );
pT = Make2DModel( [ xb2', yb2' ], TP, Opt.GeoMdl );
xb1 = pT( :, 1 )';
yb1 = pT( :, 2 )';

gBG = 0;

if ~CropP
    
    if strcmp( Opt.GeoMdl, 'LPiecewise' )
        OptPI.GeoMdl = 'Projective';
        TPpI = Fit2DModel( P1, P2, OptPI );
        pTI = Make2DModel( [ xb2, yb2 ], TPpI, OptPI.GeoMdl );
        xb1 = pTI( :, 1 );yb1 = pTI( :, 2 );
    end
    
    xmin = floor( min( [ xb1, 1 ] ) );xmax = ceil( max( [ xb1, Siz1( 2 ) ] ) );
    ymin = floor( min( [ yb1, 1 ] ) );ymax = ceil( max( [ yb1, Siz1( 1 ) ] ) );
    [ x, y ] = meshgrid( xmin:xmax, ymin:ymax );
    xr1 = find( x( 1, : ) == 1 );xr2 = find( x( 1, : ) == Siz1( 2 ) );
    yr1 = find( y( :, 1 ) == 1 );yr2 = find( y( :, 1 ) == Siz1( 1 ) );
    [ xi, yi ] = meshgrid( 1:size( x, 2 ), 1:size( y, 1 ) );
    Sizi = size( xi );
    G = single( gBG * ones( Sizi ) );
    xi = single( xi( : ) );yi = single( yi( : ) );x = single( x( : ) );y = single( y( : ) );
    R = G;
    R( yr1:yr2, xr1:xr2 ) = I1; I1 = [  ];
else
    R = I1; I1 = [  ];
    [ x, y ] = meshgrid( 1:Siz1( 2 ), 1:Siz1( 1 ) );
    x = x( : );y = y( : );xi = x;yi = y;
    G = single( gBG * ones( Siz1( 1:2 ) ) );
    Sizi = Siz1;
end
switch Opt.GeoMdl
    case 'LPiecewise'
        TRI = TP.TRI;
        
        AfnPar = zeros( size( TRI, 1 ), 6 );
        for i = 1:size( TRI, 1 )
            L = [ P2( TRI( i, 1 ), : ), P2( TRI( i, 2 ), : ), P2( TRI( i, 3 ), : ) ]';
            A = [ 1, P1( TRI( i, 1 ), : ), 0, 0, 0;0, 0, 0, 1, P1( TRI( i, 1 ), : ); ...
                1, P1( TRI( i, 2 ), : ), 0, 0, 0;0, 0, 0, 1, P1( TRI( i, 2 ), : ); ...
                1, P1( TRI( i, 3 ), : ), 0, 0, 0;0, 0, 0, 1, P1( TRI( i, 3 ), : ); ];
            AfnPar( i, : ) = ( A' * A )\( A' * L );
        end
        OptP.GeoMdl = 'Projective';
        TPp = Fit2DModel( P2, P1, OptP );
        T = pointLocation( TRI, [ double( x ), double( y ) ] );
        
        for i = 1:size( TRI, 1 )
            
            gTRIi = T == i;
            xT = AfnPar( i, 1 ) + AfnPar( i, 2 ) * x( gTRIi ) + AfnPar( i, 3 ) * y( gTRIi );
            yT = AfnPar( i, 4 ) + AfnPar( i, 5 ) * x( gTRIi ) + AfnPar( i, 6 ) * y( gTRIi );
            P = fix( xT );Q = fix( yT );a = xT - P;b = yT - Q;
            Indx1 = sub2ind( Siz2, Q, P );
            Indx2 = Indx1 + Siz2( 1 );Indx3 = Indx1 + 1;Indx4 = Indx2 + 1;
            N1 = I2( Indx1 );N2 = I2( Indx2 );N3 = I2( Indx3 );N4 = I2( Indx4 );
            
            Indx = sub2ind( Sizi, yi( gTRIi ), xi( gTRIi ) );
            G( Indx ) = ( 1 - b ) .* ( ( 1 - a ) .* N1 + a .* N2 ) + b .* ( ( 1 - a ) .* N3 + a .* N4 );

        end
        xNan = x( isnan( T ) );yNan = y( isnan( T ) );xin = xi( isnan( T ) );yin = yi( isnan( T ) );
        pTNan = Make2DModel( [ xNan, yNan ], TPp, OptP.GeoMdl );
        xT = single( pTNan( :, 1 ) );yT = single( pTNan( :, 2 ) );
        
        P = fix( xT );Q = fix( yT );
        gRei = ( P < 1 | P + 1 > Siz2( 2 ) | Q < 1 | Q + 1 > Siz2( 1 ) );
        P( gRei ) = [  ];Q( gRei ) = [  ];xT( gRei ) = [  ];yT( gRei ) = [  ];xin( gRei ) = [  ];yin( gRei ) = [  ];
        a = xT - P;b = yT - Q;Indx1 = sub2ind( Siz2, Q, P );
        Indx2 = Indx1 + Siz2( 1 );Indx3 = Indx1 + 1;Indx4 = Indx2 + 1;
        N1 = I2( Indx1 );Indx1 = [  ];%#ok<NASGU>
        N2 = I2( Indx2 );Indx2 = [  ];%#ok<NASGU>
        N3 = I2( Indx3 );Indx3 = [  ];%#ok<NASGU>
        N4 = I2( Indx4 );Indx4 = [  ];%#ok<NASGU>
        
        Indx = sub2ind( Sizi, yin, xin );
        P = [  ];Q = [  ];Indx1 = [  ];Indx2 = [  ];Indx3 = [  ];Indx4 = [  ];%#ok<NASGU>
        G( Indx ) = ( 1 - b ) .* ( ( 1 - a ) .* N1 + a .* N2 ) + b .* ( ( 1 - a ) .* N3 + a .* N4 );
        
    otherwise
        n = length( x );
        Sx = fix( n / nm );
        xT = zeros( n, 1, 'single' );yT = xT;
        
        for i = 1:nm
            if i ~= nm
                xx = x( ( i - 1 ) * Sx + 1:i * Sx );
                yy = y( ( i - 1 ) * Sx + 1:i * Sx );
                Sxi = Sx;
            else
                Sxi = n - ( i - 1 ) * Sx;
                xx = x( ( i - 1 ) * Sx + 1:end  );
                yy = y( ( i - 1 ) * Sx + 1:end  );
            end
            
            pT = Make2DModel( [ xx, yy ], TP1, Opt.GeoMdl );
            
            xTi = pT( :, 1 );
            yTi = pT( :, 2 );
            
            if i ~= nm
                xT( ( i - 1 ) * Sx + 1:i * Sx ) = xTi;
                yT( ( i - 1 ) * Sx + 1:i * Sx ) = yTi;
            else
                xT( ( i - 1 ) * Sx + 1:end  ) = xTi;
                yT( ( i - 1 ) * Sx + 1:end  ) = yTi;
            end
        end
        x = [  ]; y = [  ];
        P = fix( xT );
        Q = fix( yT );
        gRei = ( P < 1 | P + 1 > Siz2( 2 ) | Q < 1 | Q + 1 > Siz2( 1 ) );
        P( gRei ) = [  ];Q( gRei ) = [  ];xT( gRei ) = [  ];yT( gRei ) = [  ];xi( gRei ) = [  ];yi( gRei ) = [  ];
        a = xT - P;
        b = yT - Q;
        
        Indx1 = sub2ind( Siz2, Q, P );
        Indx2 = Indx1 + Siz2( 1 );
        Indx3 = Indx1 + 1;
        Indx4 = Indx2 + 1;
        
        N1 = I2( Indx1 );Indx1 = [  ];%#ok<NASGU>
        N2 = I2( Indx2 );Indx2 = [  ];%#ok<NASGU>
        N3 = I2( Indx3 );Indx3 = [  ];%#ok<NASGU>
        N4 = I2( Indx4 );Indx4 = [  ];%#ok<NASGU>
        
        Indx = sub2ind( Sizi, yi, xi );
        P = [  ];Q = [  ];Indx1 = [  ];Indx2 = [  ];Indx3 = [  ];Indx4 = [  ];%#ok<NASGU>
        G( Indx ) = ( 1 - b ) .* ( ( 1 - a ) .* N1 + a .* N2 ) + b .* ( ( 1 - a ) .* N3 + a .* N4 );
end

if Disp
    
    DCel = Opt2.DCel;
    f = figure( 'Name', 'Registered Image', 'color', 'w' );
    imshow( G, [  ] );hold on
    h = imshow( R, [  ] );
    set( h, 'AlphaData', 0.4 )
    Siz = size( R );
    RegImg = zeros( Siz, 'single' );
    PartX = fix( size( R, 2 ) / DCel );PartY = fix( size( R, 1 ) / DCel );
    if rem( PartX, 2 ) == 0;PartX = PartX - 1;end ;
    if rem( PartY, 2 ) == 0;PartY = PartY - 1;end ;
    GridSize( 1 ) = fix( size( R, 2 ) / PartX );
    GridSize( 2 ) = fix( size( R, 1 ) / PartY );
    ImgVeiw = 1;xmin = 1;ymin = 1;
    for i = 1:PartX
        if i ~= PartX
            xmax = xmin + GridSize( 1 );
        else
            xmax = Siz( 2 );
        end
        for j = 1:PartY
            if j ~= PartY
                ymax = ymin + GridSize( 2 );
            else
                ymax = Siz( 1 );
            end
            
            if ImgVeiw == 1
                RegImg( ymin:ymax, xmin:xmax ) = G( ymin:ymax, xmin:xmax );
                ImgVeiw = 2;
            else
                RegImg( ymin:ymax, xmin:xmax ) = R( ymin:ymax, xmin:xmax );
                ImgVeiw = 1;
            end
            ymin = ymax;
        end
        xmin = xmax;
        ymin = 1;
    end
    set( f, 'Units', 'normalized', 'Position', [ 0.005, .04, .495, .7 ] );
end


