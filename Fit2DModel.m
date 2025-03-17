function TP = Fit2DModel( P1, P2, Opt )
%% 参数
%  Siz1 = size( I1 );
%  Disp.Disp = 1;
%  Disp.DCel = fix( max( Siz1 ) / 10 );
%  Disp.CropP = 0;

% Model = { 'Projective', 'Polynomial', 'LPiecewise', 'MQuadrics', 'PointWise' };
% Opt.GeoMdl = Model{ 5 };
% Opt.n = 2;
% Opt.d = 0;
% Opt.k = 5;
% Opt.pw = 2;

switch Opt.GeoMdl
    case 'Similarity'
        TP = Similarity( P1, P2 );
        
    case 'Affine'
        TP = Affine( P1, P2 );
        
    case 'Bilinear'
        TP = Bilinear( P1, P2 );
        
    case 'Projective'
        TP = Projective( P1, P2 );
        
    case 'Polynomial'
        n = Opt.n;
        TP = Polynomial( P1, P2, n );
        
    case 'LPiecewise'
        TP = LPiecewise( P1, P2 );
        
    case 'MQuadrics'
        n = Opt.n;
        d = Opt.d;
        TP = MQuadrics( P1, P2, d, n );
        
    case 'TPS'
        n = Opt.n;
        d = Opt.d;
        TP = TPS( P1, P2, d, n );
        
    case 'PW'
        n = Opt.n;
        d = Opt.d;
        k = Opt.k;
        n2 = Opt.n2;
        TP = PW( P1, P2, d, n, k, n2 );
        
    case 'PointWise'
        n = Opt.n;
        d = Opt.d;
        k = Opt.k;
        pw = Opt.pw;
        TP = PointWise( P1, P2, d, n, k, pw );
        
    case 'WeightedMean'
        n = Opt.n;
        d = Opt.d;
        pw = Opt.pw;
        TP = WeightedMean( P1, P2, d, n, pw );
end

function TP = Similarity( p1, p2 )
N = size( p1, 1 );
L = zeros( 2 * N, 1 );
A = zeros( 2 * N, 4 );

A( 1:2:end  - 1, 1 ) = 1;
A( 1:2:end  - 1, 3 ) = p2( :, 1 );
A( 1:2:end  - 1, 4 ) = p2( :, 2 );
A( 2:2:end , 2 ) = 1;
A( 2:2:end , 3 ) = p2( :, 2 );
A( 2:2:end , 4 ) =  - p2( :, 1 );
L( 1:2:end  - 1 ) = p1( :, 1 );
L( 2:2:end  ) = p1( :, 2 );

Tp = ( A' * A )\( A' * L );
vi = A * Tp - L;
ri = sqrt( vi( 1:2:end  - 1 ) .^ 2 + vi( 2:2:end  ) .^ 2 );
RMSE = sqrt( sum( ri .^ 2 ) / ( N - 1 ) );

TP.Tp = Tp;
TP.RMSE = RMSE;
TP.M = 'Similarity';

function TP = Affine( p1, p2 )
N = size( p1, 1 );
L = zeros( 2 * N, 1 );
A = zeros( 2 * N, 6 );

A( 1:2:end  - 1, 1 ) = 1;
A( 1:2:end  - 1, 2:3 ) = p2( :, 1:2 );
A( 2:2:end , 4 ) = 1;
A( 2:2:end , 5:6 ) = p2( :, 1:2 );
L( 1:2:end  - 1 ) = p1( :, 1 );
L( 2:2:end  ) = p1( :, 2 );

Tp = ( A' * A )\( A' * L );
vi = A * Tp - L;
ri = sqrt( vi( 1:2:end  - 1 ) .^ 2 + vi( 2:2:end  ) .^ 2 );
RMSE = sqrt( sum( ri .^ 2 ) / ( N - 1 ) );

TP.Tp = Tp;
TP.RMSE = RMSE;
TP.M = 'Affine';

function TP = Bilinear( p1, p2 )
N = size( p1, 1 );
L = zeros( 2 * N, 1 );
A = zeros( 2 * N, 8 );

A( 1:2:end  - 1, 1 ) = 1;
A( 1:2:end  - 1, 2:3 ) = p2( :, 1:2 );
A( 1:2:end  - 1, 4 ) = p2( :, 1 ) .* p2( :, 2 );
A( 2:2:end , 5 ) = 1;
A( 2:2:end , 6:7 ) = p2( :, 1:2 );
A( 2:2:end , 8 ) = p2( :, 1 ) .* p2( :, 2 );
L( 1:2:end  - 1 ) = p1( :, 1 );
L( 2:2:end  ) = p1( :, 2 );

Tp = ( A' * A )\( A' * L );
vi = A * Tp - L;
ri = sqrt( vi( 1:2:end  - 1 ) .^ 2 + vi( 2:2:end  ) .^ 2 );
RMSE = sqrt( sum( ri .^ 2 ) / ( N - 1 ) );

TP.Tp = Tp;
TP.RMSE = RMSE;
TP.M = 'Bilinear';

function TP = Projective( p1, p2 )
N = size( p1, 1 );
L = zeros( 2 * N, 1 );
A = zeros( 2 * N, 8 );

A( 1:2:end  - 1, 1 ) = 1;
A( 1:2:end  - 1, 2:3 ) = p2( :, 1:2 );
A( 1:2:end  - 1, 7 ) =  - p1( :, 1 ) .* p2( :, 1 );
A( 1:2:end  - 1, 8 ) =  - p1( :, 1 ) .* p2( :, 2 );
A( 2:2:end , 4 ) = 1;
A( 2:2:end , 5:6 ) = p2( :, 1:2 );
A( 2:2:end , 7 ) =  - p1( :, 2 ) .* p2( :, 1 );
A( 2:2:end , 8 ) =  - p1( :, 2 ) .* p2( :, 2 );
L( 1:2:end  - 1 ) = p1( :, 1 );
L( 2:2:end  ) = p1( :, 2 );

Tp = ( A' * A )\( A' * L );
vi = A * Tp - L;
ri = sqrt( vi( 1:2:end  - 1 ) .^ 2 + vi( 2:2:end  ) .^ 2 );
RMSE = sqrt( sum( ri .^ 2 ) / ( N - 1 ) );

TP.Tp = Tp;
TP.RMSE = RMSE;
TP.M = 'Projective';

function TP = Polynomial( p1, p2, n )
N = size( p1, 1 );
X = p1( :, 1 );
Y = p1( :, 2 );
x = p2( :, 1 );
y = p2( :, 2 );
ui = [ 3, 6, 10, 15, 21 ];
u = ui( n );
L = zeros( 2 * N, 1 );
A = zeros( 2 * N, 42 );

Poly = [ ones( N, 1 ), x, y ...
    , x .* y, x .^ 2, y .^ 2 ...
    , ( x .^ 2 ) .* y, x .* ( y .^ 2 ), x .^ 3, y .^ 3 ...
    , ( x .^ 3 ) .* y, x .* ( y .^ 3 ), x .^ 4, y .^ 4, ( x .^ 2 ) .* ( y .^ 2 ) ...
    , ( x .^ 3 ) .* ( y .^ 2 ), ( x .^ 2 ) .* ( y .^ 3 ), x .^ 5, y .^ 5, ( x .^ 4 ) .* y, x .* ( y .^ 4 ) ];

A( 1:2:end  - 1, 1:21 ) = Poly;
A( 2:2:end , 22:end  ) = Poly;
A( :, u + 1:21 ) = [  ];
A( :, 2 * u + 1:end  ) = [  ];
L( 1:2:end  - 1 ) = X;
L( 2:2:end  ) = Y;

Tp = ( A' * A )\( A' * L );
vi = A * Tp - L;
ri = sqrt( vi( 1:2:end  - 1 ) .^ 2 + vi( 2:2:end  ) .^ 2 );
RMSE = sqrt( sum( ri .^ 2 ) / ( N - 1 ) );

TP.Tp = Tp;
TP.RMSE = RMSE;
TP.n = n;
TP.M = 'Polynomial';

function [ TP ] = LPiecewise( p1, p2 )

TRI = delaunayTriangulation( [ p1( :, 1 ), p1( :, 2 ) ] );
Tp = zeros( size( TRI, 1 ), 6 );
for i = 1:size( TRI, 1 )
    L = [ p1( TRI( i, 1 ), : ), p1( TRI( i, 2 ), : ), p1( TRI( i, 3 ), : ) ]';
    A = [ 1, p2( TRI( i, 1 ), : ), 0, 0, 0;0, 0, 0, 1, p2( TRI( i, 1 ), : ); ...
        1, p2( TRI( i, 2 ), : ), 0, 0, 0;0, 0, 0, 1, p2( TRI( i, 2 ), : ); ...
        1, p2( TRI( i, 3 ), : ), 0, 0, 0;0, 0, 0, 1, p2( TRI( i, 3 ), : ); ];
    Tp( i, : ) = ( A' * A )\( A' * L );
end

TP.Tp = Tp;
TP.TRI = TRI;
TP.RMSE = [  ];
TP.M = 'LPiecewise';

function [ TP ] = MQuadrics( p1, p2, d, n )

N = size( p1, 1 );
x = p2( :, 1 );
y = p2( :, 2 );

TPa = Polynomial( p1, p2, n );

xx = repmat( x, 1, N );
yy = repmat( y, 1, N );
Ri = ( ( xx - xx' ) .^ 2 + ( yy - yy' ) .^ 2 + d ) .^ .5;
pT = Make2DModel( p2, TPa, 'Polynomial' );

Tpx = ( Ri' * Ri )\( Ri' * ( p1( :, 1 ) - pT( :, 1 ) ) );
Tpy = ( Ri' * Ri )\( Ri' * ( p1( :, 2 ) - pT( :, 2 ) ) );

vx = Ri * Tpx;
vy = Ri * Tpy;

xt = pT( :, 1 ) + vx;yt = pT( :, 2 ) + vy;
dx = xt - p1( :, 1 );dy = yt - p1( :, 2 );
ri = sqrt( dx .^ 2 + dy .^ 2 );
RMSE = sqrt( sum( ri .^ 2 ) / ( N - 1 ) );

TP.Tpx = Tpx;
TP.Tpy = Tpy;
TP.Tp = TPa.Tp;
TP.RMSE = RMSE;
TP.x = x;
TP.y = y;
TP.d = d;
TP.n = n;
TP.M = 'MQuadrics';

function [ TP ] = MQuadricsI( p1, p2, d, n )

N = size( p1, 1 );
x = p2( :, 1 );
y = p2( :, 2 );

TPa = Polynomial( p1, p2, n );

xx = repmat( x, 1, N );
yy = repmat( y, 1, N );
Ri = exp(  - ( ( xx - xx' ) .^ 2 + ( yy - yy' ) .^ 2 + d + eps ) );

pT = Make2DModel( p2, TPa, 'Polynomial' );

Tpx = ( Ri' * Ri )\( Ri' * ( p1( :, 1 ) - pT( :, 1 ) ) );
Tpy = ( Ri' * Ri )\( Ri' * ( p1( :, 2 ) - pT( :, 2 ) ) );

vx = Ri * Tpx;
vy = Ri * Tpy;

xt = pT( :, 1 ) + vx;yt = pT( :, 2 ) + vy;
dx = xt - p1( :, 1 );dy = yt - p1( :, 2 );
ri = sqrt( dx .^ 2 + dy .^ 2 );
RMSE = sqrt( sum( ri .^ 2 ) / ( N - 1 ) );

TP.Tpx = Tpx;
TP.Tpy = Tpy;
TP.Tp = TPa.Tp;
TP.RMSE = RMSE;
TP.x = x;
TP.y = y;
TP.d = d;
TP.n = n;
TP.M = 'MQuadricsI';

function [ TP ] = TPS( p1, p2, d, n )

N = size( p1, 1 );
x = p2( :, 1 );
y = p2( :, 2 );

TPa = Polynomial( p1, p2, n );

xx = repmat( x, 1, N );
yy = repmat( y, 1, N );
Ri2 = ( xx - xx' ) .^ 2 + ( yy - yy' ) .^ 2 + d;
Ri = Ri2 .* log( Ri2 );
Ri( isnan( Ri ) ) = 0;

pT = Make2DModel( p2, TPa, 'Polynomial' );

Tpx = ( Ri' * Ri )\( Ri' * ( p1( :, 1 ) - pT( :, 1 ) ) );
Tpy = ( Ri' * Ri )\( Ri' * ( p1( :, 2 ) - pT( :, 2 ) ) );

vx = Ri * Tpx;
vy = Ri * Tpy;

xt = pT( :, 1 ) + vx;yt = pT( :, 2 ) + vy;
dx = xt - p1( :, 1 );dy = yt - p1( :, 2 );
ri = sqrt( dx .^ 2 + dy .^ 2 );
RMSE = sqrt( sum( ri .^ 2 ) / ( N - 1 ) );

TP.Tpx = Tpx;
TP.Tpy = Tpy;
TP.Tp = TPa.Tp;
TP.RMSE = RMSE;
TP.x = x;
TP.y = y;
TP.d = d;
TP.n = n;
TP.M = 'TPS';

function [ TP ] = TPSs( p1, p2, d, n )
N = size( p1, 1 );
x = p2( :, 1 );
y = p2( :, 2 );

TPa = Polynomial( p1, p2, n );

xx = repmat( x, 1, N );
yy = repmat( y, 1, N );
Ri2 = ( xx - xx' ) .^ 2 + ( yy - yy' ) .^ 2 + d;
Ri = Ri2 .* log( Ri2 );
Ri( isnan( Ri ) ) = 0;
Ri( end  + 1, : ) = 1;
Ri( end  + 1, : ) = x;
Ri( end  + 1, : ) = y;

pT = Make2DModel( p2, TPa, 'Polynomial' );

L1 = ( p1( :, 1 ) - pT( :, 1 ) );L1( end  + 3 ) = 0;
L2 = ( p1( :, 2 ) - pT( :, 2 ) );L2( end  + 3 ) = 0;
Tpx = ( Ri' * Ri )\( Ri' * L1 );
Tpy = ( Ri' * Ri )\( Ri' * L2 );

vx = Ri * Tpx;
vy = Ri * Tpy;

xt = pT( :, 1 ) + vx( 1:end  - 3 );yt = pT( :, 2 ) + vy( 1:end  - 3 );
dx = xt - p1( :, 1 );dy = yt - p1( :, 2 );
ri = sqrt( dx .^ 2 + dy .^ 2 );
RMSE = sqrt( sum( ri .^ 2 ) / ( N - 1 ) );

TP.Tpx = Tpx;
TP.Tpy = Tpy;
TP.Tp = TPa.Tp;
TP.RMSE = RMSE;
TP.x = x;
TP.y = y;
TP.d = d;
TP.n = n;

function [ TP ] = PW( p1, p2, d, n, k, n2 )

N = size( p1, 1 );

TPa = Polynomial( p1, p2, n );

pT = Make2DModel( p2, TPa, 'Polynomial' );

DMat = DistMatrix( p1 );
[ ~, Ind ] = sort( DMat );
KM = Ind( 1:k, : )';

Tpi = cell( size( KM, 1 ), 1 );
for i = 1:size( KM, 1 )
    p1i = p1( KM( i, : ), : );
    pTi = pT( KM( i, : ), : );
    TPaii = Polynomial( p1i, pTi, n2 );
    pT( KM( i, : ), : ) = Make2DModel( pTi, TPaii, 'Polynomial' );
    Tpi{ i } = TPaii;
end

dx = pT( :, 1 ) - p1( :, 1 );dy = pT( :, 2 ) - p1( :, 2 );
ri = sqrt( dx .^ 2 + dy .^ 2 );
RMSE = sqrt( sum( ri .^ 2 ) / ( N - 1 ) );

TP.Tpi = Tpi;
TP.Tp = TPa.Tp;
TP.RMSE = RMSE;
TP.n = n;
TP.d = d;
TP.k = k;
TP.n2 = n2;
TP.x = p1( :, 1 );
TP.y = p1( :, 2 );
TP.M = 'PW';

function [ TP ] = PointWise( p1, p2, d, n, k, pw )

N = size( p1, 1 );

TPa = Polynomial( p1, p2, n );

pT = Make2DModel( p2, TPa, 'Polynomial' );
dX = pT( :, 1 ) - p1( :, 1 );dY = pT( :, 2 ) - p1( :, 2 );

DMat = DistMatrix( p1 );
[ W, Ind ] = sort( DMat );
KM = Ind( 1:k, : )';
W = W( 1:k, : );

DX = zeros( N, 1 );DY = DX;
for i = 1:size( KM, 1 )
    dXi = dX( KM( i, : ) );
    dYi = dY( KM( i, : ) );
    Wi = 1 ./ ( W( :, i ) ) .^ pw;
    
    DX( i ) = sum( dXi .* Wi ) / sum( Wi );
    DY( i ) = sum( dYi .* Wi ) / sum( Wi );
end

pT( :, 1 ) = pT( :, 1 ) - DX;
pT( :, 2 ) = pT( :, 2 ) - DY;

dx = pT( :, 1 ) - p1( :, 1 );dy = pT( :, 2 ) - p1( :, 2 );
ri = sqrt( dx .^ 2 + dy .^ 2 );
RMSE = sqrt( sum( ri .^ 2 ) / ( N - 1 ) );

TP.dX = dX;
TP.dY = dY;
TP.Tp = TPa.Tp;
TP.RMSE = RMSE;
TP.n = n;
TP.d = d;
TP.k = k;
TP.pw = pw;
TP.x = p1( :, 1 );
TP.y = p1( :, 2 );
TP.M = 'PointWise';

function [ TP ] = WeightedMean( p1, p2, d, n, pw )

N = size( p1, 1 );

TPa = Polynomial( p1, p2, n );

pT = Make2DModel( p2, TPa, 'Polynomial' );
dX = pT( :, 1 ) - p1( :, 1 );dY = pT( :, 2 ) - p1( :, 2 );

W = DistMatrix( p1 );
W = 1 ./ ( W .^ pw + eps );
DX = zeros( N, 1 );DY = DX;
for i = 1:N
    DX( i ) = sum( dX .* W( :, i ) ) / sum( W( :, i ) );
    DY( i ) = sum( dY .* W( :, i ) ) / sum( W( i, : ) );
end

pT( :, 1 ) = pT( :, 1 ) - DX;
pT( :, 2 ) = pT( :, 2 ) - DY;

dx = pT( :, 1 ) - p1( :, 1 );dy = pT( :, 2 ) - p1( :, 2 );
ri = sqrt( dx .^ 2 + dy .^ 2 );
RMSE = sqrt( sum( ri .^ 2 ) / ( N - 1 ) );

TP.dX = dX;
TP.dY = dY;
TP.Tp = TPa.Tp;
TP.RMSE = RMSE;
TP.n = n;
TP.d = d;
TP.x = p1( :, 1 );
TP.y = p1( :, 2 );
TP.pw = pw;
TP.M = 'WeightedMean';
function DMat = DistMatrix( P )
pointsCount = size( P, 1 );
x = P( :, 1 );y = P( :, 2 );
x = repmat( x, 1, pointsCount );y = repmat( y, 1, pointsCount );
DMat = sqrt( ( x - x' ) .^ 2 + ( y - y' ) .^ 2 );

