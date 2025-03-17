function dMGB  = Dense_MGGabor(I)
[~, Nc, ~] = size(I);
gamma = 1; 
b = 1; 
Theta = 0:pi/6:pi-pi/6; 
phi = 0; 
% ----------------------------
J = (2.^(0:log2(1)) - .5) ./ Nc;
F = [ (.25 - J) (.25 + J) ]; 
F = sort(F); 
Lambda = 1 ./ F;
% ----------------------------
FeatureImages= GaborTexture(I, gamma, Lambda, b, Theta, phi);
sigma =0.8;
f = fspecial( 'gaussian', max( 1, fix( 6 * sigma + 1 ) ), sigma );
f = cat( 3, 1 * f, 3 * f, 1 * f );
dMSOG1 = convn(FeatureImages, f, 'same' );
sum1 = sum(dMSOG1, 3 );
dMGB = dMSOG1 ./ ( sum1 + 0.000000001 );
% figure, imshow(dMGB(:,:,1),[]);
end

function featuresDense = GaborTexture(I, gamma, Lambda, b, Theta, phi)
%%  Step 1. Gabor Filter bank
i = 0;
for lambda = Lambda
    for theta = Theta
        i = i + 1;
        D = gabor2(I, gamma, lambda, b, theta, phi);
        O(:, :, i) = D;
    end
end
featuresDense = O( :,:, 2:3:end  );
end

%%    GF is Gabor Filter
function [GOG, GF] = gabor2(I, gamma, lambda, b, theta, phi)
% if nargin < 7, shape = 'same'; end;
if isa(I, 'double') ~= 1, I = double(I); end
sigma = (1 / pi) * sqrt(log(2)/2) * (2^b+1) / (2^b-1) * lambda;
Sy = sigma * gamma;
for x = -fix(sigma):fix(sigma)
    for y = -fix(Sy):fix(Sy)
        xp = x * cos(theta) + y * sin(theta);
        yp = y * cos(theta) - x * sin(theta);
        % GF is Gabor Filter
        GF(fix(Sy)+y+1,fix(sigma)+x+1) = ...
            exp(-.5*(xp^2+gamma^2*yp^2)/sigma^2) * cos(2*pi*xp/lambda+phi);
    end
end
GO_x_1 =imfilter(I(:,:,1),double(GF),'replicate');
GO_y_1 =imfilter(I(:,:,1),double(GF'),'replicate');
GO_G=sqrt(GO_x_1.^2+GO_y_1.^2);
a=max(GO_G(:)); b=min(GO_G(:)); GOG=(GO_G-b)/(a-b);%¹éÒ»»¯²Ù×÷
end