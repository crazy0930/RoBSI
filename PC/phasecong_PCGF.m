% This function is revised based on the phase congruency code (http://www.peterkovesi.com/matlabfns/) of  Peter Kovesi by
% 
%Usage: [pc or] = phasecong(im, nscale, norient)
%

%
%
% Returned values:
%                 pc - Phase congruency magnitude image (values between 0 and 1)   
%                 or - Phase congruency Orientation image.  
%
% Parameters:  
%                 im - A greyscale image to be processed.
%                 nscale - the number of log gabor scale
%                 orient - the number of log gabor orientation

%


function[M,m,phaseCongruency,or]=phasecong_PCGF(im, nscale, norient)

sze = size(im);

if nargin < 2
    nscale          = 3;     % Number of wavelet scales.
end
if nargin < 3
    norient         = 6;     % Number of filter orientations.
end
if nargin < 4
    noiseMode = 1;
end
if nargin < 5
    minWaveLength   = 3;     % Wavelength of smallest scale filter.
end
if nargin < 6
    mult            = 2.0;     % Scaling factor between successive filters.

end
if nargin < 7
    sigmaOnf        = 0.55;  % Ratio of the standard deviation of the
                             % Gaussian describing the log Gabor filter's transfer function 
			     % in the frequency domain to the filter center frequency.

end
if nargin < 8
    dThetaOnSigma   = 1.7;   % Ratio of angular interval between filter orientations
			     % and the standard deviation of the angular Gaussian
			     % function used to construct filters in the
                             % freq. plane.
end
if nargin < 9
    k               = 3.0;   % No of standard deviations of the noise energy beyond the
			     % mean at which we set the noise threshold point.
			     % standard deviation to its maximum effect
                             % on Energy.
end
if nargin < 10
    cutOff          = 0.4;   % The fractional measure of frequency spread
                             % below which phase congruency values get penalized.
end
   
g               = 10;    % Controls the sharpness of the transition in the sigmoid
                         % function used to weight phase congruency for frequency
                         % spread.
epsilon         = .0001; % Used to prevent division by zero.


thetaSigma = pi/norient/dThetaOnSigma;  % Calculate the standard deviation of the
                                        % angular Gaussian function used to
                                        % construct filters in the freq. plane.

imagefft = single(fft2(im));                    % Fourier transform of image
sze = size(imagefft);
rows = sze(1);
cols = sze(2);
zero = single(zeros(sze));

totalEnergy = zero;                     % Matrix for accumulating weighted phase 
                                        % congruency values (energy).
totalSumAn  = zero;                     % Matrix for accumulating filter response
                                        % amplitude values.
% orientation = zero;                     % Matrix storing orientation with greatest
                                        % energy for each pixel.
estMeanE2n = [];

covx2 = zero;                     % Matrices for covariance data
covy2 = zero;
covxy = zero;
%accumlate the Energy for caculating the featrue orietation
EnergyV2 = zero;

EnergyV3 = zero;

% Pre-compute some stuff to speed up filter construction

x = single(ones(rows,1) * (-cols/2 : (cols/2 - 1))/(cols/2));  
y = single((-rows/2 : (rows/2 - 1))' * ones(1,cols)/(rows/2));
radius = sqrt(x.^2 + y.^2);       % Matrix values contain *normalised* radius from centre.
radius(round(rows/2+1),round(cols/2+1)) = 1; % Get rid of the 0 radius value in the middle 
                                             % so that taking the log of the radius will 
                                             % not cause trouble.
theta = single(atan2(-y,x));              % Matrix values contain polar angle.
                                  % (note -ve y is used to give +ve
                                  % anti-clockwise angles)
sintheta = sin(theta);
costheta = cos(theta);
clear x; clear y; clear theta; clear im;     % save a little memory

% The main loop...

for o = 1:norient,                   % For each orientation.
  %disp(['Processing orientation ' num2str(o)]);
  angl = (o-1)*pi/norient;           % Calculate filter angle.
  wavelength = minWaveLength;        % Initialize filter wavelength.
  sumE_ThisOrient   = zero;          % Initialize accumulator matrices.
  sumO_ThisOrient   = zero;  
  sumO_ThisOrient1   = zero;        % the odd energy expect the smallest wave scale
  sumAn_ThisOrient  = zero;      
  Energy_ThisOrient = zero;      
  EOArray = single([]);          % Array of complex convolution images - one for each scale.
  ifftFilterArray = single([]);  % Array of inverse FFTs of filters

  % Pre-compute filter data specific to this orientation
  % For each point in the filter matrix calculate the angular distance from the
  % specified filter orientation.  To overcome the angular wrap-around problem
  % sine difference and cosine difference values are first computed and then
  % the atan2 function is used to determine angular distance.

  ds = sintheta * cos(angl) - costheta * sin(angl); % Difference in sine.
  dc = costheta * cos(angl) + sintheta * sin(angl); % Difference in cosine.
  dtheta = abs(atan2(ds,dc));                           % Absolute angular distance.
  spread = exp((-dtheta.^2) / (2 * thetaSigma^2));      % Calculate the angular filter component.

  clear ds;clear dc;clear dtheta;
  for s = 1:nscale,                  % For each scale.

    % Construct the filter - first calculate the radial filter component.
    fo = 1.0/wavelength;                  % Centre frequency of filter.
    rfo = fo/0.5;                         % Normalised radius from centre of frequency plane 
    %rfo = 0.3;                                     % corresponding to fo.
    logGabor = exp((-(log(radius/rfo)).^2) / (2 * log(sigmaOnf)^2));  
    %logGabor = logGabor.*lp;
    logGabor(round(rows/2+1),round(cols/2+1)) = 0; % Set the value at the center of the filter
                                                   % back to zero (undo the radius fudge).

    filter = logGabor .* spread;          % Multiply by the angular spread to get the filter.
    filter = fftshift(filter);            % Swap quadrants to move zero frequency 
                                          % to the corners.
    clear logGabor;

    ifftFilt = real(ifft2(filter))*sqrt(rows*cols);  % Note rescaling to match power
    ifftFilterArray = single([ifftFilterArray ifftFilt]);    % record ifft2 of filter
    
    clear ifftFilt;
    
    
    % Convolve image with even and odd filters returning the result in EO
    EOfft = imagefft .* filter;           % Do the convolution.
    EO = single(ifft2(EOfft));                    % Back transform.
    
    clear EOfft;
    

    EOArray = single([EOArray, EO]);              % Record convolution result
    An = abs(EO);                         % Amplitude of even & odd filter response.

    sumAn_ThisOrient = single(sumAn_ThisOrient + An);     % Sum of amplitude responses.
    sumE_ThisOrient = single(sumE_ThisOrient + real(EO)); % Sum of even filter convolution results.
    sumO_ThisOrient = single(sumO_ThisOrient + imag(EO)); % Sum of odd filter convolution results.
    if s>1;
        sumO_ThisOrient1 = single(sumO_ThisOrient1 + imag(EO));%sum the odd amplitude response except the smallest scale
    end
    if s == 1
       maxSumO = sumO_ThisOrient; %Record the maximum odd amplitude responses
    else
        maxSumO = max(maxSumO,sumO_ThisOrient);
    end
    if s == 1                             % Record the maximum An over all scales
      maxAn = An;
    else
      maxAn = max(maxAn, An);
    end
    
    if s==1
      EM_n = sum(sum(filter.^2));           % Record mean squared filter value at smallest
    end                                     % scale. This is used for noise estimation.

    wavelength = wavelength * mult;% Finally calculate Wavelength of next filter
    
    clear An; clear filter;
  end                                       % ... and process the next scale

  % Get weighted mean filter response vector, this gives the weighted mean phase angle.

  XEnergy = sqrt(sumE_ThisOrient.^2 + sumO_ThisOrient.^2) + epsilon;   
  MeanE = sumE_ThisOrient ./ XEnergy; 
  MeanO = sumO_ThisOrient ./ XEnergy; 
  
  clear XEnergy;

  % Now calculate An(cos(phase_deviation) - | sin(phase_deviation)) | by using
  % dot and cross products between the weighted mean filter response vector and
  % the individual filter response vectors at each scale.  This quantity is 
  % phase congruency multiplied by An, which we call energy.

  for s = 1:nscale,       
      EO = submat(EOArray,s,cols);  % Extract even and odd filter 
      E = real(EO); O = imag(EO);
      Energy_ThisOrient = Energy_ThisOrient ...
        + E.*MeanE + O.*MeanO - abs(E.*MeanO - O.*MeanE);
  end
  clear EO;clear E; clear O;clear MeanE; clear MeanO;

  % Note: To calculate the phase symmetry measure replace the for loop above 
  % with the following loop. (The calculation of MeanE, MeanO, sumE_ThisOrient 
  % and sumO_ThisOrient can also be omitted). It is suggested that the value
  % of nscale is increased (to say, 5 for a 256x256 image) and that cutOff is
  % set to 0 to eliminate weighting for frequency spread.

%   for s = 1:nscale,                  
%     Energy_ThisOrient = Energy_ThisOrient ...
%      + abs(real(submat(EOArray,s,cols))) - abs(imag(submat(EOArray,s,cols)));
%   end

  % Compensate for noise
  % We estimate the noise power from the energy squared response at the smallest scale.
  % If the noise is Gaussian the energy squared will have a Chi-squared 2DOF pdf.
  % We calculate the median energy squared response as this is a robust statistic.  
  % From this we estimate the mean.  
  % The estimate of noise power is obtained by dividing the mean squared energy value
  % by the mean squared filter value

  medianE2n = median(reshape(abs(submat(EOArray,1,cols)).^2,1,rows*cols));
  meanE2n = -medianE2n/log(0.5);
  estMeanE2n = [estMeanE2n meanE2n];

  noisePower = meanE2n/EM_n;                       % Estimate of noise power.
  
  clear meanE2n;clear medianE2n; clear meanE2n;

  % Now estimate the total energy^2 due to noise
  % Estimate for sum(An^2) + sum(Ai.*Aj.*(cphi.*cphj + sphi.*sphj))

  EstSumAn2 = zero;
  for s = 1:nscale
    EstSumAn2 = EstSumAn2+submat(ifftFilterArray,s,cols).^2;
  end

  EstSumAiAj = zero;
  for si = 1:(nscale-1)
    for sj = (si+1):nscale
      EstSumAiAj = EstSumAiAj + submat(ifftFilterArray,si,cols).*submat(ifftFilterArray,sj,cols);
    end
  end

  EstNoiseEnergy2 = 2*noisePower*sum(sum(EstSumAn2)) + 4*noisePower*sum(sum(EstSumAiAj));
  
  clear EstSumAn2;
  tau = sqrt(EstNoiseEnergy2/2);                     % Rayleigh parameter
  EstNoiseEnergy = tau*sqrt(pi/2);                   % Expected value of noise energy
  EstNoiseEnergySigma = sqrt( (2-pi/2)*tau^2 );

  T =  EstNoiseEnergy + k*EstNoiseEnergySigma;       % Noise threshold
  
  clear EstNoiseEnergy; clear EstNoiseEnergySigma; clear tau;
  clear EstNoiseEnergy2;clear EstSumAiAj;clear noisePower; 

  % The estimated noise effect calculated above is only valid for the PC_1 measure. 
  % The PC_2 measure does not lend itself readily to the same analysis.  However
  % empirically it seems that the noise effect is overestimated roughly by a factor 
  % of 1.7 for the filter parameters used here.

  T = T/1.7;        % Empirical rescaling of the estimated noise effect to 
                    % suit the PC_2 phase congruency measure
  %T = 0;

  Energy_ThisOrient = max(Energy_ThisOrient - T, zero);  % Apply noise threshold

  % Form weighting that penalizes frequency distributions that are particularly
  % narrow.
  % Calculate fractional 'width' of the frequencies present by taking
  % the sum of the filter response amplitudes and dividing by the maximum 
  % amplitude at each point on the image.

  width = sumAn_ThisOrient ./ (maxAn + epsilon) / nscale;    

  % Now calculate the sigmoidal weighting function for this orientation.

  weight = 1.0 ./ (1 + exp( (cutOff - width)*g)); 

  % Apply weighting

  Energy_ThisOrient =   weight.*Energy_ThisOrient;
  
  clear weight;clear width;

  % Update accumulator matrix for sumAn and totalEnergy

  totalSumAn  = totalSumAn + sumAn_ThisOrient;%分母
  totalEnergy = totalEnergy + Energy_ThisOrient;%分子
  
  %caculate the orientated pc
  PC{o} = Energy_ThisOrient./sumAn_ThisOrient;
  
  % Build up covariance data for every point
  covx = PC{o}*cos(angl);
  covy = PC{o}*sin(angl);
  covx2 = covx2 + covx.^2;
  covy2 = covy2 + covy.^2;
  covxy = covxy + covx.*covy;

   % 生成HARRIS函数
%     harris_function=covx2 .*covy2-covxy.*covxy-0.04*(covx2+covy2).^2;
%     harris_function=harris_function*(180/norient);

  % Update orientation matrix by finding image points where the energy in this
  % orientation is greater than in any previous orientation (the change matrix)
  % and then replacing these elements in the orientation matrix with the
  % current orientation number.

EnergyV2 = EnergyV2 + cos(angl)*sumO_ThisOrient;
EnergyV3 = EnergyV3 + sin(angl)*sumO_ThisOrient;
  
  clear sumAn_ThisOrient; clear Energy_ThisOrient; clear sumO_ThisOrient;
  clear sumO_ThisOrient; clear spread; clear EOArray; clear ifftFilterArray;

  if(o == 1),
   % maxEnergy = Energy_ThisOrient;
   % featType = E + i*O;
  else
    %change = Energy_ThisOrient > maxEnergy;
    %orientation = (o - 1).*change + orientation.*(~change);
    %featType = (E+i*O).*change + featType.*(~change);
    %maxEnergy = max(maxEnergy, Energy_ThisOrient);
  end

end  % For each orientation

    % First normalise covariance values by the number of orientations/2
    %  加权最大力矩图
    covx2 = covx2/(norient/2);
    covy2 = covy2/(norient/2);
    covxy = 4*covxy/norient;   % This gives us 2*covxy/(norient/2)
    denom = sqrt(covxy.^2 + (covx2-covy2).^2)+epsilon;
    M = (covy2+covx2 + denom)/2;          % Maximum moment
    m = (covy2+covx2 - denom)/2;           % Minimum moment

phaseCongruency = double(totalEnergy ./ (totalSumAn + epsilon));
%used the atan (EnergyV3,-EnergyV2) to abain phase cogruency magnitude 

or = double(atan(EnergyV3./(-EnergyV2)));%note the y direction is reverse


function a = submat(big,i,cols)

a = big(:,((i-1)*cols+1):(i*cols));

