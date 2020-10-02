function [dODFilt,dOD] = edgeReferenceTransient(Ion,Ioff,Eedge)
%
% This code implements edge-pixel referencing of a transient absorption
% experiment, as described in R. Geneaux et al., "Source noise suppression 
% in attosecond transient absorption spectroscopy by edge-pixel referencing"
% (to be published)
% -----------------
%
% This function performs the following:
% (1) computes the optical density from a set of pump-off and
% pump-on spectra;
% (2) Uses the pump-off spectra to characterize the correlation of the
% light spectrum, by calculating the B matrix using the user-defined
% edge-pixels;
% (3) Applies the B matrix to suppress source in the dataset.
%
% Inputs (required)
%   Ion - set of spectra taken with the pump on, of size (Nscan x Ntime x Nenergy)
%              Nscan = Number of averages
%              Ntime = Number of time points
%              Nenergy = Number of pixels in the energy domain;
%   Ioff - same with the pump off;
%   Eedge: vector of boolean of size Nenergy wherein Eedge == 1 corresponds
%   to the edge-pixel region (where there is no transient signal).
%
% Outputs:
%   dODFilt - Filtered transient absorption(size Ntime x Nenergy);
%   dOD - Raw change in optical absorption
% ---------------------------------------------
%
%   Example:
%   % E is the energy axis in eV
%   % There is no signal between 17.5 and 20 eV, and between 48 and 50 eV,
%   % these are the edge-pixel regions
%   Eedge = ((E<20)&(E>17.5))|((E<50)&(E>48));
%
%   [dODFilt,dOD] = edgeReferenceTransient(Ion,Ioff,Eedge)
%
% Citation for this code or some of its parts: R. Geneaux et al. "Source noise 
% suppression in attosecond transient absorption spectroscopy by edge-pixel 
% referencing" (to be published)
%==========================================================================

%size checks
if any(size(Ion) ~= size(Ioff))
    error('Ion and Ioff do not have the same size');
elseif size(Eedge) ~= size(Ion,3)
    error('Eedge does not have the same size as Ion along the energy dimension')
else
    
    % Define raw dOD
    if size(Ion,1) == 1  %if there is only one average
        dOD = squeeze(-real(log10(Ion./Ioff)));
    else
        dOD = squeeze(mean(-real(log10(Ion./Ioff))));
    end
    dOD(isnan(dOD)) = 0; dOD(isinf(dOD)) = 0;
    
    % Use all shots with pump off to construct an OD calibration set for the noise
    allShotsPumpOff = reshape(permute(Ioff,[2 1 3]),size(Ioff,1)*size(Ioff,2),[]);
    calibOD = -real(log10(allShotsPumpOff(1:2:(end-1),:)./allShotsPumpOff(2:2:end,:)));
    
    % define signal and edge zones
    edgeZone = calibOD(:,Eedge);
    signalZone = calibOD;
    
    %check if there is more calibration points than edge pixels
    [p,m] = size(edgeZone);
    if p <= m
        msg = sprintf(['There are less calibration measurements (' num2str(p) ') than edge-pixels (' num2str(m)...
            '). Reduce size of edge-pixel region or use more calibration points. \nThe covariance matrix does not have full rank (RCond = ' num2str(rcond(cov(edgeZone))) '), results will be inaccurate.']);
        warning(msg);
        warning('off','MATLAB:nearlySingularMatrix') %turn off redundant default warning
    end
            
    
    %compute B matrix
    B=bsxfun(@minus,signalZone,mean(signalZone))'*bsxfun(@minus,edgeZone,mean(edgeZone))./(size(signalZone,1)-1)/cov(edgeZone);
    
    %Apply B
    Correction = dOD(:,Eedge)*B';
    dODFilt = dOD - Correction;
    warning('on','MATLAB:nearlySingularMatrix') %turn default warning back on
end
end

