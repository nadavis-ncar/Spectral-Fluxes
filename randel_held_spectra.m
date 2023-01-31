function [C WAVENUM LAT RH_SPECTRA] = randel_held_spectra(X,Y,LAT,NumberDays,Windowing,time_res,varargin)

%This function, and all associated sub-functions, were given to me by Elizabeth Barnes
%as part of the Objective Analysis course at Colorado State University. They have been adapted
%to operate on non-cyclic grids (with windowing), parallelized, etc.

if nargin>6
   lon_fraction=varargin{1};
else
   lon_fraction=1;
end

[EAST3 WEST3] = loop_estimates(X,Y,NumberDays,Windowing);

C =[0.0:4:400]';
FREQ=[0:size(EAST3,1)-1]/(size(EAST3,1)-1)*0.5;
WAVENUM=[0:size(EAST3,2)-1];

EAST3_c(:,:,:) = nan(length(LAT),length(C),length(WAVENUM));
WEST3_c(:,:,:) = nan(length(LAT),length(C),length(WAVENUM));

%Parallelize the spectral analysis
parfor ilat = 1:length(LAT)

    if(abs(LAT(ilat))==90)
        continue
    end
    
    [K_e_new1, K_w_new1] = convert_fk_to_ck(FREQ, WAVENUM, cosd(LAT(ilat)), C, squeeze(EAST3(:,:,ilat)), squeeze(WEST3(:,:,ilat)), time_res, lon_fraction);
    
    EAST3_c(ilat,:,:) = K_e_new1;
    WEST3_c(ilat,:,:) = K_w_new1;
    
end
EAST3_c=permute(EAST3_c,[2 3 1]);
WEST3_c=permute(WEST3_c,[2 3 1]);

C = [-flipud(C);C];
RH_SPECTRA = cat(1,flipdim(WEST3_c,1),EAST3_c);
