function [EAST WEST] = loop_estimates(UVEL,VVEL,CHUNK_VEC,Windowing)

[TIME_SIZE LAT_SIZE LONG_SIZE] = size(UVEL);

noverlap = 0;

nrealizations = size(CHUNK_VEC,1);
Nt = length(CHUNK_VEC(1,1):CHUNK_VEC(1,2));
ham=0.54-0.46*cos(2*pi*[0:Nt-1]/(Nt-1));

UVEL = permute(UVEL,[3 1 2]);
VVEL = permute(VVEL,[3 1 2]);

EAST = nan(nrealizations,(LONG_SIZE/2+1),Nt/2+1,LAT_SIZE);
WEST = nan(nrealizations,(LONG_SIZE/2+1),Nt/2+1,LAT_SIZE);

disp('------running RH loop_estimates.m-------');

for ireal=1:nrealizations

    disp_text = sprintf('RH spectrum realization = %d/%d',ireal,nrealizations);
    disp(disp_text); 
    
    istart = CHUNK_VEC(ireal,1);
    iend = CHUNK_VEC(ireal,2);
    Nt = length(istart:iend);
    
    Ut=UVEL(:,istart:iend,:);
    Vt=VVEL(:,istart:iend,:);

    Utham=Ut.*repmat(ham,[LONG_SIZE 1 LAT_SIZE]);
    Vtham=Vt.*repmat(ham,[LONG_SIZE 1 LAT_SIZE]);

    [EAST_1 WEST_1] = space_time_new(Utham,Vtham,Windowing);
    
    EAST(ireal,:,:,:) = EAST_1;
    WEST(ireal,:,:,:) = WEST_1;

end


EAST=squeeze(mean(EAST,1,'omitnan'));
WEST=squeeze(mean(WEST,1,'omitnan'));

EAST = permute(EAST,[2 1 3]);
WEST = permute(WEST,[2 1 3]);

fprintf('\n');
