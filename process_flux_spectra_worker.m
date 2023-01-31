
function process_flux_spectra_worker(file_structure,worker_input)

year=file_structure.date(1)/10000;
month=(file_structure.date(1)-year*10000)/100;
day=(file_structure.date(1)-year*10000-month*100);

%Set up nc file and then dump each constituent as we go
filestring=[worker_input.experiment,'_',sprintf('%0.4i',year),'-',sprintf('%0.2i',month),'-',sprintf('%0.2i',day),'_',...
  worker_input.var1,worker_input.var2,'_space_time_spectra.nc'];

if exist(filestring,'file')==2
   delete(filestring);
end
ncid=netcdf.create(filestring,'NETCDF4');
titleid=netcdf.getConstant('GLOBAL');
authid=netcdf.getConstant('GLOBAL');
dateid=netcdf.getConstant('GLOBAL');

npres=netcdf.defDim(ncid,'lev',file_structure.nlev);
nreg=netcdf.defDim(ncid,'region',2);
nnlat=netcdf.defDim(ncid,'lat',file_structure.numlat);
ncspeed=netcdf.defDim(ncid,'phasespeed',file_structure.nc);
nwavenum=netcdf.defDim(ncid,'wavenum',file_structure.nwave);
nwavenum_local=netcdf.defDim(ncid,'wavenum_local',file_structure.nwave_expansion);

vpres=netcdf.defVar(ncid,'lev','double',npres);
vcspeed=netcdf.defVar(ncid,'phasespeed','double',ncspeed);
vwavenum=netcdf.defVar(ncid,'wavenum','double',nwavenum);
vwavenum_local=netcdf.defVar(ncid,'wavenum_local','double',nwavenum_local);
vlat=netcdf.defVar(ncid,'lat','double',nnlat);
var=netcdf.defVar(ncid,[worker_input.var1,worker_input.var2],'double',[npres ncspeed nwavenum nnlat]);
var_local=netcdf.defVar(ncid,[worker_input.var1,worker_input.var2,'_local'],'double',[npres ncspeed nwavenum_local nnlat]);
netcdf.endDef(ncid);

%Parallel collectors - global and local wavenumbers
flux_spectra_collector_full=zeros(file_structure.nlev,file_structure.nc,file_structure.nwave,file_structure.numlat);
flux_spectra_collector_local=zeros(file_structure.nlev,file_structure.nc,file_structure.nwave_expansion,file_structure.numlat);

for l=1:file_structure.nlev
 clc
 l

 in_var_1=permute(squeeze(ncread(worker_input.file,worker_input.var1,[1 1 file_structure.lev_ind(l) 1],[Inf Inf 1 Inf])),[3 2 1]);
 if (worker_input.copy_var==1)
    in_var_2=in_var_1;
 else
    in_var_2=permute(squeeze(ncread(worker_input.file,worker_input.var2,[1 1 file_structure.lev_ind(l) 1],[Inf Inf 1 Inf])),[3 2 1]);
 end

 %Each spectrum
 for reg=1:2
     % Full zonal or limited area
     if reg==1
        lon_range=1:file_structure.nlon;
     else
        lon_range=file_structure.wlon:file_structure.elon;
     end

     if reg==1
        [cspeed_out,wavenum_out,latout,flux_spectra_collector_full(l,:,:,:)] = ...
              randel_held_spectra(in_var_1(:,file_structure.slat:file_structure.nlat,lon_range),...
              in_var_2(:,file_structure.slat:file_structure.nlat,lon_range),...
              file_structure.lat(file_structure.slat:file_structure.nlat),[1 size(in_var_1,1)],file_structure.windowing{reg},file_structure.L,...
              length(lon_range)/file_structure.nlon);
     else 
        [cspeed_out,wavenum_out,latout,flux_spectra_collector_local(l,:,:,:)] = ...
              randel_held_spectra(in_var_1(:,file_structure.slat:file_structure.nlat,lon_range),...
              in_var_2(:,file_structure.slat:file_structure.nlat,lon_range),...
              file_structure.lat(file_structure.slat:file_structure.nlat),[1 size(in_var_1,1)],file_structure.windowing{reg},file_structure.L,...
              length(lon_range)/file_structure.nlon);
     end

     %Save dimension info
     %Doesn't matter from which core
     if reg==1
        wavenum=wavenum_out;
        cspeed=cspeed_out;
     else
        wavenum_local=wavenum_out;
     end
  end
end

disp('Saving dims')
%Save dimensions
netcdf.putVar(ncid,vlat,file_structure.lat(file_structure.slat:file_structure.nlat));
netcdf.putVar(ncid,vpres,file_structure.lev);
netcdf.putVar(ncid,vcspeed,cspeed);
netcdf.putVar(ncid,vwavenum,wavenum);
netcdf.putVar(ncid,vwavenum_local,wavenum_local);

%Dump spectra to file, close
disp('Writing to file')
netcdf.putVar(ncid,var,flux_spectra_collector_full);
netcdf.putVar(ncid,var_local,flux_spectra_collector_local);
netcdf.close(ncid);
