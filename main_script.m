clear all
close all
clc

%Wrapper function to send each flux/variance calculation to a seperate node

%First grabs info about the input file, sets up dimensional info, and 
%then loops through variables to send individual batch jobs 

%Manages the jobs and checks for job failures, resubmitting as needed

%Clears the temp folder and the model output after completion of the script

%Variables to process, vertical fluxes and variances
fluxvars={'U';'T';'CO';'O';'CO2';'H2O';'O3';'NO';'HOX'};
varvars={'T';'OMEGA'};
flux_var_switchover=length(fluxvars);
for i=1:length(fluxvars)
   processvars{i}=fluxvars{i};
end
for i=1:length(varvars)
   processvars{i+flux_var_switchover}=varvars{i};
end

%Get basic information needed to run
experiment=getenv('POSTPROC_RUN');
days=str2num(getenv('DAYS'));

%Grab the most recent file
maindir=['/glade/scratch/nadavis/archive/',experiment,'/atm/hist/'];
maindir
filelist=getfilenames(maindir,'*cam.h2*00000.nc');
file=filelist{1};
file_structure=get_file_structure(file,days);

%Variables for job management
count=0;                         %Currently running jobs
outer_lim=-1;                    %While loop constraint; bootstrap to allow loop to begin
varsleft=length(processvars);    %How many jobs left to do
varstotal=varsleft;              %How many jobs in total
varnum_current=[];               %Current job list
count_limit=4;                   %Maximum number of jobs

%Main loop
while count>outer_lim

	%Monitor each job, resubmit failures/submit new job
	while length(varnum_current) == count_limit && count>0
           vars_to_remove=[];

	   for c=1:length(varnum_current)
	      varnum=varnum_current(c);
              state=jobs{varnum}.State;

	      if strcmp(state,'finished')
	      
	         %Move on to new variable (in main loop) if successfully processed
	         count=count-1;
		 vars_to_remove=cat(1,vars_to_remove,varnum);

	      elseif strcmp(state,'failed')
	      
		 jobs{varnum}
		 %Resubmit the same variable if the job failed
		 [jobs{varnum},files{varnum}]=submit_batch_job(varnum,file_structure,file,experiment,maindir,processvars{varnum},flux_var_switchover);
 	      end
	   end

           %Remove variables if successfully processed
	   if ~isempty(vars_to_remove)
	      for i=1:length(vars_to_remove)
	         varnum_current=remove_var_from_jobs(varnum_current,vars_to_remove(i));
	      end
	   end

	   %Adjust maximum number of jobs once there are no new vars to process
	   if varsleft==0
	      count_limit=length(varnum_current);
	   end
	end

	%Only submit a new variable if theres a var left to process 
	if varsleft>0

	   %Unique identifier for each job
	   varnum=get_var_num(varsleft,varstotal);
		
           [jobs{varnum},files{varnum}]=submit_batch_job(varnum,file_structure,file,experiment,maindir,processvars{varnum},flux_var_switchover);
	   varnum_current=add_var_to_jobs(varnum_current,varnum);
	
	   count=count+1;
	   varsleft=varsleft-1;
 
	end

   	outer_lim=0;
	count=length(varnum_current);
end

%Clear out the temporary directory
cd('output')
delete('*.mat') 
delete('*.txt')
rmdir('Job*','s')

%Delete files
for f=1:length(files)
   delete(files{f})
end

exit;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function to generalize submitting a job to process a given variable
function [j,input_file]=submit_batch_job(varnum,file_structure,file,experiment,maindir,var,flux_var_switchover)

%Set up worker information
worker_input=struct;
worker_input.experiment=experiment;
if varnum>flux_var_switchover
   worker_input.copy_var=1;
   worker_input.var1=var;
   worker_input.var2=var;
   worker_input.file=[strrep(file,'.nc',''),'.',var,'.nc'];
else
   worker_input.copy_var=0;
   worker_input.var2=var;
   worker_input.var1='OMEGA';
   worker_input.file=[strrep(file,'.nc',''),'.',var,'OMEGA.nc'];
end
worker_input.file

% Start PBS cluster
% Add cluster profile if not already present
if ~any(strcmp(parallel.clusterProfiles, 'ncar_mps'))
    ncar_mps = parallel.importProfile('/glade/u/apps/opt/matlab/parallel/ncar_mps.mlsettings');
end

% Start PBS cluster and submit job with custom number of workers
c = parcluster('ncar_mps');

% Setup cluster attributes - 1 node (35 workers) per cluster
jNodes = 1;
jTasks = getenv('MPSTASKS');
jWorkers = jNodes * str2num(jTasks)-1;
jAccount = getenv('MPSACCOUNT');
jQueue = getenv('MPSQUEUE');
jWalltime = getenv('MPSWALLTIME');
%c.NumWorkers=str2num(jTasks);

% Initiate job
c.ClusterMatlabRoot = getenv('NCAR_ROOT_MATLAB');
c.ResourceTemplate = append('-l select=',num2str(jNodes),':ncpus=',num2str(jWorkers),':mpiprocs=', num2str(jWorkers),':mem=109GB');
c.SubmitArguments = append('-A ', jAccount, ' -q ', jQueue, ' -l walltime=', jWalltime);
c.JobStorageLocation = append(getenv('PWD'), '/output');

% Output cluster settings
c

% Submit batch job, ensure worker has access to file so it doesn't have to create a local copy
j = batch(c, @process_flux_spectra_worker, 0, {file_structure,worker_input}, 'pool', jWorkers, 'AdditionalPaths', {maindir;'/glade/scratch/nadavis/nsc/'});

input_file=worker_input.file

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function to get the next variable number 
function varnum=get_var_num(varsleft,varstotal)
	varnum=varstotal-varsleft+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Add varnum to currently running jobs
function varnum_current=add_var_to_jobs(varnum_current,varnum)
	varnum_current=cat(1,varnum_current,varnum);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Remove varnum from currently running jobs
function varnum_current=remove_var_from_jobs(varnum_current,varnum)
	ind=find(varnum_current==varnum,1,'first');
	varnum_current(ind)=[];
end
