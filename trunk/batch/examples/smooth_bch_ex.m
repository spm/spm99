
%---------------------------------------------------------------
% user variables defined here 
%---------------------------------------------------------------

% enter the file names in a variable (called F here) with 
% F = ['/home/spm/spm99b/jba_devel/test/beta_0001.img'];
% or with
% F = spm_get_files('DIRECTORY','FILTER')
% as in 
% F = spm_get_files('/home/shallice/Rik/SPM99course/Smooth_snraf','sn*.img');

%---------------------------------------------------------------
% batch variables defined here for analysis 'smooth'
%---------------------------------------------------------------

smooth = struct(...
'FWHMmm', 4, ...
'files', F ...
);

