

%---------------------------------------------------------------
% user variables defined here 
%---------------------------------------------------------------


p = '/home/shallice/Rik/SPM99course/CatStats/params.mat'
load(p);
F1 = spm_get_files('/home/shallice/Rik/SPM99course/Smooth_snraf','sn*.img');

z=zeros(1,8);
c1 = [1 0 1 0 1 0 1 0 0 0   0];
c2 = -c1;

cond_names = {'n1','n2','f1','f2','target'}

o1 = [ 4.7451    6.7451   10.9902   13.2353   14.5000   20.5000   26.7451];
o2 = [28.7451   37.2353   46.7451   54.9902   56.9902   59.2353   68.7451];



%---------------------------------------------------------------
% batch variables defined here for analysis 'model'
%---------------------------------------------------------------


model(1) = struct( ...
 'types',          1, ...
 'global_effects', {'Scaling'}, ...
 'burst_mode',     0, ...
 'HF_fil',         'specify',  ...
 'HF_cut',         [232 232], ...
 'LF_fil',         'Gaussian', ...
 'LF_cut',         4, ...
 'int_corr',       'none', ... 
 'now_later',      1, ...
 'stop_writing',   0, ...
 'trial_fcon',     0, ...
 'RT',             4.1, ...
 'replicated',     0, ...
 'nsess',          2, ...
 'nscans',         [100 200], ...
 'files',          {{F1 F1}}, ...
 'conditions_nb',  [5 2], ...     
 'conditions',     [1 2], ...
 'regressors_nb',  [0 0], ...
 'regressors',     [], ...
 'parametrics_type', {{'none','none'}}, ...
 'parametrics',    [] ...
 'stochastics_flag', [0 0], ...
 'stochastics',    [], ...
);

model(2) =  model(1);
model(2).HF_cut = [232];
model(2).nscans = [100];
model(2).files = {{F1}};
model(2).conditions_nb = [5];     
model(2).conditions = [1];
model(2).regressors_nb = [0];
model(2).parametrics_type = {{'none'}};
model(2).stochastics_flag = [0];
model(2).stochastics = [];

%-------------------------------------------
 
conditions(1) = struct( ...
 'names',   {cond_names}, ...
 'onsets',  {sot}, ...   
 'types',    {{'events','events','events','events','events'}}, ... 
 'bf_ev',   [1 1 1 1 1], ...
 'bf_ep',   [0 0 0 0 0], ...
 'volterra',  0, ...
 'variable_dur', 0 ...
);

conditions(2) = conditions(1);
conditions(2).names = {'n1','n2'};
conditions(2).onsets = {o1, o2};
conditions(2).types = {{'events','events'}};
conditions(2).bf_ev = [1 0];
conditions(2).bf_ep = [0 1];

%-------------------------------------------

bf_ev(1) = struct( ...
  'ev_type', 2, ...
  'win_len', [] ...
);
bf_ev(2) = bf_ev(1) ;
bf_ev(2).ev_type = 1;

bf_ep(1) = struct( ...
  'ep_type', 4, ...
  'length',3, ...
  'conv',  1, ...  
  'deriv', 0 ...
);

%---------------------------------------------------------------
% batch variables defined here for analysis 'contrasts'
%---------------------------------------------------------------

contrasts(1) = struct( ...
   'names',  {{'c1','c2'}}, ...
   'types',   {{'T','F'}}, ...
   'values',  {{c1,c2}} ...
);
