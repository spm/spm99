function x = spm_mireg(varargin)
% Between modality coregistration using Mutual Information
% FORMAT x = spm_mireg(VG,VF,params)
% VG - handle for first image (see spm_vol).
% VF - handle for second image.
% x - the parameters describing the rigid body rotation.
%     such that a mapping from voxels in G to voxels in F
%     is attained by:  VF.mat\spm_matrix(x(:)')*VG.mat
% params - a cell array.
%          params{1} - optimisation sampling steps (mm)
%          params{2} - starting estimates (6 elements)
%
% The registration method used here is based on the work described in:
% A Collignon, F Maes, D Delaere, D Vandermeulen, P Suetens & G Marchal
% (1995) "Automated Multi-modality Image Registration Based On
% Information Theory". In the proceedings of Information Processing in
% Medical Imaging (1995).  Y. Bizais et al. (eds.).  Kluwer Academic
% Publishers.
%
% The original interpolation method described in this paper has been
% changed in order to give a smoother cost function.  The images are
% also smoothed slightly, as is the histogram.  This is all in order to
% make the cost function as smooth as possible, to give faster
% convergence and less chance of local minima.
%
% The mutual Information is essentially given by:
% H  = H/(sum(H(:)));
% s1 = sum(H,1);
% s2 = sum(H,2);
% H  = H.*log2((H+eps)./(s2*s1+eps));
% mi = sum(H(:));
%
% where H is a smoothed 256x256 joint histogram, and mi is the mutual information.
%
% This was subsequently refined to use the Entropy Correlation Coefficient:
% ecc = 2*mi/(-sum(s1.*log(s1))-sum(s2.*log(s2)));
%
% according to:
% F Maes, A Collignon, D Vandermeulen, G Marchal & P Suetens (1997).
% "Multimodality image registration by maximisation of mutual
% information". IEEE Transactions on Medical Imaging 16(2):187-198
%
% At the end, the voxel-to-voxel affine transformation matrix is
% displayed, along with the histograms for the images in the original
% orientations, and the final orientations.  The registered images are
% displayed at the bottom.
%_______________________________________________________________________
% @(#)spm_mireg.m	2.7 John Ashburner 01/09/28

if nargin>=4,
	x = optfun(varargin{:});
	return;
end;

if nargin < 1,
	VG = spm_vol(spm_get(1,'*.img','Select reference image'));
else,
	VG = varargin{1};
end;
if nargin < 2,
	VF = spm_vol(spm_get(1,'*.img','Select moved image'));
else,
	VF = varargin{2};
end;
if nargin <3, params = cell(0); else, params = varargin{3}; end;

if length(params)<1 | length(params{1}) < 1,
	params{1} = [4 2];
	if isglobal('sptl_MIStps'),
		global sptl_MIStps;
		params{1} = sptl_MIStps;
	end;
end;

if length(params)<2 | length(params{2}) > 12 | length(params{2}) == 0,
	params{2} = [0 0 0  0 0 0]';
end;
x  = params{2}(:);

% Voxel sizes (mm)
vxg   = sqrt(sum(VG.mat(1:3,1:3).^2));
vxf   = sqrt(sum(VF.mat(1:3,1:3).^2));

% Smoothnesses
fwhmg = sqrt(max([1 1 1]*params{1}(end)^2 - vxg.^2, [0 0 0]))./vxg;
fwhmf = sqrt(max([1 1 1]*params{1}(end)^2 - vxf.^2, [0 0 0]))./vxf;

if ~isfield(VG, 'uint8'),
	VG.uint8 = loaduint8(VG);
	VG       = smooth_uint8(VG,fwhmg); % Note side effects
end;
if ~isfield(VF, 'uint8'),
	VF.uint8 = loaduint8(VF);
	VF       = smooth_uint8(VF,fwhmf); % Note side effects
end;


sc = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001]'; % Required accuracy
sc = sc(1:length(params{2}));
xi = diag(sc*20);

for samp=params{1}(:)',
	[x,fval] = spm_powell(x(:), xi,sc,mfilename,VG,VF,samp);
	x        = x(:)';
end;
display_results(VG,VF,x);
return;
%_______________________________________________________________________

%_______________________________________________________________________
function o = optfun(x,VG,VF,s)
% The function that is minimised.  - (Normalised) Mutual Information
if nargin<4, s=[1 1 1]; end;

% Voxel sizes
vxg   = sqrt(sum(VG.mat(1:3,1:3).^2));sg = s./vxg;

% Create the joint histogram
H =   spm_hist2(VG.uint8,VF.uint8,    VF.mat\spm_matrix(x(:)')*VG.mat ,sg);

% Smooth the histogram
krn = exp(-([-6:6].^2)/8); krn = krn/sum(krn);
%krn = exp(-([-12:12].^2)/32); krn = krn/sum(krn);
H   = conv2(H,krn );
H   = conv2(H,krn');
%d   = 64;
%H   = sum(reshape(H,[256/d d 256]),1);
%H   = reshape(sum(reshape(H,[d 256/d d]),2),[d d]);

% Compute mutual information from histogram
H  = H+0.1;
H  = H/(sum(H(:)));
s1 = sum(H,1);
s2 = sum(H,2);
H  = H.*log2(H./(s2*s1));
mi = sum(H(:));

% actually use Entropy Correlation Coefficient of:
% Maes, Collignon, Vandermeulen, Marchal & Suetens (1997).
% "Multimodality image registration by maximisation of mutual
% information". IEEE Transactions on Medical Imaging 16(2):187-198
ecc = 2*mi/(-sum(s1.*log(s1))-sum(s2.*log(s2)));

% Negate, as the optimisation function minimises the cost function
o   = -ecc;
return;
%_______________________________________________________________________

%_______________________________________________________________________
function udat = loaduint8(V)
% Load data from file indicated by V into an array of unsigned bytes.
if size(V.pinfo,2)==1 & V.pinfo(1) == 2,
	mx = 255*V.pinfo(1) + V.pinfo(2);
	mn = V.pinfo(2);
else,
	spm_progress_bar('Init',V.dim(3),...
		['Computing max/min of ' spm_str_manip(V.fname,'t')],...
		'Planes complete');
	mx = -Inf; mn =  Inf;
	for p=1:V.dim(3),
		img = spm_slice_vol(V,spm_matrix([0 0 p]),V.dim(1:2),1);
		mx  = max([max(img(:))+paccuracy(V,p) mx]);
		mn  = min([min(img(:)) mn]);
		spm_progress_bar('Set',p);
	end;
end;
spm_progress_bar('Init',V.dim(3),...
	['Loading ' spm_str_manip(V.fname,'t')],...
	'Planes loaded');

udat = uint8(0);
udat(V.dim(1),V.dim(2),V.dim(3))=0;
rand('state',100);
for p=1:V.dim(3),
	img = spm_slice_vol(V,spm_matrix([0 0 p]),V.dim(1:2),1);
	acc = paccuracy(V,p);
	if acc==0,
		udat(:,:,p) = uint8(round((img-mn)*(255/(mx-mn))));
	else,
		% Add random numbers before rounding to reduce aliasing artifact
		r = rand(size(img))*acc;
		udat(:,:,p) = uint8(round((img+r-mn)*(255/(mx-mn))));
	end;
	spm_progress_bar('Set',p);
end;
spm_progress_bar('Clear');
return;

function acc = paccuracy(V,p)
if ~spm_type(V.dim(4),'intt'),
	acc = 0;
else,
	if size(V.pinfo,2)==1,
		acc = abs(V.pinfo(1,1));
	else,
		acc = abs(V.pinfo(1,p));
	end;
end;
%_______________________________________________________________________
%_______________________________________________________________________
function V = smooth_uint8(V,fwhm)
% Convolve the volume in memory (fwhm in voxels).
s  = fwhm/sqrt(8*log(2));
x  = round(6*s(1)); x = [-x:x];
y  = round(6*s(2)); y = [-y:y];
z  = round(6*s(3)); z = [-z:z];
x  = exp(-x.^2/(2*s(1).^2+eps));
y  = exp(-y.^2/(2*s(2).^2+eps));
z  = exp(-z.^2/(2*s(3).^2+eps));
x  = x/sum(x);
y  = y/sum(y);
z  = z/sum(z);

i  = (length(x) - 1)/2;
j  = (length(y) - 1)/2;
k  = (length(z) - 1)/2;
spm_conv_vol(V.uint8,V.uint8,x,y,z,-[i j k]);
return;
%_______________________________________________________________________

%_______________________________________________________________________
function display_results(VG,VF,x)
fig = spm_figure('FindWin','Graphics');
if isempty(fig), return; end;
set(0,'CurrentFigure',fig);
spm_figure('Clear','Graphics');

% Display text
%-----------------------------------------------------------------------
ax = axes('Position',[0.1 0.8 0.8 0.15],'Visible','off','Parent',fig);
text(0.5,0.7, 'Mutual Information Coregistration','FontSize',16,...
	'FontWeight','Bold','HorizontalAlignment','center','Parent',ax);

Q = inv(VF.mat\spm_matrix(x(:)')*VG.mat);
text(0,0.5, sprintf('X1 = %0.3f*X %+0.3f*Y %+0.3f*Z %+0.3f',Q(1,:)),'Parent',ax);
text(0,0.3, sprintf('Y1 = %0.3f*X %+0.3f*Y %+0.3f*Z %+0.3f',Q(2,:)),'Parent',ax);
text(0,0.1, sprintf('Z1 = %0.3f*X %+0.3f*Y %+0.3f*Z %+0.3f',Q(3,:)),'Parent',ax);

% Display scatter-plots
%-----------------------------------------------------------------------
ax  = axes('Position',[0.1 0.5 0.35 0.3],'Visible','off','Parent',fig);
H   = spm_hist2(VG.uint8,VF.uint8,VF.mat\VG.mat,[1 1 1]);
tmp = log(H+1);
image(tmp*(64/max(tmp(:))),'Parent',ax');
set(ax,'DataAspectRatio',[1 1 1],...
	'PlotBoxAspectRatioMode','auto','XDir','normal','YDir','normal',...
	'XTick',[],'YTick',[]);
title('Original Histogram','Parent',ax);
xlabel(spm_str_manip(VG.fname,'k22'),'Parent',ax);
ylabel(spm_str_manip(VF.fname,'k22'),'Parent',ax);

H   = spm_hist2(VG.uint8,VF.uint8,VF.mat\spm_matrix(x(:)')*VG.mat,[1 1 1]);
ax  = axes('Position',[0.6 0.5 0.35 0.3],'Visible','off','Parent',fig);
tmp = log(H+1);
image(tmp*(64/max(tmp(:))),'Parent',ax');
set(ax,'DataAspectRatio',[1 1 1],...
	'PlotBoxAspectRatioMode','auto','XDir','normal','YDir','normal',...
	'XTick',[],'YTick',[]);
title('Final Histogram','Parent',ax);
xlabel(spm_str_manip(VG.fname,'k22'),'Parent',ax);
ylabel(spm_str_manip(VF.fname,'k22'),'Parent',ax);

% Display ortho-views
%-----------------------------------------------------------------------
spm_orthviews('Reset');
h1 = spm_orthviews('Image',VG.fname,[0.01 0.01 .48 .49]);
h2 = spm_orthviews('Image',VF.fname,[.51 0.01 .48 .49]);
global st
st.vols{h2}.premul = inv(spm_matrix(x(:)'));
spm_orthviews('Space',h1);

return;

