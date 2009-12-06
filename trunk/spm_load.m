function [x] = spm_load(f)
% function to load ascii file data as matrix
% FORMAT [x] = spm_load(f)
% f  - file {ascii file containing a regular array of numbers
% x  - corresponding data matrix
%_______________________________________________________________________
% 98/11/17 Karl Friston, Andrew Holmes @(#)spm_load.m	2.1


%-Get a filename if none was passed
%-----------------------------------------------------------------------
if nargin==0
	f = spm_get(1,'*','Select ASCII data file');
end


%-Load the data file into double precision matrix x
%-----------------------------------------------------------------------
x = load(f,'-ascii');
