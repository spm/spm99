function varargout = spm_clusters(varargin)
% Returns the cluster index for a point list - a compiled routine
% FORMAT [A] = spm_clusters(L)
% L     - locations [x y x]' {in voxels}
%
% A     - cluster index or region number
%_______________________________________________________________________
%
% spm_clusters characterizes a point list of voxel values defined with
% their locations (L) in terms of edge, face and vertex connected
% subsets, returning a list of indices in A, such that the ith location
% belongs to cluster A(i) (using an 18 connectivity scheme).
%
%_______________________________________________________________________
% @(#)spm_clusters.m	2.2 FIL 99/04/19

%-This is merely the help file for the compiled routine
error('spm_clusters.c not compiled - see spm_MAKE.sh')
