function [m] = ncc(win1,win2)
%NCC Summary of this function goes here
%   Detailed explanation goes here
    s = size(win1,1);
    nc = normxcorr2(win1,win2);
    m = nc(s,s);
end

