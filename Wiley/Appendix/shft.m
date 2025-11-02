function [out]=shft(A,n)
%This function shifts (circularly) the vector A with
% an amount of n
 
% Inputs:
% A  : the vector to be shifted
% n  : shift amount
 
% Output:
% out  : shifted vector 
 
out = A(l-n+2:l);
out(n:l) = A(1:l-n+1);
