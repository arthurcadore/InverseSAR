function out=cevir2(a,nx,ny)
% This function converts a 1D vector to a 2D matrix
 
% Inputs:
% a   : 1D vector of length (nx*ny)
% nx  : column length of the output matrix
% ny  : row length of the output matrix
 
% Output:
% out  : 2D matrix of size nx by ny
 
 
for p = 1:nx;
      out(p,1:ny) = a(1,(p-1)*ny+1:ny*p);
end;
