% This function rotates the mass matrix and the stiffness matrix of the
% assembled beam in order to assemble different beams with different
% orientations
%%
% DCFA swept wing assignement
%
% Teamwork
% Team members: Venti Edoardo         944421
%               Zemello Matteo        942003
%               Zucchelli Umberto     952952
%               
%           
%
function b = b_rotate(b)
% global reference 
X1 = [1,0,0]'; 
X2 = [0,1,0]'; 
X3 = [0,0,1]'; 
% local reference b.vx b.vy
X1Prime = b.vx; 
X2Prime = b.vy; 
X3Prime = cross(b.vx,b.vy);

% assembly of the rotation matrix 3x3
A = zeros(3);   
A(1,1) = dot(X1, X1Prime);
A(1,2) = dot(X1, X2Prime);
A(1,3) = dot(X1, X3Prime);
A(2,1) = dot(X2, X1Prime);
A(2,2) = dot(X2, X2Prime);
A(2,3) = dot(X2, X3Prime);
A(3,1) = dot(X3, X1Prime);
A(3,2) = dot(X3, X2Prime);
A(3,3) = dot(X3, X3Prime);

% assembly of the big rotation matrix
N = size(b.M,1)/3;
Ar = repmat(A, 1, N);
Ac = mat2cell(Ar, size(A,1), repmat(size(A,2),1,N));
R = sparse(blkdiag(Ac{:}));                                   

% rotate matrices
b.M = R*b.M*R';
b.K = R*b.K*R';
b.navier=b.navier*R';