%example script of use of ineria for Toeplitz block matrix
%Generation of a square Toeplitz Block matrix A, with 3*3 blocks. In this
%case we generate a symmetric matrix.
v1=rand(1,4);
T11=toeplitz(v1);
v2=rand(1,4);
T12=toeplitz (v2);
T22=toeplitz([3:6]);
T13=toeplitz([-4:-1]);
T33=toeplitz([20:23]);
v3=randn(1,4);
T23=toeplitz(v3);
A=[T11,T12, T13;T12',T22, T23; T13',T23',T33]

%Computation of generator matrices. The second argument is the number of
%blocks (the same for row blocks than for column blocks).

[Gn,Hn] = generadores_toep_block(A,3);

% The mex inertia_toep_mex returns the same vector than algorithm, i.e.,
% the number of negative eigenvalues of A.
v_inercia=inercia_toep_mex(A,Gn,Hn)

%The mex inertia_toep_mex2 returns 1 if all the eigenvalues of A are
%positive, and 0 if there is any eigenvalue negative or cero. This mex is
%at least as fast as inercia_toep_mex, and very often it is substantially
%faster.
estable_inerciam=inercia_toep_mex2(A,Gn,Hn)

%check the results using matlab's function eig.
eig(A)

