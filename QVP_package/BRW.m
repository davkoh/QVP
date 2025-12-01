function [output]=BRW(y,x,taus)
%Coded by Tibor Szendrei following Bondell, Reich, Wang (2010). 


[n,p]=size(x);
m=max(size(taus));
pp=p+1; %for constant
Y=repmat(y,m,1); %repeat y for each quantiles
xtemp=x;

%SCALING
x=zeros(n,p);
shifts=min(xtemp);
x=xtemp-repmat(shifts,n,1);
scalings=max(x);
x=x./repmat(scalings,n,1);

%Save the scaling and shift for later: will need to do the inverse xform on
%the betas (Fine to do since QR is nonlinear transformation equivariant)

%MAKE X MATRIX FOR JOINT ESTIMATION
x=[ones(n,1),x];
D=ones(m);
D=tril(D);
X=kron(D,x);

X2=[X,-X]; %This is needed for NC constraints. Esseintally we restrict parameters to be above 0 so to get negative params we need the same variables but in negative
%Converting X2 to sparse matrix helps in computation speed

K=m*pp;

%RESTRICTIONS
%Restrictions will be Rb>=r
%We want constant to be larger than the sum of negative differences. R1 is
%for positive constants (this should be positive), R2 is for negative coefficients

%Adaptive NC constraint
kronmat=eye(m);

negshift=ones(1,p); %Negative shifters are the maximum of x (which is 1 after scaling)
posshift=zeros(1,p); %Positive shifters are minimum of x (which is 0 after scaling)

R1temp=[1,posshift];
R1=kron(kronmat,R1temp);

R2temp=[1,negshift];
R2=kron(kronmat,R2temp);

R1=R1(2:end,:); %We drop the constraint for first quantile: estimate without restriction
R2=R2(2:end,:); %We drop the constraint for first quantile: estimate without restriction

NC=[R1,-R2];

%All constraints together
R=[eye(2*K);NC]; %Add identity matrix to the top
%Converting R to sparse matrix helps in computation speed

%r
r=zeros(2*K+(m-1),1);

%b
x1=reshape(ones(n,1)*taus,[],1);
b=X2'*(1-x1);

%Other parameters needed for QR
u=ones(size(x1,1),1);
x2=ones(size(R,1),1);

%ESTIMATION
%Uses interior point method
coeff1=-lp_fnm(X2',-Y',R',-r',b,u,x1,x2); %Functions from Roger Koenker's rqic codefile
%Calculate coefficients using solution of interior method
coeff=coeff1(1:K)-coeff1((K+1):end);
gamma=reshape(coeff,[],m);
D=ones(m);
D=triu(D);
bhat_temp=gamma*D;

%Calculate transformation matrix
transform=[1,-shifts./scalings];
transform=[transform;zeros(p,1),diag(1./scalings)];

output=transform*bhat_temp; %final betas

end