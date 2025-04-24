function fun  = funSThM(mu,beta)
	if iscell(mu)
		mu = mu{1};
	end
% Set the number of constraints
np = 3; %three calibration parameters
nz = length(beta)-np;
nx = int32((length(mu)-nz)/2);

q = nx+nz;
fun = zeros(q,1);

%Calibration samples
for i=1:nx
  fun(i) = beta(1)*mu(i)/(beta(2)+ mu(i))+beta(3) - mu(i+nx);
end
%Unknown samples
for i=1:nz
    fun(nx+i) = beta(1)*beta(i+np)/(beta(2)+ beta(i+np))+beta(3) - mu(i+2*nx);
end
end

