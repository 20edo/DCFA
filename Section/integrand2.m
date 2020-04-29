function I=integrand2(x,y,sc,result,rho,th)

if sc.pol==true
    x=rho.*cos(th);
    y=rho.*sin(th);
end
psi=interpolateSolution(result,x,y);
psi = reshape(psi,size(x));
I=psi.*y;