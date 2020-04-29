function I=integrand(x,y,sc,result,rho,th)

if sc.pol==true
    x=rho.*cos(th);
    y=rho.*sin(th);
end
[gradx,grady] = evaluateGradient(result,x,y);
gradx=reshape(gradx,size(x));
grady=reshape(grady,size(y));
I=sc.G(x,y).*((grady+x).*x-(gradx-y).*y);
end
