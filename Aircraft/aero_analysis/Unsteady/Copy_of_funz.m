function res=Copy_of_funz(M,Cs,K,v,q,X,e)

%% Compute loads
chord=7.72;
k=chord/2*imag(e)/norm(v);

%% Reduce matrices
chord = 9;        % Chord
e_ca = 1.6;          % e -> distance between the shear center (elastic axis) and AC
m = 1.95;         % mass of the element
b=chord/2;
span = 1;  

Q_th = [[                                                                                                                                                                                                                                                                                                                                                                                                                                                         - (3*pi*k^2)/(4*m) + (pi*k*2i)/m,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           (2*((3*3^(1/2))/4 - (1014840697911479*k^2)/9007199254740992))/m + (pi*k*2i)/(3*m),                                                                                                                                                                                                                                                -(pi*k^2)/(b*m)]
[ (((39782776467233*imag(besselh(1, 2, k)/(besselh(1, 2, k) + besselh(0, 2, k)*1i)))/562949953421312 + (39782776467233*k*real(besselh(1, 2, k)/(besselh(1, 2, k) + besselh(0, 2, k)*1i)))/562949953421312 + (3^(1/2)*k)/8)*2i)/m - (2*((1014840697911479*k^2)/9007199254740992 - (39782776467233*real(besselh(1, 2, k)/(besselh(1, 2, k) + besselh(0, 2, k)*1i)))/562949953421312 + (39782776467233*k*imag(besselh(1, 2, k)/(besselh(1, 2, k) + besselh(0, 2, k)*1i)))/562949953421312))/m, (2*((1217813138399833*real(besselh(1, 2, k)/(besselh(1, 2, k) + besselh(0, 2, k)*1i)))/9007199254740992 - (479206101361749*k^2)/9007199254740992 - (1653738962577673*k*imag(besselh(1, 2, k)/(besselh(1, 2, k) + besselh(0, 2, k)*1i)))/36028797018963968 + 529960012437871/2251799813685248))/(m*pi) + (((7186389630543293*k)/18014398509481984 + (1217813138399833*imag(besselh(1, 2, k)/(besselh(1, 2, k) + besselh(0, 2, k)*1i)))/9007199254740992 + (1653738962577673*k*real(besselh(1, 2, k)/(besselh(1, 2, k) + besselh(0, 2, k)*1i)))/36028797018963968)*2i)/(m*pi), - (2*((567094513656589*k^2)/4503599627370496 + (39782776467233*k*imag(besselh(1, 2, k)/(besselh(1, 2, k) + besselh(0, 2, k)*1i)))/562949953421312))/(b*m) + (k*real(besselh(1, 2, k)/(besselh(1, 2, k) + besselh(0, 2, k)*1i))*39782776467233i)/(281474976710656*b*m)]
[                                                                                                                                                                   - (2*pi*(k^2/2 - 2*real(besselh(1, 2, k)/(besselh(1, 2, k) + besselh(0, 2, k)*1i)) + 2*k*imag(besselh(1, 2, k)/(besselh(1, 2, k) + besselh(0, 2, k)*1i))))/m + (pi*(k + 2*imag(besselh(1, 2, k)/(besselh(1, 2, k) + besselh(0, 2, k)*1i)) + 2*k*real(besselh(1, 2, k)/(besselh(1, 2, k) + besselh(0, 2, k)*1i)))*2i)/m,                                                                                           - (2*((567094513656589*k^2)/4503599627370496 - (8616390187129275*real(besselh(1, 2, k)/(besselh(1, 2, k) + besselh(0, 2, k)*1i)))/2251799813685248 + (3*3^(1/2)*k*imag(besselh(1, 2, k)/(besselh(1, 2, k) + besselh(0, 2, k)*1i)))/4))/m + (((5532085316927607*k)/9007199254740992 + (8616390187129275*imag(besselh(1, 2, k)/(besselh(1, 2, k) + besselh(0, 2, k)*1i)))/2251799813685248 + (3*3^(1/2)*k*real(besselh(1, 2, k)/(besselh(1, 2, k) + besselh(0, 2, k)*1i)))/4)*2i)/m,                                                                                          - (2*pi*(k^2 + 2*k*imag(besselh(1, 2, k)/(besselh(1, 2, k) + besselh(0, 2, k)*1i))))/(b*m) + (k*pi*real(besselh(1, 2, k)/(besselh(1, 2, k) + besselh(0, 2, k)*1i))*4i)/(b*m)]]; 
Q_th = Q_th([1,3],[1,3]);
Q_th(2,1) = -Q_th(2,1); 
Q_th(1,2) = -Q_th(1,2); 
C = [chord^2, 0; 0, chord]; 
J = [1, 0; -e, -1];
Ham = J'*C*Q_th*J; 
%% Calculate residual
res=(e^2*M+e*Cs+K-q*Ham)*X;

end