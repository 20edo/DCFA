%%
% This script builds the typical section matrices for the unsteady
% aerodynamics.

% DCFA swept wing assignement
%
% Teamwork
% Team members: Venti Edoardo         944421
%               Zemello Matteo        942003
%               Zucchelli Umberto     952952
%
%
%
clear all,close all, clc

%Add to path the necessary folders
cd ..
cd ..
cd ..
init
cd Aircraft\aero_analysis\Unsteady

%% Generate model of a section clamped for all gdl except rotation around x and torsion

chord = 9;        % Chord
e_ca = 1.6;          % e -> distance between the shear center (elastic axis) and AC
m = 1.95;         % mass of the element
b=chord/2;
span = 1;         % span of the half wing

K = [14536517779.8877, 0; 0, 2032828194795.00];
M = [130.972022554927, -16.0315875316925; -16.0315875316925, 318.094945246663];
Cs = [0,0; 0,0];

[X_old, e_old] = polyeig(K,Cs,M);


[T,a,P,rho] = atmosisa(10000);
v = [0:1:200];
q = 1/2*rho.*v.^2;

eig_ = zeros(length(v),size(X_old,2));
eig_(1,:) = e_old;


%Following iterations
for i=2:length(v)
    for index=1:length(e_old)
        tic
        e_old(k)
        Build matrices to solve the non-linear system
        k = b*imag(e_old(index))/v(i);
        
        Q_th = [[                                                                                                                                                                                                                                                                                                                                                                                                                                                         - (3*pi*k^2)/(4*m) + (pi*k*2i)/m,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           (2*((3*3^(1/2))/4 - (1014840697911479*k^2)/9007199254740992))/m + (pi*k*2i)/(3*m),                                                                                                                                                                                                                                                -(pi*k^2)/(b*m)]
            [ (((39782776467233*imag(besselh(1, 2, k)/(besselh(1, 2, k) + besselh(0, 2, k)*1i)))/562949953421312 + (39782776467233*k*real(besselh(1, 2, k)/(besselh(1, 2, k) + besselh(0, 2, k)*1i)))/562949953421312 + (3^(1/2)*k)/8)*2i)/m - (2*((1014840697911479*k^2)/9007199254740992 - (39782776467233*real(besselh(1, 2, k)/(besselh(1, 2, k) + besselh(0, 2, k)*1i)))/562949953421312 + (39782776467233*k*imag(besselh(1, 2, k)/(besselh(1, 2, k) + besselh(0, 2, k)*1i)))/562949953421312))/m, (2*((1217813138399833*real(besselh(1, 2, k)/(besselh(1, 2, k) + besselh(0, 2, k)*1i)))/9007199254740992 - (479206101361749*k^2)/9007199254740992 - (1653738962577673*k*imag(besselh(1, 2, k)/(besselh(1, 2, k) + besselh(0, 2, k)*1i)))/36028797018963968 + 529960012437871/2251799813685248))/(m*pi) + (((7186389630543293*k)/18014398509481984 + (1217813138399833*imag(besselh(1, 2, k)/(besselh(1, 2, k) + besselh(0, 2, k)*1i)))/9007199254740992 + (1653738962577673*k*real(besselh(1, 2, k)/(besselh(1, 2, k) + besselh(0, 2, k)*1i)))/36028797018963968)*2i)/(m*pi), - (2*((567094513656589*k^2)/4503599627370496 + (39782776467233*k*imag(besselh(1, 2, k)/(besselh(1, 2, k) + besselh(0, 2, k)*1i)))/562949953421312))/(b*m) + (k*real(besselh(1, 2, k)/(besselh(1, 2, k) + besselh(0, 2, k)*1i))*39782776467233i)/(281474976710656*b*m)]
            [                                                                                                                                                                   - (2*pi*(k^2/2 - 2*real(besselh(1, 2, k)/(besselh(1, 2, k) + besselh(0, 2, k)*1i)) + 2*k*imag(besselh(1, 2, k)/(besselh(1, 2, k) + besselh(0, 2, k)*1i))))/m + (pi*(k + 2*imag(besselh(1, 2, k)/(besselh(1, 2, k) + besselh(0, 2, k)*1i)) + 2*k*real(besselh(1, 2, k)/(besselh(1, 2, k) + besselh(0, 2, k)*1i)))*2i)/m,                                                                                           - (2*((567094513656589*k^2)/4503599627370496 - (8616390187129275*real(besselh(1, 2, k)/(besselh(1, 2, k) + besselh(0, 2, k)*1i)))/2251799813685248 + (3*3^(1/2)*k*imag(besselh(1, 2, k)/(besselh(1, 2, k) + besselh(0, 2, k)*1i)))/4))/m + (((5532085316927607*k)/9007199254740992 + (8616390187129275*imag(besselh(1, 2, k)/(besselh(1, 2, k) + besselh(0, 2, k)*1i)))/2251799813685248 + (3*3^(1/2)*k*real(besselh(1, 2, k)/(besselh(1, 2, k) + besselh(0, 2, k)*1i)))/4)*2i)/m,                                                                                          - (2*pi*(k^2 + 2*k*imag(besselh(1, 2, k)/(besselh(1, 2, k) + besselh(0, 2, k)*1i))))/(b*m) + (k*pi*real(besselh(1, 2, k)/(besselh(1, 2, k) + besselh(0, 2, k)*1i))*4i)/(b*m)]];
        Q_th = Q_th([1,3],[1,3]); 
        C = [chord^2, 0; 0, chord];
        J = [1, 0; -e_ca, -1];
        Ham = J'*C*Q_th*J;
        
        [X,e] = polyeig(K-q(i)*Ham,Cs,M);
        e
        d_e = abs(e-e_old(index));
        d_e
        if 1 %sum(d_e < 0.1) == 1
            e = e(find(d_e == min(d_e)));
            X = X(:,find(d_e == min(d_e)));
        else
            guess = find(d_e < 0.1);
            d_v = vecnorm(X_old(:,guess)-X);
            index = guess(find(d_v == min(d_v)));
            e = e(index);
            X = X(:,index);
        end
        X_old(:,index) = X;
        e_old(index) = e;
        phrase = ['Eig number ',num2str(index),' out of ',num2str(length(e_old)),'; Velocity ',num2str(i),' out of ',num2str(length(v))];
        disp(phrase)
        toc
        
    end
    eig_(i,:) = e_old;
end