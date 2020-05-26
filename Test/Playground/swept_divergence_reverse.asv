% Problem of divergence and control reversal for a swept wing
% given a model containing only the wing



function element = el_staticaero_assebl(element,right_wing)

Lambda x GJ EJ Le e c CL0 CLa CMac CLb CMb real
GJ = element.sc.GJ; 
EJ = element.sc.EJy; 
e = element.sc.e % TROVARE UN MODO PER ESPRIMERE e 
c = element.sc.c % CORDA
CL0 = element.sc.CL0; % CL0, cl a zero incidenza 
CLa = element.sc.CLa; % CL derivato alpha (2pi???) 
CLb = element.sc. CLb; % CM derivato alpha 
CMb = element.sc. CMb; % comefficiente di momento derivato beta 

% K_thth = [[  GJ/Le, -GJ/Le];
%         [ -GJ/Le,  GJ/Le]];
%     
% K_ww = [[  (12*EJ)/Le^3, -(6*EJ)/Le^2, -(12*EJ)/Le^3, -(6*EJ)/Le^2];
%         [  -(6*EJ)/Le^2,    (4*EJ)/Le,   (6*EJ)/Le^2,    (2*EJ)/Le];
%         [ -(12*EJ)/Le^3,  (6*EJ)/Le^2,  (12*EJ)/Le^3,  (6*EJ)/Le^2];
%         [  -(6*EJ)/Le^2,    (2*EJ)/Le,   (6*EJ)/Le^2,    (4*EJ)/Le]];

%% Assembly of the aerodynamics matrices
fa_th = [(Le*cos(Lambda)^2*(CMac*c^2 + CL0*e*c))/2; 
    (Le*cos(Lambda)^2*(CMac*c^2 + CL0*e*c))/2]; 

fa_w = [- (sin(2*Lambda)*(CMac*c^2 + CL0*e*c))/2 - (CL0*Le*c*cos(Lambda))/2;
                                         (CL0*Le^2*c*cos(Lambda))/12;
   (sin(2*Lambda)*(CMac*c^2 + CL0*e*c))/2 - (CL0*Le*c*cos(Lambda))/2;
                                        -(CL0*Le^2*c*cos(Lambda))/12];

Ka_thth = [[ (CLa*Le*c*e*cos(Lambda)^3)/3, (CLa*Le*c*e*cos(Lambda)^3)/6];
          [ (CLa*Le*c*e*cos(Lambda)^3)/6, (CLa*Le*c*e*cos(Lambda)^3)/3]];

Ka_thw = [[ (CLa*c*e*sin(Lambda)*(sin(Lambda)^2 - 1))/2,  (CLa*Le*c*e*sin(Lambda)*(sin(Lambda)^2 - 1))/12, -(CLa*c*e*sin(Lambda)*(sin(Lambda)^2 - 1))/2, -(CLa*Le*c*e*sin(Lambda)*(sin(Lambda)^2 - 1))/12];
[ (CLa*c*e*sin(Lambda)*(sin(Lambda)^2 - 1))/2, -(CLa*Le*c*e*sin(Lambda)*(sin(Lambda)^2 - 1))/12, -(CLa*c*e*sin(Lambda)*(sin(Lambda)^2 - 1))/2,  (CLa*Le*c*e*sin(Lambda)*(sin(Lambda)^2 - 1))/12]];

Ka_wth = [[       (CLa*c*sin(Lambda)*(sin(Lambda)^2 - 1))/2 - (7*CLa*Le*c*cos(Lambda)^2)/20,     (CLa*c*sin(Lambda)*(sin(Lambda)^2 - 1))/2 - (3*CLa*Le*c*cos(Lambda)^2)/20];
[   (CLa*Le^2*c*cos(Lambda)^2)/20 + (CLa*Le*c*sin(Lambda)*(sin(Lambda)^2 - 1))/12, (CLa*Le^2*c*cos(Lambda)^2)/30 - (CLa*Le*c*sin(Lambda)*(sin(Lambda)^2 - 1))/12];
[     - (CLa*c*sin(Lambda)*(sin(Lambda)^2 - 1))/2 - (3*CLa*Le*c*cos(Lambda)^2)/20,   - (CLa*c*sin(Lambda)*(sin(Lambda)^2 - 1))/2 - (7*CLa*Le*c*cos(Lambda)^2)/20];
[ - (CLa*Le^2*c*cos(Lambda)^2)/30 - (CLa*Le*c*sin(Lambda)*(sin(Lambda)^2 - 1))/12, (CLa*Le*c*sin(Lambda)*(sin(Lambda)^2 - 1))/12 - (CLa*Le^2*c*cos(Lambda)^2)/20]];
 
Ka_ww = [[ (CLa*c*sin(2*Lambda))/4 - (6*CLa*c*cos(Lambda)*(cos(Lambda)^2 - 1))/(5*Le),       (CLa*c*cos(Lambda)*(cos(Lambda)^2 - 1))/10 + (CLa*Le*c*sin(2*Lambda))/20,   (6*CLa*c*cos(Lambda)*(cos(Lambda)^2 - 1))/(5*Le) - (CLa*c*sin(2*Lambda))/4,       (CLa*c*cos(Lambda)*(cos(Lambda)^2 - 1))/10 - (CLa*Le*c*sin(2*Lambda))/20];
[   (CLa*c*cos(Lambda)*(cos(Lambda)^2 - 1))/10 - (CLa*Le*c*sin(2*Lambda))/20,                               -(2*CLa*Le*c*cos(Lambda)*(cos(Lambda)^2 - 1))/15,     (CLa*Le*c*sin(2*Lambda))/20 - (CLa*c*cos(Lambda)*(cos(Lambda)^2 - 1))/10, (CLa*c*sin(2*Lambda)*Le^2)/120 + (CLa*c*cos(Lambda)*(cos(Lambda)^2 - 1)*Le)/30];
[ (CLa*c*sin(2*Lambda))/4 + (6*CLa*c*cos(Lambda)*(cos(Lambda)^2 - 1))/(5*Le),     - (CLa*c*cos(Lambda)*(cos(Lambda)^2 - 1))/10 - (CLa*Le*c*sin(2*Lambda))/20, - (CLa*c*sin(2*Lambda))/4 - (6*CLa*c*cos(Lambda)*(cos(Lambda)^2 - 1))/(5*Le),       (CLa*Le*c*sin(2*Lambda))/20 - (CLa*c*cos(Lambda)*(cos(Lambda)^2 - 1))/10];
[   (CLa*c*cos(Lambda)*(cos(Lambda)^2 - 1))/10 + (CLa*Le*c*sin(2*Lambda))/20, (CLa*Le*c*cos(Lambda)*(cos(Lambda)^2 - 1))/30 - (CLa*Le^2*c*sin(2*Lambda))/120,   - (CLa*c*cos(Lambda)*(cos(Lambda)^2 - 1))/10 - (CLa*Le*c*sin(2*Lambda))/20,                               -(2*CLa*Le*c*cos(Lambda)*(cos(Lambda)^2 - 1))/15]];

if(right_wing = 0)
    fa_th = -fa_th; 
    fa_w = -fa_w; 
end

element.sc.fa_th = fa_th; 
element.sc.fa_w = fa_w; 
element.sc.Ka_thth = Ka_thth; 
element.sc.Ka_thw = Ka_thw; 
element.sc.Ka_wth = Ka_wth; 
element.sc.Ka_ww = Ka_ww; 

end 

function beam = beam_staticaero_assembl(beam) 

b.fa_th = zeros(b.nel+1,1); 
b.fa_w = zeros(2*(b.nel+1),1); 
b.Ka_thth = zeros(b.nel+1,b.nel+1); 
b.Ka_thw = zeros(b.nel+1,2*(b.nel+1));
b.Ka_wth = zeros(2*(b.nel+1),b.nel+1);
b.Ka_ww = zeros(2*(b.nel+1),2*(b.nel+1)); 

for i=1:b.nel
    b.el(i) = el_staticaero_assebl(b.el(i),-dot(b.vx,[0, 1, 0]))
    b.fa_th(i:i+1,1) = b.fa_th(i:i+1,1) + b.el(i).fa_th; 
    b.fa_w(2*(i-1)+1:2*(i+1),1) = b.fa_w(2*(i-1)+1:2*(i+1),1) + b.el(i).fa_w;
    b.Ka_thth(i:i+1,i:i+1) = b.Ka_thth(i:i+1,i:i+1) + b.el(i).Ka_thth; 
    b.Ka_thw(i:i+1,2*(i-1)+1:2*(i+1)) = b.Ka_thw(i:i+1,2*(i-1)+1:2*(i+1)) + b.el(i).Ka_thw;
    b.Ka_wth(2*(i-1)+1:2*(i+1),i:i+1) = b.Ka_wth(2*(i-1)+1:2*(i+1),i:i+1) + b.el(i).Ka_wth;
    b.Ka_ww(2*(i-1)+1:2*(i+1),2*(i-1)+1:2*(i+1)) = b.Ka_ww(2*(i-1)+1:2*(i+1),2*(i-1)+1:2*(i+1)) + b.el(i).Ka_ww; 
end
end

