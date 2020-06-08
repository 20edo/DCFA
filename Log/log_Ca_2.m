clear all, clc
syms lambda e c CLa L real 

Ca(:,3) = [0; 
    0; 
    -(13*CLa*c*L*cos(lambda))/35 - (3*pi*c^2*cos(lambda)*sin(lambda)^2)/(20*L) - (CLa*c*e*cos(lambda)*sin(lambda))/2;
    -(pi*c^2*cos(lambda)^2*sin(lambda))/16 - (7*CLa*c*e*L*cos(lambda)^2)/20; 
    (pi*c^2*cos(lambda)*sin(lambda)^2)/80 + (11*CLa*c*L^2*cos(lambda))/210 - (CLa*c*e*L*cos(lambda)*sin(lambda))/10;
    0;
    0;
    0;
    (3*pi*c^2*cos(lambda)*sin(lambda)^2)/(20*L) - (9*CLa*c*L*cos(lambda))/70 + (CLa*c*e*cos(lambda)*sin(lambda))/2; 
    -(pi*c^2*cos(lambda)^2*sin(lambda))/16 - (3*CLa*c*e*L*cos(lambda)^2)/20; 
    (pi*c^2*cos(lambda)*sin(lambda)^2)/80 - (13*CLa*c*L^2*cos(lambda))/420 + (CLa*c*e*L*cos(lambda)*sin(lambda))/10;
    0];

Ca(:,4) = [0; 
    0; 
    -(pi*c^2*cos(lambda)^2*sin(lambda))/16; 
    -(pi*c^2*L*cos(lambda)^3)/24; 
    -(pi*c^2*L*cos(lambda)^2*sin(lambda))/96; 
    0; 
    0;
    0; 
    (pi*c^2*cos(lambda)^2*sin(lambda))/16; 
    -(pi*c^2*L*cos(lambda)^3)/48; 
    (pi*c^2*L*cos(lambda)^2*sin(lambda))/96; 
    0]; 

Ca(:,5) = [0; 
    0; 
    (pi*c^2*cos(lambda)*sin(lambda)^2)/80 + (11*CLa*c*L^2*cos(lambda))/210 + (CLa*c*e*L*cos(lambda)*sin(lambda))/10; 
    (CLa*c*e*L^2*cos(lambda)^2)/20 - (pi*c^2*L*cos(lambda)^2*sin(lambda))/96; 
    -(CLa*c*L^3*cos(lambda))/105 - (pi*c^2*L*cos(lambda)*sin(lambda)^2)/60; 
    0; 
    0; 
    0; 
    (13*CLa*c*cos(lambda)*L^2)/420 - (CLa*c*e*cos(lambda)*sin(lambda)*L)/10 - (c^2*pi*cos(lambda)*sin(lambda)^2)/80; 
    (CLa*c*e*L^2*cos(lambda)^2)/30 + (pi*c^2*L*cos(lambda)^2*sin(lambda))/96; 
    (CLa*c*L^3*cos(lambda))/140 + (pi*c^2*L*cos(lambda)*sin(lambda)^2)/240 - (CLa*c*e*L^2*cos(lambda)*sin(lambda))/60; 
    0]; 

Ca(:,6:8) = zeros(12,3); 

Ca(:,9) = [0; 
    0; 
    (3*pi*c^2*cos(lambda)*sin(lambda)^2)/(20*L) - (9*CLa*c*L+cos(lambda))/70 - (CLa*c*e*cos(lambda)*sin(lambda))/2; 
    (pi*c^2^cos(lambda)^2*sin(lambda))/16 - (3*CLa*c*e*L*cos(lambda)^2)/20; 
    (13*CLa*c*cos(lambda)*L^2)/420 + (CLa*c*e*cos(lambda)*sin(lambda)*L)/10 - (c^2*pi*cos(lambda)*sin(lambda)^2)/80;
    0; 
    0; 
    0; 
    (CLa*c*e*cos(lambda)*sin(lambda))/2 - (3*pi*c^2*cos(lambda)*sin(lambda)^2)/(20*L) - (13*CLa*c*L*cos(lambda))/35; 
    (pi*c^2*cos(lambda)^2*sin(lambda))/16 - (7*CLa*c*e*L*cos(lambda)^2)/20; 
    -(pi*c^2*cos(lambda)*sin(lambda)^2)/80 - (11*CLa*c*L^2*cos(lambda))/210 - (CLa*c*e*L*cos(lambda)*sin(lambda))/10; 
    0]; 

Ca(:,10) = [0; 
    0; 
    -(pi*c^2*cos(lambda)^2*sin(lambda))/16; 
    -(pi*c^2*L*cos(lambda)^3)/48; 
    (pi*c^2*L*cos(lambda)^2*sin(lambda))/96; 
    0; 
    0;
    0; 
    (pi*c^2*cos(lambda)^2*sin(lambda))/16; 
    -(pi*c^2*L*cos(lambda)^3)/24; 
    -(pi*c^2*L*cos(lambda)^2*sin(lambda))/96; 
    0];

Ca(:,11) = [0; 
    0; 
    (pi*c^2*cos(lambda)*sin(lambda)^2)/80 - (13*CLa*c*L^2*cos(lambda))/420 - (CLa*c*e*L*cos(lambda)*sin(lambda))/10; 
    -(CLa*c*e*L^2*cos(lambda)^2)/30 + (pi*c^2*L*cos(lambda)^2*sin(lambda))/96; 
    (CLa*c*L^3*cos(lambda))/140 + (pi*c^2*L*cos(lambda)*sin(lambda)^2)/240 + (CLa*c*e*L^2*cos(lambda)*sin(lambda))/60; 
    0; 
    0; 
    0; 
    -(11*CLa*c*cos(lambda)*L^2)/210 + (CLa*c*e*cos(lambda)*sin(lambda)*L)/10 - (c^2*pi*cos(lambda)*sin(lambda)^2)/80; 
    -(CLa*c*e*L^2*cos(lambda)^2)/20 - (pi*c^2*L*cos(lambda)^2*sin(lambda))/96; 
    -(CLa*c*L^3*cos(lambda))/105 - (pi*c^2*L*cos(lambda)*sin(lambda)^2)/60; 
    0]; 

Ca(:,12) = zeros(12,1); 


    
    

    











