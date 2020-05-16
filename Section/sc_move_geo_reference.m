function sc = sc_move_geo_reference(sc)
% This function moves the reference geometry to obtain yct=zct=0. No
% property of the section is changed, except geometrical definition

%%
% DCFA swept wing assignement
%
% Teamwork
% Team members: Pasturenzi Lorenzo    944610
%               Tacchi Alberto        944579
%               Venti Edoardo         944421
%               Zemello Matteo        942003
%               Zucchelli Umberto     952952
%               
%           
%

if sc.cart==true
    sc.Ymax=sc.Ymax-sc.yct;
    sc.Ymin=sc.Ymin-sc.yct;
    sc.Zmax=@(y) sc.Zmax(y)-sc.zct;
    sc.Zmin=@(y) sc.Zmin(y)-sc.zct;
    
elseif sc.pol==true
    
    % Due to definition for each th in the 0,0 reference there is one and
    % only one rho
    sc.Thmin=atan2(sc.Rhomax(sc.Thmin)*cos(sc.Thmin)+sc.yct,sc.Rhomax(sc.Thmin)*sin(sc.Thmin)+sc.zct);
    sc.Thmax=atan2(sc.Rhomax(sc.Thmax)*cos(sc.Thmax)+sc.yct,sc.Rhomax(sc.Thmax)*sin(sc.Thmax)+sc.zct);
    
    rhoct=norm([sc.yct,sc.zct],2);
    thct=atan2(sc.yct,sc.zct);
    
    sc.Rhomax=@(th) abs(sc.Rhomax(th)*exp(1j*th)-rhoct*exp(1j*thct));
    sc.Rhomin=@(th) abs(sc.Rhomin(th)*exp(1j*th)-rhoct*exp(1j*thct));
else
    error('Geometry not well defined');
end

end

