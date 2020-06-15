%%
% DCFA swept wing assignement
%
% Teamwork
% Team members: Venti Edoardo         944421
%               Zemello Matteo        942003
%               Zucchelli Umberto     952952
%               
%         
%
cd ..
cd generate_model\
generate_model
cd Eigenmodes\

if 1
    options.plot_original          = 1;
    options.plot_deformed          = 1;
    options.plotColor              = 'green';
    options.saveSTL                = 0;
    options.point_section          = 8;
    options.N                      = 31;
    options.saveImages             = 1;
    m_Modes3d(aircraft,options);
end

beep