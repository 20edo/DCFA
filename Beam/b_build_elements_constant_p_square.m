function b=b_build_elements_constant_p_square(b)
% This function builds the elements of the beam b.
% Required fields are:
% The geometry properly defined
% b.E
% b.G
% b.rho
% b.nel
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
%% Check if the input beam satisfies the requirements
nin=b.nel+1;
if isa(b.E,'function_handle')==false
     error('E not well defined')
end
if isa(b.G,'function_handle')==false
     error('G not well defined')
end
if isa(b.rho,'function_handle')==false
     error('rho not well defined')
end
if isnan(b.in)==true
    b.in=[];
    for j=linspace(0,b.L,nin)
        node=in_init();
        node.x=j;
        b.in=[b.in node];
    end    
end
%% Build elements
b.el=[];
for i=1:length(b.in)-1
    b.el=[b.el el_build_constant_p_square(b,b.in(i),b.in(i+1))];
end
end