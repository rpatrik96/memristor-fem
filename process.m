%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Memristor FEM simulation - script only
% Reizinger Patrik (2018)
%
% Functionality:
% this function creates evaulates the results obtained from a PDE problem
% solved with memristor_pde.m, and also it creates animations.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% %% Epsilon charakteristics
% syms x;
% fplot((1-1./(1+exp(18e7*(x)))))
% xlim([-1e-7; 1e-7]);
% title('Az anyagjellemzõ súlyozására használt függvény');
% xlabel('Hossz [m]');
% ylabel('Függvényérték');

%% Animation init
n=10;                  % number of pictures
% Bounding rectangle
% A - along x axis
% B - along y axis
rect_x = 4e-08;
rect_y = 2.5e-08;
A = 2*rect_x;
B = 2*rect_y;

% maximua placeholder
maxu = -Inf;
maxE = -Inf;

% material threshold offset
material_boundary = -rect_x:A/n:rect_x;

%% Solve for all material distributions
% and get maxima
for i=1:n
   R(i) = memristor_pde(20e-9, 10e-9, material_boundary(i), 0);
   
   % get max of u
   if R(i).NodalSolution > maxu
       maxu = R(i).NodalSolution;
   end
   
   % get max of E
   maxE_abs = max(norm([R(i).XGradients, R(i).YGradients]));
   if maxE_abs > maxE
       maxE = maxE_abs;
   end
end

% newplot;
% title('Set size then press any key!')
% disp('Set size then press any key!');
% pause;
% maxJ=max(abs_grad);       % maximum a szinskalahoz
% maxu=max(abs(u));       % maximum az indukciovonalakhoz
% blackcmap=zeros(64,3);  % fekete paletta a konturvonalakhoz
% Jmap=load('Jmap.dat');  % szines paletta az aramsuruseghez
% 
% % Kockarol kockara
% for k=1:n
%    pdeplot(p,e,t,'xydata',real(u*exp(j*(k-1)*2*pi/n)),'xystyle','off', ...
%        'contour','on','levels',-maxu:maxu/20:maxu,'colormap',blackcmap);
%    hold;
%    pdeplot(p,e,t,'xydata',real(abs_grad*exp(j*(k-1)*2*pi/n)),'xystyle','flat', ...
%        'contour','off','colormap',Jmap);
%    caxis([-maxJ maxJ]);
%    colorbar;
%    axis([-7.5e-08 7.5e-08 -5e-08 5e-08]), set(gca,'DataAspectRatio',[1 1 1]); axis off 
%    M(k)=getframe(gcf); 
% end
% clf;
% axes('Position',[0 0 1 1]);
% movie(M,50);