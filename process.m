%% Epszilon karakterisztika
syms x;
fplot((1-1./(1+exp(18e7*(x)))))
xlim([-1e-7; 1e-7]);
title('Az anyagjellemzõ súlyozására használt függvény');
xlabel('Hossz [m]');
ylabel('Függvényérték');

%% Animáció készítése

% gradiens
[ux,uy] = pdegrad(p,t,u);
abs_grad = sqrt(ux.^2+uy.^2);

% ide az kellene, hogy c-t kiértékeljük a háromszögeken, amihez a
% háromszöghálót griddé kell konvertálni
%HA MINDEN IGAZ, VAN RÁ FV, TEGNAP LÁTTAM

% Az animacio elokeszitese
n=10;                  % a kepek szama
newplot;
title('Allitsa be az animacio meretet, majd nyomjon le egy billentyut!')
disp('Allitsa be az animacio meretet, majd nyomjon le egy billentyut!');
pause;
maxJ=max(abs_grad);       % maximum a szinskalahoz
maxu=max(abs(u));       % maximum az indukciovonalakhoz
blackcmap=zeros(64,3);  % fekete paletta a konturvonalakhoz
Jmap=load('Jmap.dat');  % szines paletta az aramsuruseghez

% Kockarol kockara
for k=1:n
   pdeplot(p,e,t,'xydata',real(u*exp(j*(k-1)*2*pi/n)),'xystyle','off', ...
       'contour','on','levels',-maxu:maxu/20:maxu,'colormap',blackcmap);
   hold;
   pdeplot(p,e,t,'xydata',real(abs_grad*exp(j*(k-1)*2*pi/n)),'xystyle','flat', ...
       'contour','off','colormap',Jmap);
   caxis([-maxJ maxJ]);
   colorbar;
   axis([-7.5e-08 7.5e-08 -5e-08 5e-08]), set(gca,'DataAspectRatio',[1 1 1]); axis off 
   M(k)=getframe(gcf); 
end
clf;
axes('Position',[0 0 1 1]);
movie(M,50);