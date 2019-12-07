t1=load('./fort.901');
figure;
pcolor(1/(24*365):1/(24*365.):50,1:45,t1');set(gca,'ydir','rev');shading interp;colorbar
ylabel('layers');

print  -depsc 'T1.ps'
figure;

s1=load('./fort.911');
figure;
pcolor(1/(24*365):1/(24*365.):50,1:45,s1');set(gca,'ydir','rev');shading interp;colorbar
ylabel('layers');
xlabel('year');ylabel('layers');
print  -depsc 'S1.ps'

figure;

x1=load('fort.105');
figure;plot(1/(24*365):1/(24.*365):50,x1(:,1));
xlabel('year');ylabel('Vice(m)');
print  -depsc 'Vice.ps'
