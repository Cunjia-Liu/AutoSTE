function g = drawPlume(f, s, domain)

figure(f);

xmin = domain(1);
xmax = domain(2);

ymin = domain(3);
ymax = domain(4);

zmin = domain(5);
zmax = domain(6);

nGird = 200;

x_coord = linspace(xmin, xmax, nGird);
y_coord = linspace(ymin, ymax, nGird);
z_coord = linspace(zmin, zmax, nGird);

% Create 3D grid
[X,Y,Z] = meshgrid(x_coord,y_coord,z_coord);

ex.x_matrix = X;
ex.y_matrix = Y;
ex.z_matrix = Z;

conc = plumeModel(s,ex);

g = pcolor(ex.x_matrix(:,:,1),ex.y_matrix(:,:,1),conc(:,:,1));
axis([xmin xmax ymin ymax]);

shading interp

xlabel('x');
ylabel('y');
colorbar
hold on

plot(s.x,s.y,'k.','MarkerSize',20)


end