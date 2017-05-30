[X, Y] = meshgrid(0:.1:20, 0:.1:20);
x = [ 5 15 10 7 9 2];
y = [ 10 10 17 4 13 18];
q = [ 2 1 -2 2 -3 2];

U = zeros(201, 201, 2);
UU = zeros(201, 201);

for i = 1:6
    r = sqrt((X-x(i)).^2 + (Y-y(i)).^2 + 0.01);   % to avoid NaN
    U(:,:,i) = q(i)./r;
    UU = UU + U(:,:,i);
end

[dx, dy] = gradient(-UU); 
%quiver(X, Y, dx, dy,2 );


% define the starting point of stremlines (electric field lines)
% originating from positive charges

th = linspace(0,360,19)*pi/180;
x0 = 0.1 * cos(th); 
y0 = 0.1 * sin(th); 

x1 = 1 * cos(th);
y1 = 1 * sin(th);

hold on 

for i = 1:6
    if q(i) > 0 
    streamline(X, Y, dx, dy, x0 + x(i), y0 + y(i));
    end       
end

% plot the boundaries, starting points picked by trial and error

streamline(X, Y, dx, dy, 13.20, 6.73)
streamline(X, Y, dx, dy, 13.19, 6.75)

streamline(X, Y, dx, dy, 1.85, 13.36)
streamline(X, Y, dx, dy, 1.84, 13.36)

streamline(X, Y, dx, dy, 5.58, 7.05)
streamline(X, Y, dx, dy, 5.58, 7.04)

hold off