clear;
close all;

Size = 1;
n = 100;

Z = 0;
[X,Y] = meshgrid(-Size/2:Size/n:Size/2);
OBJ = cos(2*pi*X/Size).* sin(2*pi*Y/Size) + cos(2*pi*Y/Size)* sin(2*pi*Z)...
        + cos(2*pi*Z).* sin(2*pi*X/Size);

levels = (-0.95:0.1:0.95);
indices = zeros(1, length(levels));
x_start = -0.2;
for i = (1:length(levels))
    value = find_level(x_start,0.5,Z,levels(i),0.001,1,0);
    indices(i) = value(1);
end

n_points = 3000; %Arbitrary Value
x = zeros(length(levels), n_points);
y = zeros(length(levels), n_points);
x(:,1) = transpose(indices);
y(:,1) = 0.5;
z = Z;
for i = (2:n_points)
    P = next_point(x(1,i-1),y(1,i-1),z,levels(1),0.01,0);
    x(1,i) = P(1);
    y(1,i) = P(2);
    if abs(x(1,i)) > Size/2
        x(1,i) = sign(x(1,i))*Size/2;
        x = x(:,1:i);
        y = y(:,1:i);
        break
    end
    if abs(y(1,i)) > Size/2
        y(1,i) = sign(y(1,i))*Size/2;
        x = x(:,1:i);
        y = y(:,1:i);
        break
    end
end

for i = (1:width(x))
    for j = (2:length(levels))
        P = perpendicular_move(x(j-1,i),y(j-1,i),z,levels(j),0.001);
        x(j,i) = P(1);
        y(j,i) = P(2);
    end
end
%%
figure(1)
contourf(X,Y,(OBJ),30)

figure(2)
contour(X,Y, (OBJ), [-1,-0.5,0,0.5,1], "ShowText", "on")
xlim([-0.5 0.5])
ylim([-0.5 0.5])
hold on;
for j = (1:height(x))
    plot(x(j,:), y(j,:), "-", "Color","b")
end
for i = (1:width(x))
    line(x(:,i), y(:,i), "Color","k")
end
hold off;
%%

coords = zeros(3,width(x)*height(x));
for i = (1:height(x))
    coords(:,(i-1)*width(x)+1:i*width(x)) = [x(i,:);y(i,:);z*ones(1,width(x))];
end
rectangle_array = zeros((width(x)-1)*(height(x)-1),5);
for i = (1:height(x)-1)
    for j = (1:width(x)-1)
        count = (i-1)*(width(x)-1)+j;
        rectangle_array(count,:) = [count,(i-1)*width(x)+j,(i-1)*width(x)+j+1,...
            i*width(x)+j+1,i*width(x)+j];
    end
end

coords = transpose(coords);
writematrix(transpose(1:height(coords)),"coords_indices","Delimiter","space")
writematrix(coords,"coords","Delimiter","space")
writematrix(rectangle_array,"rectangles","Delimiter","space")

function value = f(x,y,z)
    value = cos(2*pi*x)*sin(2*pi*y)+cos(2*pi*y)*sin(2*pi*z)...
        +cos(2*pi*z)*sin(2*pi*x);
end

function value = fx(x,y,z) 
    value = -2*pi*sin(2*pi*x)*sin(2*pi*y) + 2*pi*cos(2*pi*z)*cos(2*pi*x);
end

function value = fy(x,y,z) 
    value = 2*pi*cos(2*pi*x)*cos(2*pi*y) - 2*pi*sin(2*pi*y)*sin(2*pi*z);
end

function value = dy_dx(x,y,z)
    value = -fx(x,y,z)/fy(x,y,z);
end

function value = find_level(x,y,z,level,ds,i,j)
    epsilon = 1e-6;
    fval = f(x,y,z) - level;
    while abs(fval) > epsilon
        f_dot = (fx(x,y,z)*i + fy(x,y,z)*j);
        dx = -ds*fval*i/f_dot;
        dy = -ds*fval*j/f_dot;
        fval = f(x+dx,y+dy,z) - level;
        x = x+dx;
        y = y+dy;
    end
    value = [x, y]; 
end

function value = x_next(x,y,z,dy)
    const = f(x,y,z);
    epsilon = 1e-6;
    xdash = x + dy/dy_dx(x,y,z);
    fxdash = f(xdash,y+dy,z)-const;
    while (abs(fxdash) > epsilon)
        xdash = xdash - fxdash/fx(xdash,y+dy,z);
        fxdash = f(xdash,y+dy,z)-const;
    end
    value = xdash;
end

function value = next_point(x,y,z,const,ds,dir)
    epsilon = 1e-6;
    dx = ds*fy(x,y,z)/sqrt(fx(x,y,z)^2 + fy(x,y,z)^2);
    dy = dx*dy_dx(x,y,z);
    xdash = x + ((-1)^dir)*dx;
    ydash = y + ((-1)^dir)*dy;
    fdash = f(xdash,ydash,z)-const;
    while (abs(fdash) > epsilon)
        xdash_new = xdash - fdash*fx(x,y,z)/(fx(x,y,z)^2 + fy(x,y,z)^2);
        ydash_new = ydash - fdash*fy(x,y,z)/(fx(x,y,z)^2 + fy(x,y,z)^2);
        xdash = xdash_new;
        ydash = ydash_new;
        fdash = f(xdash,ydash,z)-const;
    end
    value = [xdash, ydash];
end

function value = perpendicular_move(x,y,z,level,ds)
    epsilon = 1e-6;
    fval = f(x,y,z) - level;
    while abs(fval) > epsilon
        f_dot = (fx(x,y,z)^2 + fy(x,y,z)^2);
        dx = -ds*fval*fx(x,y,z)/f_dot;
        dy = -ds*fval*fy(x,y,z)/f_dot;
        fval = f(x+dx,y+dy,z) - level;
        x = x+dx;
        y = y+dy;
    end
    value = [x, y];
end