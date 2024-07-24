clc;
clear;
close all;

%User Inputs
Size = 1;
ds = 0.001;
x_start = 0;
y_start = -0.5;
z_start = 0;
dir1 = [0;1;0];
dir2 = [0;1;0];
start_2 = 1;
start_3 = 6;
levels = 0;

n_points = 3000; %Arbitrary Value
x = zeros(n_points, n_points);
y = zeros(n_points, n_points);
z = zeros(n_points, n_points);
x(1:2,1) = x_start;
y(1:2,1) = y_start;
z(1:2,1) = z_start;
level_init = f(x_start, y_start, z_start);

count_h = n_points*ones(1,2);
for j = (1:2)  
    for i = (2:n_points)
        P = next_point_plane(x(j,i-1),y(j,i-1),z(j,i-1),level_init,ds,dir1(1),dir1(2),dir1(3),j+1);
        x(j,i) = P(1);
        y(j,i) = P(2);
        z(j,i) = P(3);
        if abs(x(j,i)) > Size/2
            x(j,i) = sign(x(j,i))*Size/2;
            count_h(j) = i-1;
            break
        end
        if abs(y(j,i)) > Size/2
            y(j,i) = sign(y(j,i))*Size/2;
            count_h(j) = i-1;
            break
        end
        if abs(z(j,i)) > Size/2
            z(j,i) = sign(z(j,i))*Size/2;
            count_h(j) = i-1;
            break
        end
    end
end
x(1,1:count_h(1)+count_h(2)-1) = [flip(x(1,2:count_h(1))),x(2,1:count_h(2))];
x = x(:,1:count_h(1)+count_h(2)-1);
x(2,:) = 0;
y(1,1:count_h(1)+count_h(2)-1) = [flip(y(1,2:count_h(1))),y(2,1:count_h(2))];
y = y(:,1:count_h(1)+count_h(2)-1);
y(2,:) = 0;
z(1,1:count_h(1)+count_h(2)-1) = [flip(z(1,2:count_h(1))),z(2,1:count_h(2))];
z = z(:,1:count_h(1)+count_h(2)-1);
z(2,:) = 0;

count_v = n_points*ones(1,width(x));
for i = (1:width(x))
    dir = dir2;
    for j = (2:n_points)
        curv_vect = principal_directions(x(j-1,i),y(j-1,i),z(j-1,i));
        [~,I] = max(abs(transpose(dir)*curv_vect));
        dir = sign(transpose(dir)*curv_vect(:,I))*curv_vect(:,I);
        P = next_point_vector(x(j-1,i),y(j-1,i),z(j-1,i),level_init,ds,dir(1),dir(2),dir(3));
        x(j,i) = P(1);
        y(j,i) = P(2);
        z(j,i) = P(3);
        dir = [x(j,i)-x(j-1,i);y(j,i)-y(j-1,i);z(j,i)-z(j-1,i)];
        dir = dir/sqrt(transpose(dir)*dir);
         if abs(x(j,i)) > Size/2
            count_v(i) = j-1;
            x(j:end,i) = x(j-1,i);
            y(j:end,i) = y(j-1,i);
            z(j:end,i) = z(j-1,i);
            break
        end
        if abs(y(j,i)) > Size/2
            count_v(i) = j-1;
            x(j:end,i) = x(j-1,i);
            y(j:end,i) = y(j-1,i);
            z(j:end,i) = z(j-1,i);
            break
        end
        if abs(z(j,i)) > Size/2
            count_v(i) = j-1;
            x(j:end,i) = x(j-1,i);
            y(j:end,i) = y(j-1,i);
            z(j:end,i) = z(j-1,i);
            break
        end
    end
end
x = x(1:max(count_v),:);
y = y(1:max(count_v),:);
z = z(1:max(count_v),:);

%%
curv_max = ones(1,width(x));
for i = 1:width(x)
    dt_max = 0;
    for j = 2:(count_v(i)-1)
        n = [fx(x(j,i),y(j,i),z(j,i));fy(x(j,i),y(j,i),z(j,i));fz(x(j,i),y(j,i),z(j,i))];
        n = n/sqrt(transpose(n)*n);
        t1 = [x(j,i)-x(j-1,i);y(j,i)-y(j-1,i);z(j,i)-z(j-1,i)];
        t1 = t1/sqrt(transpose(t1)*t1);
        t2 = [x(j+1,i)-x(j,i);y(j+1,i)-y(j,i);z(j+1,i)-z(j,i)];
        t2 = t2/sqrt(transpose(t2)*t2);
        dt = (t2-t1)/ds;
        dt = dt - (transpose(dt)*n)*n;
        dt = sqrt(transpose(dt)*(dt));
        if dt > dt_max
            dt_max = dt;
            curv_max(i) = j;
        end
    end
end

div_indices = zeros(1,width(x));
for i = (2:width(x))
    for j = (2:height(x))
        t1 = [x(j,i-1)-x(j-1,i-1);y(j,i-1)-y(j-1,i-1);z(j,i-1)-z(j-1,i-1)];
        t1 = t1/sqrt(transpose(t1)*t1);
        t2 = [x(j,i)-x(j-1,i);y(j,i)-y(j-1,i);z(j,i)-z(j-1,i)];
        t2 = t2/sqrt(transpose(t2)*t2);
        if transpose(t1)*t2 < 0.5
            div_indices(i-1) = 1;
            div_indices(i) = 1;
            break
        end
    end
end

n = 0;
edge_indices = zeros(1, div_indices*transpose(div_indices)+2-div_indices(1)-div_indices(width(x)));
x_edge = zeros(height(x), width(edge_indices));
y_edge = zeros(height(y),width(x_edge));
z_edge = zeros(height(z),width(x_edge));
if div_indices(1) == 0
    edge_indices(1,1) = 1;
    x_edge(:,1) = x(:,1);
    y_edge(:,1) = y(:,1);
    z_edge(:,1) = z(:,1);
    n = n+1; 
end
for i = 1:width(div_indices)
    if div_indices(i) == 1
        n = n+1;
        edge_indices(1,n) = i;
        x_edge(:,n) = x(:,i);
        y_edge(:,n) = y(:,i);
        z_edge(:,n) = z(:,i);
    end
end
if div_indices(width(x)) == 0
    edge_indices(1,width(edge_indices)) = width(x);
    x_edge(:,width(x_edge)) = x(:,width(x));
    y_edge(:,width(y_edge)) = y(:,width(y));
    z_edge(:,width(z_edge)) = z(:,width(z));
end
clear n;

%%
width_2 = count_v(start_2);
x2 = zeros(n_points, width_2);
x2(1,:) = transpose(x(1:width_2,start_2));
y2 = zeros(n_points, width_2);
y2(1,:) = transpose(y(1:width_2,start_2));
z2 = zeros(n_points, width_2);
z2(1,:) = transpose(z(1:width_2,start_2));
count_v = [count_v, zeros(1,width_2)];
for i = (1:width(x2))
    n = [fx(x(i,start_2),y(i,start_2),z(i,start_2));...
        fy(x(i,start_2),y(i,start_2),z(i,start_2));fz(x(i,start_2),y(i,start_2),z(i,start_2))];
    n = n/sqrt(transpose(n)*n);
    if i < width(x2)
        t = [x(i+1,start_2)-x(i,start_2);y(i+1,start_2)-y(i,start_2);...
        z(i+1,start_2)-z(i,start_2)];
        t = t/sqrt(transpose(t)*t);
    end
    dir = [t(2,1)*n(3,1)-t(3,1)*n(2,1);t(3,1)*n(1,1)-t(1,1)*n(3,1);...
        t(1,1)*n(2,1)-t(2,1)*n(1,1)];
    for j = (2:n_points)
        curv_vect = principal_directions(x2(j-1,i),y2(j-1,i),z2(j-1,i));
        [~,I] = max(abs(transpose(dir)*curv_vect));
        dir = sign(transpose(dir)*curv_vect(:,I))*curv_vect(:,I);
        P = next_point_vector(x2(j-1,i),y2(j-1,i),z2(j-1,i),level_init,ds,dir(1),dir(2),dir(3));
        x2(j,i) = P(1);
        y2(j,i) = P(2);
        z2(j,i) = P(3);
        dir = [x2(j,i)-x2(j-1,i);y2(j,i)-y2(j-1,i);z2(j,i)-z2(j-1,i)];
        dir = dir/sqrt(transpose(dir)*dir);
        if and(j>2, is_intersecting(j-1,i,x2,y2,z2,x_edge,y_edge,z_edge,count_v(edge_indices),ds))
            count_v(width(x)+i) = j-1;
            x2(j:end,i) = x2(j-1,i);
            y2(j:end,i) = y2(j-1,i);
            z2(j:end,i) = z2(j-1,i);
            break
        end
        if abs(x2(j,i)) > Size/2
            count_v(width(x)+i) = j-1;
            x2(j:end,i) = x2(j-1,i);
            y2(j:end,i) = y2(j-1,i);
            z2(j:end,i) = z2(j-1,i);
            break
        end
        if abs(y2(j,i)) > Size/2
            count_v(width(x)+i) = j-1;
            x2(j:end,i) = x2(j-1,i);
            y2(j:end,i) = y2(j-1,i);
            z2(j:end,i) = z2(j-1,i);
            break
        end
        if abs(z2(j,i)) > Size/2
            count_v(width(x)+i) = j-1;
            x2(j:end,i) = x2(j-1,i);
            y2(j:end,i) = y2(j-1,i);
            z2(j:end,i) = z2(j-1,i);
            break
        end
    end
end
x2 = x2(1:max(count_v(width(x)+1:end)),:);
y2 = y2(1:max(count_v(width(x)+1:end)),:);
z2 = z2(1:max(count_v(width(x)+1:end)),:);

%%
if ~(start_3 == 0)
    id1 = edge_indices(start_3);
    if start_3+1 < width(edge_indices)
        id2 = edge_indices(start_3+1);
    else
        id2 = 0;
    end
    if id2 == id1+1
        width_3 = count_v(id1)+count_v(id2)-curv_max(id1)-curv_max(id2)+2;
        x3 = zeros(n_points, width_3);
        x3(1,:) = transpose([flip(x(curv_max(id1):count_v(id1),id1));...
            x(curv_max(id2):count_v(id2),id2)]);
        y3 = zeros(n_points, width_3);
        y3(1,:) = transpose([flip(y(curv_max(id1):count_v(id1),id1));...
           y(curv_max(id2):count_v(id2),id2)]);
        z3 = zeros(n_points, width_3);
        z3(1,:) = transpose([flip(z(curv_max(id1):count_v(id1),id1));...
            z(curv_max(id2):count_v(id2),id2)]);
    else
        width_3 = count_v(id1)-curv_max(id1)+1;
        x3 = zeros(n_points, width_3);
        x3(1,:) = transpose(x(curv_max(id1):count_v(id1),id1));
        y3 = zeros(n_points, width_3);
        y3(1,:) = transpose(y(curv_max(id1):count_v(id1),id1));
        z3 = zeros(n_points, width_3);
        z3(1,:) = transpose(z(curv_max(id1):count_v(id1),id1));
    end
    
    count_v = [count_v, zeros(1,width_3)];
    for i = (1:width(x3))
        n = [fx(x3(1,i),y3(1,i),z3(1,i));fy(x3(1,i),y3(1,i),z3(1,i));fz(x3(1,i),y3(1,i),z3(1,i))];
        n = n/sqrt(transpose(n)*n);
        if i < width(x3)
            t = [x3(1,i+1)-x3(1,i);y3(1,i+1)-y3(1,i);...
            z3(1,i+1)-z3(1,i)];
            t = t/sqrt(transpose(t)*t);
        end
        dir = [t(2,1)*n(3,1)-t(3,1)*n(2,1);t(3,1)*n(1,1)-t(1,1)*n(3,1);...
            t(1,1)*n(2,1)-t(2,1)*n(1,1)];
        for j = (2:n_points)
            curv_vect = principal_directions(x3(j-1,i),y3(j-1,i),z3(j-1,i));
            [~,I] = max(abs(transpose(dir)*curv_vect));
            dir = sign(transpose(dir)*curv_vect(:,I))*curv_vect(:,I);
            P = next_point_vector(x3(j-1,i),y3(j-1,i),z3(j-1,i),level_init,ds,dir(1),dir(2),dir(3));
            x3(j,i) = P(1);
            y3(j,i) = P(2);
            z3(j,i) = P(3);
            dir = [x3(j,i)-x3(j-1,i);y3(j,i)-y3(j-1,i);z3(j,i)-z3(j-1,i)];
            dir = dir/sqrt(transpose(dir)*dir);
            if and(j>2, is_intersecting(j-1,i,x3,y3,z3,x_edge,y_edge,z_edge,count_v(edge_indices),ds))
                count_v(width(x)+width(x2)+i) = j-1;
                x3(j:end,i) = x3(j-1,i);
                y3(j:end,i) = y3(j-1,i);
                z3(j:end,i) = z3(j-1,i);
                break
            end
            if abs(x3(j,i)) > Size/2
                count_v(width(x)+width(x2)+i) = j-1;
                x3(j:end,i) = x3(j-1,i);
                y3(j:end,i) = y3(j-1,i);
                z3(j:end,i) = z3(j-1,i);
                break
            end
            if abs(y3(j,i)) > Size/2
                count_v(width(x)+width(x2)+i) = j-1;
                x3(j:end,i) = x3(j-1,i);
                y3(j:end,i) = y3(j-1,i);
                z3(j:end,i) = z3(j-1,i);
                break
            end
            if abs(z3(j,i)) > Size/2
                count_v(width(x)+width(x2)+i) = j-1;
                x3(j:end,i) = x3(j-1,i);
                y3(j:end,i) = y3(j-1,i);
                z3(j:end,i) = z3(j-1,i);
                break
            end
        end
    end
end
x3 = x3(1:max(count_v(width(x)+width(x2)+1:end)),:);
y3 = y3(1:max(count_v(width(x)+width(x2)+1:end)),:);
z3 = z3(1:max(count_v(width(x)+width(x2)+1:end)),:);

%%
indices = zeros(1, length(levels));
for i = (1:length(levels))
    value = find_level(x_start,y_start,z_start,levels(i),ds,1,0,0);
    indices(i) = value(1);
end

for i = (1:width(x))
    for j = (2:length(levels))
        P = gradiant_move(x(j-1,i),y(j-1,i),z(j-1,i),levels(j),ds);
        x(j,i) = P(1);
        y(j,i) = P(2);
        z(j,i) = P(3);
    end
end

%%
id1 = 1;
id2 = 1;
vectors = 0.05*principal_directions(x(id1,id2),y(id1,id2),z(id1,id2));
n = [fx(x(id1,id2),y(id1,id2),z(id1,id2));...
    fy(x(id1,id2),y(id1,id2),z(id1,id2));fz(x(id1,id2),y(id1,id2),z(id1,id2))];
n = 0.05*n/sqrt(transpose(n)*n);
vectors = [vectors,n];

n = 100;
[X, Y, Z] = meshgrid(-Size/2:Size/n:Size/2);
OBJ = cos(2*pi*X/Size).* sin(2*pi*Y/Size) + cos(2*pi*Y/Size).* sin(2*pi*Z)...
        + cos(2*pi*Z).* sin(2*pi*X/Size);

figure(1)
isosurface(X,Y,Z, (OBJ), 0);

xlim([-0.5*Size 0.5*Size])
ylim([-0.5*Size 0.5*Size])
zlim([-0.5*Size 0.5*Size])
hold on;
for i = 1:3
    line([x(id1,id2),x(id1,id2)+vectors(1,i)],...
        [y(id1,id2),y(id1,id2)+vectors(2,i)],[z(id1,id2),z(id1,id2)+vectors(3,i)])
end
plot3(x(1,:), y(1,:), z(1,:), "-o", "Color","b")
for i = (1:width(x))
    c = "k";
    if div_indices(i) == 1
        j = curv_max(i);
        scatter3(x(j,i),y(j,i),z(j,i),6,[0 1 0])
        %p1 = find_divergence_point(x(j,i),y(j,i),z(j,i),level_init,0.001);
        %scatter3(p1(1),p1(2),p1(3),6,[1 1 0])
        c = "r";
    end
    plot3(x(:,i), y(:,i), z(:,i), "-", "Color",c)
end
plot3(x2(1,:), y2(1,:), z2(1,:), "-o", "Color","b")
for i = (1:width(x2))
    plot3(x2(:,i), y2(:,i), z2(:,i), "-", "Color","r")
end
plot3(x3(1,:), y3(1,:), z3(1,:), "-o", "Color","b")
for i = (1:width(x3))
    plot3(x3(:,i), y3(:,i), z3(:,i), "-", "Color","g")
end
hold off;
%%

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

function value = fz(x,y,z) 
    value = 2*pi*cos(2*pi*y)*cos(2*pi*z) - 2*pi*sin(2*pi*z)*sin(2*pi*x);
end

function value = fxx(x,y,z) 
    value = -4*pi^2*cos(2*pi*x)*sin(2*pi*y) - 4*pi^2*cos(2*pi*z)*sin(2*pi*x);
end

function value = fyy(x,y,z) 
    value = -4*pi^2*cos(2*pi*x)*sin(2*pi*y) - 4*pi^2*cos(2*pi*y)*sin(2*pi*z);
end

function value = fzz(x,y,z) 
    value = -4*pi^2*cos(2*pi*y)*sin(2*pi*z) - 4*pi^2*cos(2*pi*z)*sin(2*pi*x);
end

function value = fxy(x,y,~) 
    value = -4*pi^2*sin(2*pi*x)*cos(2*pi*y);
end

function value = fyz(~,y,z) 
    value = -4*pi^2*sin(2*pi*y)*cos(2*pi*z);
end

function value = fzx(x,~,z) 
    value = -4*pi^2*sin(2*pi*z)*cos(2*pi*x);
end

function value = find_level(x,y,z,const,ds,i,j,k)
    epsilon = 1e-12;
    fval = f(x,y,z) - const;
    while abs(fval) > epsilon
        f_dot = (fx(x,y,z)*i + fy(x,y,z)*j + fz(x,y,z)*k);
        dx = -ds*fval*i/f_dot;
        dy = -ds*fval*j/f_dot;
        dz = -ds*fval*k/f_dot;
        fval = f(x+dx,y+dy,z) - const;
        x = x+dx;
        y = y+dy;
        z = z+dz;
    end
    value = [x, y, z]; 
end

function value = Hessian(x,y,z)
    x_ = fx(x,y,z);
    y_ = fy(x,y,z);
    z_ = fz(x,y,z);
    xx = fxx(x,y,z);
    yy = fyy(x,y,z);
    zz = fzz(x,y,z);
    xy = fxy(x,y,z);
    yz = fyz(x,y,z);
    zx = fzx(x,y,z);
    denominator = (x_^2+y_^2+z_^2)^(3/2);
    H11 = xx*y_^2 - x_*y_*xy + xx*z_^2 - x_*z_*zx; 
    H12 = xy*y_^2 - x_*y_*yy + xy*z_^2 - x_*z_*yz; 
    H13 = zx*y_^2 - x_*y_*yz + zx*z_^2 - x_*z_*zz;
    H21 = xy*x_^2 - x_*y_*xx + xy*z_^2 - y_*z_*zx;
    H22 = yy*x_^2 - x_*y_*xy + yy*z_^2 - y_*z_*yz;
    H23 = yz*x_^2 - x_*y_*zx + yz*z_^2 - y_*z_*zz;
    H31 = zx*x_^2 - x_*z_*xx + zx*y_^2 - y_*z_*xy;
    H32 = yz*x_^2 - x_*z_*xy + yz*y_^2 - y_*z_*yz; 
    H33 = zz*x_^2 - x_*z_*zx + zz*y_^2 - y_*z_*yz;
    value = [H11,H12,H13;H21,H22,H23;H31,H32,H33]/denominator;
end

function value = principal_directions(x,y,z)
    H = Hessian(x,y,z);
    n = [fx(x,y,z);fy(x,y,z);fz(x,y,z)];
    n = n/sqrt(transpose(n)*n);
    [eig_vect,~] = eig(H);
    [~,I] = max(abs(transpose(n)*eig_vect));
    value = eig_vect(:,[mod(I+1,3)+1,mod(I,3)+1]);
    for i = 1:2
        value(:,i) = value(:,i) - (transpose(value(:,i))*n)*n;
        value(:,i) = value(:,i)/sqrt(transpose(value(:,i))*value(:,i)); 
    end
end

function value = next_point_plane(x,y,z,const,ds,i,j,k,dir)
    epsilon = 1e-12;
    n = [i,j,k]/sqrt(i^2+j^2+k^2);
    g = [fx(x,y,z),fy(x,y,z),fz(x,y,z)];
    g = g/sqrt(g*transpose(g));
    ds_cap = [n(2)*g(3)-n(3)*g(2),n(3)*g(1)-n(1)*g(3),...
        n(1)*g(2)-n(2)*g(1)];
    ds_cap = ds_cap/sqrt(ds_cap*transpose(ds_cap));
    dx = ds*ds_cap(1);
    dy = ds*ds_cap(2);
    dz = ds*ds_cap(3);
    xdash = x + ((-1)^dir)*dx;
    ydash = y + ((-1)^dir)*dy;
    zdash = z + ((-1)^dir)*dz;
    fdash = f(xdash,ydash,zdash)-const;
    while (abs(fdash) > epsilon)
        g = [fx(x,y,z),fy(x,y,z),fz(x,y,z)];
        g_plane = g - (g*transpose(n))*n;
        g_plane = g_plane/(g_plane*transpose(g_plane));
        xdash_new = xdash - fdash*g_plane(1);
        ydash_new = ydash - fdash*g_plane(2);
        zdash_new = zdash - fdash*g_plane(3);
        xdash = xdash_new;
        ydash = ydash_new;
        zdash = zdash_new;
        fdash = f(xdash,ydash,zdash)-const;
    end
    value = [xdash, ydash, zdash];
end

function value = next_point_vector(x,y,z,const,ds,i,j,k)
    epsilon = 1e-12;
    v = [i,j,k]/sqrt(i^2+j^2+k^2);
    n = [fx(x,y,z),fy(x,y,z),fz(x,y,z)];
    n = n/sqrt(n*transpose(n));
    ds_cap = v - (n*transpose(v))*n;
    ds_cap = ds_cap/sqrt(ds_cap*transpose(ds_cap));
    dx = ds*ds_cap(1);
    dy = ds*ds_cap(2);
    dz = ds*ds_cap(3);
    x = x + dx;
    y = y + dy;
    z = z + dz;
    fdash = f(x,y,z)-const;
    while (abs(fdash) > epsilon)
        n = [fx(x,y,z),fy(x,y,z),fz(x,y,z)];
        n = n/(n*transpose(n));
        x = x - fdash*n(1);
        y = y- fdash*n(2);
        z = z - fdash*n(3);
        fdash = f(x,y,z)-const;
    end
    value = [x, y, z];
end

function value = gradiant_move(x,y,z,level,ds)
    epsilon = 1e-12;
    fval = f(x,y,z) - level;
    while abs(fval) > epsilon
        f_dot = (fx(x,y,z)^2 + fy(x,y,z)^2 + fz(x,y,z)^2);
        dx = -ds*fval*fx(x,y,z)/f_dot;
        dy = -ds*fval*fy(x,y,z)/f_dot;
        dz = -ds*fval*fz(x,y,z)/f_dot;
        fval = f(x+dx,y+dy,z) - level;
        x = x+dx;
        y = y+dy;
        z = z+dz;
    end
    value = [x, y, z];
end

function p1 = find_divergence_point(x,y,z,const,ds)
    epsilon_fraction = 1e-3;
    epsilon_absolute = 1e-15;
    p1 = [x,y,z];
    p1 = gradiant_move(p1(1),p1(2),p1(3),const,0.001);
    H1 = det(Hessian(p1(1),p1(2),p1(3)));
    convergence = false;
    while ~convergence
        H1_old = H1;
        n_ = [fx(p1(1),p1(2),p1(3));fy(p1(1),p1(2),p1(3));...
        fz(p1(1),p1(2),p1(3))];
        n_ = n_/sqrt(transpose(n_)*n_);
        d_ = principal_directions(p1(1),p1(2),p1(3));
        p2 = p1 + ds*transpose(d_(:,1));
        p2 = gradiant_move(p2(1,1),p2(1,2),p2(1,3),const,ds);
        H2 = det(Hessian(p2(1),p2(2),p2(3)));
        p3 = p1 + ds*transpose(d_(:,2));
        p3 = gradiant_move(p3(1,1),p3(1,2),p3(1,3),const,ds);
        H3 = det(Hessian(p3(1),p3(2),p3(3)));
        g_ = [0,H2-H1,H3-H1]/[n_,transpose(p2-p1),transpose(p3-p1)];
        p1 = p1 - H1*g_/(g_*transpose(g_));
        p1 = gradiant_move(p1(1),p1(2),p1(3),const,0.001);
        H1 = det(Hessian(p1(1),p1(2),p1(3)));
        convergence = or(abs((H1-H1_old)/H1_old) < epsilon_fraction,...
            abs(H1) < epsilon_absolute);
    end
end

function value = is_intersecting(j,i,x2,y2,z2,x,y,z,count_v,ds)
    p1 = [x2(j,i),y2(j,i),z2(j,i)];
    p2 = [x2(j+1,i),y2(j+1,i),z2(j+1,i)];
    p12 = (p1+p2)/2;
    n_ = [fx(p12(1),p12(2),p12(3));fy(p12(1),p12(2),p12(3));...
        fz(p12(1),p12(2),p12(3))];
    t_ = [(p2(1)-p1(1)),(p2(2)-p1(2)),(p2(3)-p1(3))];
    nt_ = [n_(2)*t_(3)-n_(3)*t_(2),n_(3)*t_(1)-n_(1)*t_(3),...
        n_(1)*t_(2)-n_(2)*t_(1)];
    nt_ = nt_/sqrt(nt_*transpose(nt_));
    mid_plane = @(p_) (nt_(1)*(p_(1)-p12(1))...
        + nt_(3)*(p_(2)-p12(2)) + nt_(2)*(p_(3)-p12(3)));
    x_arr = abs(x - p12(1));
    y_arr = abs(y - p12(2));
    z_arr = abs(z - p12(3));
    d_arr = sqrt(x_arr.^2+y_arr.^2+z_arr.^2);
    [M,~] =  min(d_arr,[],"all","linear");
    if M > ds
        value = false;
        return 
    end
    for i_ = 1:width(x)
        x_arr = abs(x(:,i_) - p12(1));
        y_arr = abs(y(:,i_) - p12(2));
        z_arr = abs(z(:,i_) - p12(3));
        d_arr = sqrt(x_arr.^2+y_arr.^2+z_arr.^2);
        [M,I] =  min(d_arr,[],"all","linear");
        if M > ds
            continue
        end
        for j_ = max(I-2,1):min(I+2,count_v(i_)-1)
            if mid_plane([x(j_,i_),y(j_,i_),z(j_,i_)])...
                    *mid_plane([x(j_+1,i_),y(j_+1,i_),z(j_+1,i_)]) <= 0
                value = true;
                return 
            end
        end
    end
    value = false;
end
