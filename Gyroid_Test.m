clc;
clear;

Size = 1;
P_start = [0,-0.5,0];
f_init = 0;
n_points = 15000;
ds = 0.001;

x = zeros(1,n_points);
x(1) = P_start(1);
y = zeros(1,n_points);
y(1) = P_start(2);
z = zeros(1,n_points);
z(1) = P_start(3);
P_dash_array = zeros(n_points,3);
Px = @(x,y,z) (x-Size/2)*(x+Size/2);
Py = @(x,y,z) (y-Size/2)*(y+Size/2);
Pz = @(x,y,z) (z-Size/2)*(z+Size/2);

dir = 1;
init_normal = [0,1,0];
normal = init_normal;
edge_loop = false;
normal_new = false;
find_P_dash = true;
for i = 1:width(x)-1
    if find_P_dash
        P = [x(i),y(i),z(i)];
        for j = 1:3
            if ~and(abs(P(mod(j-1,3)+1)) >= Size/2, abs(P(mod(j,3)+1)) >= Size/2)
                P_dash = P-sign(P*transpose(normal))*Size*normal/sqrt(normal*transpose(normal));
                P_dash_array(i,:) = P_dash;
                break
            end
        end
    end
    find_P_dash = true;
    if ~normal_new
        P_next = next_point_plane(x(i),y(i),z(i),f_init,ds,normal(1),normal(2),normal(3),dir);
    else
        P_next = next_point_plane(x(i),y(i),z(i),f_init,ds,normal(1),normal(2),normal(3),dir);
        if max(abs(P_next)) > Size/2
            dir = 1-dir;
        end
        P_next = next_point_plane(x(i),y(i),z(i),f_init,ds,normal(1),normal(2),normal(3),dir);
        normal_new = false;
    end
    x(i+1) = P_next(1);
    y(i+1) = P_next(2);
    z(i+1) = P_next(3);
    if and(abs(normal) == abs(init_normal),...
            [x(1),y(1),z(1)]*transpose(abs(init_normal)) == [x(i+1),y(i+1),z(i+1)]*transpose(abs(normal)))
        if and(i > 2, (x(i+1)-x(1))^2+(y(i+1)-y(1))^2+(z(i+1)-z(1))^2 <= ds^2)
            x(i+1) = x(1);
            x = x(1:i+1);
            y(i+1) = y(1);
            y = y(1:i+1);
            z(i+1) = z(1);
            z = z(1:i+1);
            break
        end
    end
    if abs(x(i+1)) > Size/2
        edge_loop = true;
        x(i+1) = Size/2*sign(x(i+1));
    end
    if abs(y(i+1)) > Size/2
        edge_loop = true;
        y(i+1) = Size/2*sign(y(i+1));
    end
    if abs(z(i+1)) > Size/2
        edge_loop = true;
        z(i+1) = Size/2*sign(z(i+1));
    end
    if edge_loop
        face_array = zeros(1,3);
        face_array(1) = Px(x(i+1),y(i+1),z(i+1));
        face_array(2) = Py(x(i+1),y(i+1),z(i+1));
        face_array(3) = Pz(x(i+1),y(i+1),z(i+1));
        edge = zeros(1,3);
        [M,I] = max(abs(face_array));
        if M ~= 0
            edge(I) = 1;
            P_next = find_level(x(i+1),y(i+1),z(i+1),f_init,ds,edge(1),edge(2),edge(3));
            x(i+1) = P_next(1);
            y(i+1) = P_next(2);
            z(i+1) = P_next(3);
            normal = [normal(2)*edge(3)-normal(3)*edge(2),...
                normal(3)*edge(1)-normal(1)*edge(3),normal(1)*edge(2)-normal(2)*edge(1)];
        else
            [~,I] = max(abs(normal));
            grad = [fx(x(i+1),y(i+1),z(i+1)),fy(x(i+1),y(i+1),z(i+1)),fz(x(i+1),y(i+1),z(i+1))];
            grad = grad/sqrt(grad*transpose(grad));
            n1 = zeros(1,3);
            n1(1,mod(I,3)+1) = 1;
            if det([x(i+1),y(i+1),z(I+1);grad(1),grad(2),grad(3);n1(1),n1(2),n1(3)]) ~= 0
                normal = n1;
            else
                normal = zeros(1,3);
                normal(1,mod(I,3)+1) = 1;
            end
        end
        normal_new = true;
        edge_loop = false;
    else
        P_dash = (P_dash_array(1:i,1)-P_next(1)).^2 ...
            +(P_dash_array(1:i,2)-P_next(2)).^2+(P_dash_array(1:i,3)-P_next(3)).^2;
        [M,I] = min(P_dash);
        if M < ds^2
            [~,I_] = max(abs([x(i+1),y(i+1),z(i+1)]));
            [~,I_dash] = max(abs(P_dash_array(I,:)));
            if I_ == I_dash
                x(i+1) = P_dash_array(I,1);
                y(i+1) = P_dash_array(I,2);
                z(i+1) = P_dash_array(I,3);
                find_P_dash = false;
            end
        end
    end
end

%% 
n_points = 1000;
x = x(:,1:end);
y = y(:,1:end);
z = z(:,1:end);

x = [x;zeros(n_points,width(x))];
y = [y;zeros(n_points,width(y))];
z = [z;zeros(n_points,width(z))];
count_v_all = (n_points+1)*ones(1,width(x));
for i = (1:width(x))
    p_array = [x(1,i);y(1,i);z(1,i)];
    [~,I] = max(abs(p_array));
    dir = zeros(3,1);
    dir(I,1) = -sign(p_array(I));
    for j = (2:n_points+1)
        curv_vect = principal_directions(x(j-1,i),y(j-1,i),z(j-1,i));
        [~,I] = max(abs(transpose(dir)*curv_vect));
        dir = sign(transpose(dir)*curv_vect(:,I))*curv_vect(:,I);
        P = next_point_vector(x(j-1,i),y(j-1,i),z(j-1,i),f_init,ds,dir(1),dir(2),dir(3));
        x(j,i) = P(1);
        y(j,i) = P(2);
        z(j,i) = P(3);
        dir = [x(j,i)-x(j-1,i);y(j,i)-y(j-1,i);z(j,i)-z(j-1,i)];
        dir = dir/sqrt(transpose(dir)*dir);
         if abs(x(j,i)) > Size/2
            count_v_all(i) = j-1;
            x(j:end,i) = x(j-1,i);
            y(j:end,i) = y(j-1,i);
            z(j:end,i) = z(j-1,i);
            break
        end
        if abs(y(j,i)) > Size/2
            count_v_all(i) = j-1;
            x(j:end,i) = x(j-1,i);
            y(j:end,i) = y(j-1,i);
            z(j:end,i) = z(j-1,i);
            break
        end
        if abs(z(j,i)) > Size/2
            count_v_all(i) = j-1;
            x(j:end,i) = x(j-1,i);
            y(j:end,i) = y(j-1,i);
            z(j:end,i) = z(j-1,i);
            break
        end
    end
end
%%

is_div_line = zeros(1,width(x));
for i = 1:width(x)-1
    i_ = i+1;
    if and(count_v_all(i) == 1, count_v_all(i_) >= 10)
        is_div_line(i_) = 1; 
        continue
    end
    for j = (2:min(count_v_all(i),count_v_all(i_)))
        t1 = [x(j,i)-x(j-1,i);y(j,i)-y(j-1,i);z(j,i)-z(j-1,i)];
        t1 = t1/sqrt(transpose(t1)*t1);
        t2 = [x(j,i_)-x(j-1,i_);y(j,i_)-y(j-1,i_);z(j,i_)-z(j-1,i_)];
        t2 = t2/sqrt(transpose(t2)*t2);
        if transpose(t1)*t2 < 0.5
            is_div_line(i) = 1;
            is_div_line(i_) = 1;
            break
        end
    end
end

n_edge = is_div_line*ones(width(x),1);
div_indices = zeros(1,n_edge);
count_v = (n_points+1)*ones(1,n_edge);
x_edge = zeros(n_points+1,n_edge);
y_edge = zeros(n_points+1,n_edge);
z_edge = zeros(n_points+1,n_edge);
count = 0;
for i = 1:width(x)
    if is_div_line(i) == 1
        count = count+1;
        div_indices(1,count) = i;
        count_v(1,count) = count_v_all(1,i);
        x_edge(:,count) = x(:,i);
        y_edge(:,count) = y(:,i);
        z_edge(:,count) = z(:,i);
        if count_v(1,count) < 5
            count_v(1,count) = 1;
            count_v_all(1,i) = 1;
            x_edge(1:end,count) = x(1,i);
            y_edge(1:end,count) = y(1,i);
            z_edge(1:end,count) = z(1,i);
        end
    end
end

curv_max = ones(1,n_edge);
for i_ = 1:n_edge
    i = div_indices(1,i_);
    dt_max = 0;
    for j = 2:(count_v(i_)-1)
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
            curv_max(i_) = j;
        end
    end
end
%%

count = 1;
reaches_div_point = zeros(1,n_edge);
for i_ = 1:n_edge-1
    if curv_max(1,i_) > 2
        increment_count = false;
        i1 = div_indices(1,i_);
        j1 = curv_max(i_);
        x1 = x(j1,i1);
        y1 = y(j1,i1);
        z1 = z(j1,i1);
        for j_ = i_+1:n_edge
            if and(curv_max(1,j_) > 2, reaches_div_point(j_) == 0)
                i2 = div_indices(1,j_);
                j2 = curv_max(j_);
                x2 = x(j2,i2);
                y2 = y(j2,i2);
                z2 = z(j2,i2);
                if (x2-x1)^2+(y2-y1)^2+(z2-z1)^2 < (10*ds)^2
                    if and(reaches_div_point(i_) == 0, reaches_div_point(j_) == 0)
                        reaches_div_point(i_) = count;
                        reaches_div_point(j_) = count;
                        increment_count = true;
                    elseif reaches_div_point(i_)*reaches_div_point(j_) == 0
                        reaches_div_point(i_) = max(reaches_div_point(i_),reaches_div_point(j_));
                        reaches_div_point(j_) = reaches_div_point(i_);
                    else
                        reaches_div_point(i_) = min(reaches_div_point(i_),reaches_div_point(j_));
                        reaches_div_point(j_) = reaches_div_point(i_);
                    end
                    count_v(i_) = j1;
                    count_v_all(div_indices(i_)) = j1;
                    count_v(j_) = j2;
                    count_v_all(div_indices(j_)) = j2;
                    x_edge(j1:end,i_) = x1;
                    y_edge(j1:end,i_) = y1;
                    z_edge(j1:end,i_) = z1;
                    x_edge(j2:end,j_) = x2;
                    y_edge(j2:end,j_) = y2;
                    z_edge(j2:end,j_) = z2;
                end
            end
        end
        if increment_count
            count = count+1;
        end
    end
end

count = 0;
p = zeros(max(reaches_div_point),3);
for i = 1:n_edge
    if reaches_div_point(i) > count
        count = reaches_div_point(i);
        p(count,:) = find_divergence_point(x_edge(count_v(i),i),y_edge(count_v(i),i),...
            z_edge(count_v(i),i),f_init,ds);
    end
    if reaches_div_point(i)
        count_v(i) = count_v(i)+1;
        count_v_all(div_indices(i)) = count_v(i);
        x_edge(count_v(i):end,i) = p(reaches_div_point(i),1);
        y_edge(count_v(i):end,i) = p(reaches_div_point(i),2);
        z_edge(count_v(i):end,i) = p(reaches_div_point(i),3);
    end
end

clear count increment_count p i i_ j j_ i1 j1 x1 x2 y1 y2 z1 z2 t1 t2 dt dt_max
%%

intersect_index = zeros(n_points+1,n_edge);
for count = 1:max(reaches_div_point)
    n_ = 0;
    indices = zeros(1,n_edge);
    for i = 1:n_edge
        if reaches_div_point(i) == count
            n_ = n_+1;
            indices(1,n_) = i;
        end
    end
    indices = indices(1,1:n_);
    x_ = x_edge(:,indices);
    y_ = y_edge(:,indices);
    z_ = z_edge(:,indices);
    for i = 1:n_
        for j = count_v(indices(i))-2:-1:1
            I = is_intersecting(j,i,x_,y_,z_,x_(:,[1:i-1,i+1:n_]), y_(:,[1:i-1,i+1:n_]),...
               z_(:,[1:i-1,i+1:n_]),count_v(indices([1:i-1,i+1:n_]))-1,ds);
            if I
                if I(2) >= i
                    I(2) = I(2)+1;
                end
                intersect_index(j,indices(i)) = indices(I(2));
            end
        end
    end
    p = [x_edge(count_v(indices(1)),indices(1)),y_edge(count_v(indices(1)),indices(1)),...
        z_edge(count_v(indices(1)),indices(1))];
    for i = 1:n_
        i_ = indices(i);
        for j = count_v(i_)-1:-1:1
            if intersect_index(j,i_)
                k = intersect_index(j,i_);
                common_end = true;
                for j_ = count_v(k):-1:1
                    if intersect_index(j_,k) == i_
                        break
                    elseif intersect_index(j_,k)
                        common_end = false;
                    end
                end
                if and(common_end, i_ > k)
                    break
                end
                p1 = [x_edge(j,i_),y_edge(j,i_),z_edge(j,i_)];
                p2 = [x_edge(j+1,i_),y_edge(j+1,i_),z_edge(j+1,i_)];
                p12 = (p1+p2)/2;
                q1 = [x_edge(j_,k), y_edge(j_,k), z_edge(j_,k)];
                q2 = [x_edge(j_+1,k), y_edge(j_+1,k), z_edge(j_+1,k)];
                n = [fx(p12(1),p12(2),p12(3));fy(p12(1),p12(2),p12(3));...
                    fz(p12(1),p12(2),p12(3))];
                t_ = [(p2(1)-p1(1)),(p2(2)-p1(2)),(p2(3)-p1(3))];
                nt_ = [n(2)*t_(3)-n(3)*t_(2),n(3)*t_(1)-n(1)*t_(3),...
                    n(1)*t_(2)-n(2)*t_(1)];
                nt_ = nt_/sqrt(nt_*transpose(nt_));
                q_cap = ((q2-q1)/sqrt((q2-q1)*transpose(q2-q1)));
                q_plane = q1 + q_cap*((p12-q1)*transpose(nt_))/(q_cap*transpose(nt_));
                q_plane = gradiant_move(q_plane(1),q_plane(2),q_plane(3),f_init,ds);
                count_v(i_) = j+2;
                count_v_all(div_indices(i_)) = j+2;
                x_edge(j+1,i_) = q_plane(1);
                y_edge(j+1,i_) = q_plane(2);
                z_edge(j+1,i_) = q_plane(3);
                x_edge(j+2:end,i_) = p(1);
                y_edge(j+2:end,i_) = p(2);
                z_edge(j+2:end,i_) = p(3);
                x(:,div_indices(i_)) = x_edge(:,i_);
                y(:,div_indices(i_)) = y_edge(:,i_);
                z(:,div_indices(i_)) = z_edge(:,i_);
                if common_end
                    count_v(k) = j_+2;
                    count_v_all(div_indices(k)) = j_+2;
                    x_edge(j_+1,k) = q_plane(1);
                    y_edge(j_+1,k) = q_plane(2);
                    z_edge(j_+1,k) = q_plane(3);
                    x_edge(j_+2:end,k) = p(1);
                    y_edge(j_+2:end,k) = p(2);
                    z_edge(j_+2:end,k) = p(3);
                    x(:,div_indices(k)) = x_edge(:,k);
                    y(:,div_indices(k)) = y_edge(:,k);
                    z(:,div_indices(k)) = z_edge(:,k);
                end
                break
            end
        end
    end
    for i = 1:n_
        i_ = indices(i);
        if max(intersect_index(:,i_)) == 0
            x_ = x_edge(:,indices);
            y_ = y_edge(:,indices);
            z_ = z_edge(:,indices);
            j = count_v(i_)-1;      
            p1 = [x_(j,i),y_(j,i),z_(j,i)];
            p2 = [x_(j+1,i),y_(j+1,i),z_(j+1,i)];
            p12 = (p1+p2)/2;
            n = [fx(p12(1),p12(2),p12(3));fy(p12(1),p12(2),p12(3));...
                fz(p12(1),p12(2),p12(3))];
            t_ = [(p2(1)-p1(1)),(p2(2)-p1(2)),(p2(3)-p1(3))];
            nt_ = [n(2)*t_(3)-n(3)*t_(2),n(3)*t_(1)-n(1)*t_(3),...
                n(1)*t_(2)-n(2)*t_(1)];
            nt_ = nt_/sqrt(nt_*transpose(nt_));
            mid_plane = @(p_) (nt_(1)*(p_(1)-p12(1))...
                + nt_(2)*(p_(2)-p12(2)) + nt_(3)*(p_(3)-p12(3)));
            x_arr = abs(x_-p12(1));
            y_arr = abs(y_-p12(2));
            z_arr = abs(z_-p12(3));
            d_arr = (x_arr.^2+y_arr.^2+z_arr.^2);
            k = 0;
            for i__ = [1:i-1,i+1:n_]
                [~,I] =  min(d_arr(1:count_v(indices(i__))-1,i__));
                for j_ = max(I-2,1):min(I+2,count_v(indices(i__))-2)
                    p1_ = [x_(j_,i__),y_(j_,i__),z_(j_,i__)];
                    p2_ = [x_(j_+1,i__),y_(j_+1,i__),z_(j_+1,i__)];
                    p12_ = (p1_+p2_)/2;
                    if mid_plane(p1_)*mid_plane(p2_) <= 0
                        n = [fx(p12_(1),p12_(2),p12_(3)),...
                            fy(p12_(1),p12_(2),p12_(3)),fz(p12_(1),p12_(2),p12_(3))];
                        if det([(p1-p12_);(p1_-p2_);n])*det([(p2-p12_);(p1_-p2_);n]) <= 0
                            k = i__;
                            break
                        end
                    end
                end
                if k
                    [~,I] = min((x_(1:count_v(i_)-1,i)-x_(count_v(indices(k))-1,k)).^2 ...
                        +(y_(1:count_v(i_)-1,i)-y_(count_v(indices(k))-1,k)).^2 ...
                        +(z_(1:count_v(i_)-1,i)-z_(count_v(indices(k))-1,k)).^2);
                    count_v(i_) = I+1;
                    count_v_all(div_indices(i_)) = I+1;
                    x_edge(I+1:end,i_) = x_(count_v(indices(k))-1,k);
                    y_edge(I+1:end,i_) = y_(count_v(indices(k))-1,k);
                    z_edge(I+1:end,i_) = z_(count_v(indices(k))-1,k);
                    x(:,div_indices(i_)) = x_edge(:,i_);
                    y(:,div_indices(i_)) = y_edge(:,i_);
                    z(:,div_indices(i_)) = z_edge(:,i_);
                    break
                end
            end
        end
    end
end
%%    

for i = 1:n_edge
    if ~reaches_div_point(i)
        for j = 1:count_v(i)-1
            p1 = [x_edge(j,i),y_edge(j,i),z_edge(j,i)];
            p2 = [x_edge(j+1,i),y_edge(j+1,i),z_edge(j+1,i)];
            p12 = (p1+p2)/2;
            n_ = [fx(p12(1),p12(2),p12(3));fy(p12(1),p12(2),p12(3));...
                fz(p12(1),p12(2),p12(3))];
            t_ = [(p2(1)-p1(1)),(p2(2)-p1(2)),(p2(3)-p1(3))];
            nt_ = [n_(2)*t_(3)-n_(3)*t_(2),n_(3)*t_(1)-n_(1)*t_(3),...
                n_(1)*t_(2)-n_(2)*t_(1)];
            nt_ = nt_/sqrt(nt_*transpose(nt_));
            mid_plane = @(p_) (nt_(1)*(p_(1)-p12(1))...
                + nt_(2)*(p_(2)-p12(2)) + nt_(3)*(p_(3)-p12(3)));
            x_arr = abs(x_edge - p12(1));
            x_arr = [x_arr(:,1:i-1),x_arr(:,i+1:n_edge)];
            y_arr = abs(y_edge - p12(2));
            y_arr = [y_arr(:,1:i-1),y_arr(:,i+1:n_edge)];
            z_arr = abs(z_edge - p12(3));
            z_arr = [z_arr(:,1:i-1),z_arr(:,i+1:n_edge)];
            d_arr = sqrt(x_arr.^2+y_arr.^2+z_arr.^2);
            [M,~] =  min(d_arr,[],"all","linear");
            if M > ds
                continue
            end
            for i_ = [1:i-1,i+1:n_edge]
                x_arr = abs(x_edge(:,i_) - p12(1));
                y_arr = abs(y_edge(:,i_) - p12(2));
                z_arr = abs(z_edge(:,i_) - p12(3));
                d_arr = sqrt(x_arr.^2+y_arr.^2+z_arr.^2);
                [M,I] =  min(d_arr,[],"all","linear");
                if M > ds
                    continue
                end
                for j_ = max(I-2,1):min(I+2,count_v(i_)-1)
                    p1_ = [x_edge(j_,i_),y_edge(j_,i_),z_edge(j_,i_)];
                    p2_ = [x_edge(j_+1,i_),y_edge(j_+1,i_),z_edge(j_+1,i_)];
                    p12_ = (p1_+p2_)/2;
                    cross_p = [(p2(2)-p1(2))*(p2_(3)-p1_(3)) - (p2(3)-p1(3))*(p2_(2)-p1_(2)),...
                        (p2(3)-p1(3))*(p2_(1)-p1_(1)) - (p2(1)-p1(1))*(p2_(3)-p1_(3)),...
                        (p2(1)-p1(1))*(p2_(2)-p1_(2)) - (p2(2)-p1(2))*(p2_(1)-p1_(1))];
                    cross_p = cross_p/sqrt((p2-p1)*transpose(p2-p1)*(p2_-p1_)*transpose(p2_-p1_));
                    if and(sqrt((p12_-p12)*transpose(p12_-p12)) < ds, sqrt(cross_p*transpose(cross_p)) > 0.1)
                        if mid_plane(p1_)*mid_plane(p2_) <= 0
                            n_ = [fx(p12_(1),p12_(2),p12_(3)),...
                                fy(p12_(1),p12_(2),p12_(3)),fz(p12_(1),p12_(2),p12_(3))];
                            if det([(p1-p12_);(p1_-p2_);n_])*det([(p2-p12_);(p1_-p2_);n_]) <= 0          
                                intersect_index(j,i) = i_;
                                break
                            end
                        end
                    end
                end
            end
        end
    end
end
%%
%Run from here after loading variables

break_indices = zeros(2,n_edge);
for i = 1:n_edge
    if ~reaches_div_point(i)
        for j = 1:count_v(i)
            k = intersect_index(j,i);
            if k
                if reaches_div_point(1,k)
                    for j_ = j+1:count_v(i)
                        k_ = intersect_index(j_,i);
                        if k_
                            for j__ = 1:count_v(k_)
                                if intersect_index(j__,k_) == i
                                    intersect_index(j__,k_) = 0;
                                end
                            end
                        end
                    end
                    break_indices(1,i) = j;
                    break_indices(2,i) = k;
                    break
                end
            end
        end
    end
end

for i = 1:n_edge
    if and(break_indices(1,i),break_indices(2,i))
        j = break_indices(1,i);
        k = break_indices(2,i);
        for j_ = j-1:-1:1
            if intersect_index(j_,i)
                j = j_;
                k = intersect_index(j,i);
                break
            end
        end
        [~,I] = min(sqrt((x_edge(1:count_v(k),k)-x_edge(j,i)).^2 ...
            + (y_edge(1:count_v(k),k)-y_edge(j,i)).^2 + (z_edge(1:count_v(k),k)-z_edge(j,i)).^2));
        q1 = [x_edge(I,k), y_edge(I,k), z_edge(I,k)];
        p1 = [x_edge(j,i),y_edge(j,i),z_edge(j,i)];
        p2 = [x_edge(j+1,i),y_edge(j+1,i),z_edge(j+1,i)];
        p12 = (p1+p2)/2;
        n_ = [fx(p12(1),p12(2),p12(3));fy(p12(1),p12(2),p12(3));...
            fz(p12(1),p12(2),p12(3))];
        t_ = [(p2(1)-p1(1)),(p2(2)-p1(2)),(p2(3)-p1(3))];
        nt_ = [n_(2)*t_(3)-n_(3)*t_(2),n_(3)*t_(1)-n_(1)*t_(3),...
            n_(1)*t_(2)-n_(2)*t_(1)];
        nt_ = nt_/sqrt(nt_*transpose(nt_));
        if I > 1
            if (p2-q1)*transpose([x_edge(I+1,k), y_edge(I+1,k), z_edge(I+1,k)]-q1)...
                    >= (p2-q1)*transpose([x_edge(I-1,k), y_edge(I-1,k), z_edge(I-1,k)]-q1)
                q2 = [x_edge(I+1,k), y_edge(I+1,k), z_edge(I+1,k)];
            else
                q2 = [x_edge(I-1,k), y_edge(I-1,k), z_edge(I-1,k)];
            end
        else
            q2 = [x_edge(I+1,k), y_edge(I+1,k), z_edge(I+1,k)];
        end
        q_cap = ((q2-q1)/sqrt((q2-q1)*transpose(q2-q1)));
        q_plane = q1 + q_cap*((p12-q1)*transpose(nt_))/(q_cap*transpose(nt_));
        q_plane = gradiant_move(q_plane(1),q_plane(2),q_plane(3),f_init,ds);
        x_edge(j+1:end,i) = q_plane(1);
        y_edge(j+1:end,i) = q_plane(2);
        z_edge(j+1:end,i) = q_plane(3);
        x(j+1:end,div_indices(i)) = q_plane(1);
        y(j+1:end,div_indices(i)) = q_plane(2);
        z(j+1:end,div_indices(i)) = q_plane(3);
        count_v(i) = j+1;
        count_v_all(div_indices(i)) = j+1;
        for j_ = j+1:count_v(i)
            k_ = intersect_index(j_,i);
            if k_
                for j__ = 1:count_v(k_)
                    if intersect_index(j__,k_) == i
                        intersect_index(j__,k_) = 0;
                    end
                end
            end
        end
    end
end

%%
frac = 20;
for i = 1:frac:width(x)
    if ~is_div_line(i)
        for j = 1:count_v_all(i)-1
            I = is_intersecting(j,i,x,y,z,x_edge,y_edge,z_edge,count_v,ds);
            if I
                q1 = [x_edge(I(1),I(2)), y_edge(I(1),I(2)), z_edge(I(1),I(2))];
                p1 = [x(j,i),y(j,i),z(j,i)];
                p2 = [x(j+1,i),y(j+1,i),z(j+1,i)];
                p12 = (p1+p2)/2;
                n_ = [fx(p12(1),p12(2),p12(3));fy(p12(1),p12(2),p12(3));...
                    fz(p12(1),p12(2),p12(3))];
                t_ = [(p2(1)-p1(1)),(p2(2)-p1(2)),(p2(3)-p1(3))];
                nt_ = [n_(2)*t_(3)-n_(3)*t_(2),n_(3)*t_(1)-n_(1)*t_(3),...
                    n_(1)*t_(2)-n_(2)*t_(1)];
                nt_ = nt_/sqrt(nt_*transpose(nt_));
                if I > 1
                    if (p2-q1)*transpose([x_edge(I(1)+1,I(2)), y_edge(I(1)+1,I(2)), z_edge(I(1)+1,I(2))]-q1)...
                            >= (p2-q1)*transpose([x_edge(I(1)-1,I(2)), y_edge(I(1)-1,I(2)), z_edge(I(1)-1,I(2))]-q1)
                        q2 = [x_edge(I(1)+1,I(2)), y_edge(I(1)+1,I(2)), z_edge(I(1)+1,I(2))];
                    else
                        q2 = [x_edge(I(1)-1,I(2)), y_edge(I(1)-1,I(2)), z_edge(I(1)-1,I(2))];
                    end
                else
                    q2 = [x_edge(I(1)+1,I(2)), y_edge(I(1)+1,I(2)), z_edge(I(1)+1,I(2))];
                end
                q_cap = ((q2-q1)/sqrt((q2-q1)*transpose(q2-q1)));
                q_plane = q1 + q_cap*((p12-q1)*transpose(nt_))/(q_cap*transpose(nt_));
                q_plane = gradiant_move(q_plane(1),q_plane(2),q_plane(3),f_init,ds);
                x(j+1:end,i) = q_plane(1);
                y(j+1:end,i) = q_plane(2);
                z(j+1:end,i) = q_plane(3);
                count_v_all(i) = j+1;
                break
            end
        end
    end
end
%%

for i = 1:frac:width(x)
    if and(~is_div_line(1,i), count_v_all(i) > 1)
        x_ = x(count_v_all(i),i);
        y_ = y(count_v_all(i),i);
        z_ = z(count_v_all(i),i);
        if min(abs(abs([x_,y_,z_])-Size/2)) <= ds
            lin_indices = 1 + height(x)*(([1:frac:i-frac,i+frac:frac:width(x)])-1);
            if ~isempty(lin_indices)
                [M,I] = min((x(lin_indices)-x_).^2 + (y(lin_indices)-y_).^2 + (z(lin_indices)-z_).^2);
                if I < (i-1)/frac + 1
                    I = (I-1)*frac + 1;
                else
                    I = I*frac + 1;
                end
                if is_div_line(I) || count_v_all(I) == 1 || sqrt(M) > frac*ds
                    continue
                else
                    [~,M] = max(abs([x(1,i),y(1,i),z(1,i)]));
                    [~,M_] = max(abs([x(count_v_all(I),I),y(count_v_all(I),I),z(count_v_all(I),I)]));
                    if M ~= M_
                        continue
                    end
                end
                old_count = count_v_all(i);
                count_v_all(i) = max(count_v_all(i),count_v_all(I));
                x_new = zeros(n_points+1,1);
                y_new = zeros(n_points+1,1);
                z_new = zeros(n_points+1,1);
                for j = 1:count_v_all(i)
                    frc = (j-1)/(count_v_all(i)-1);
                    j_ = [floor(frc*(old_count-1))+1,floor(frc*(old_count-1))+2,...
                        floor((1-frc)*(count_v_all(I)-1))+1,floor((1-frc)*(count_v_all(I)-1))+2];
                    a = [j_(2)-(frc*(old_count-1)+1),(frc*(old_count-1)+1)-j_(1),...
                        j_(4)-((1-frc)*(count_v_all(I)-1)+1),((1-frc)*(count_v_all(I)-1))+1-j_(3)];
                    if j_(2) > old_count
                        j_(2) = old_count;
                    end
                    if j_(4) > count_v_all(I)
                        j_(4) = count_v_all(I);
                    end
                    x_new(j,1) = (1-frc)*(a(1)*x(j_(1),i)+a(2)*x(j_(2),i))...
                        + frc*(a(3)*x(j_(3),I)+a(4)*x(j_(4),I));
                    y_new(j,1) = (1-frc)*(a(1)*y(j_(1),i)+a(2)*y(j_(2),i))...
                        + frc*(a(3)*y(j_(3),I)+a(4)*y(j_(4),I));
                    z_new(j,1) = (1-frc)*(a(1)*z(j_(1),i)+a(2)*z(j_(2),i))...
                        + frc*(a(3)*z(j_(3),I)+a(4)*z(j_(4),I));
                    p = gradiant_move(x_new(j,1),y_new(j,1),z_new(j,1),f_init,ds);
                    x_new(j,1) = p(1);
                    y_new(j,1) = p(2);
                    z_new(j,1) = p(3);
                end
                x_new(count_v_all(i):end,1) = x_new(count_v_all(i),1);
                x(:,i) = x_new(:,1);
                y_new(count_v_all(i):end,1) = y_new(count_v_all(i),1);
                y(:,i) = y_new(:,1);
                z_new(count_v_all(i):end,1) = z_new(count_v_all(i),1);
                z(:,i) = z_new(:,1);
                x(1:end,I) = x(1,I);
                y(1:end,I) = y(1,I);
                z(1:end,I) = z(1,I);
                count_v_all(I) = 1;
            end
        end
    end
end
clear lin_indices M M_ I a j_ fra old_count p x_new y_new z_new

%%
close all;

n = 100;
[X, Y, Z] = meshgrid(-Size/2:Size/n:Size/2);
OBJ = cos(2*pi*X/Size).* sin(2*pi*Y/Size) + cos(2*pi*Y/Size).* sin(2*pi*Z/Size)...
        + cos(2*pi*Z/Size).* sin(2*pi*X/Size);

figure(1)
isosurface(X,Y,Z, (OBJ), 0); 

xlim([-0.5*Size 0.5*Size])
ylim([-0.5*Size 0.5*Size])
zlim([-0.5*Size 0.5*Size])
hold on;
plot3(x(1,:),y(1,:),z(1,:),"-","Color","b")
for i = 1:frac:width(x)
    if ~is_div_line(i)
        plot3(x(1:1:end,i), y(1:1:end,i), z(1:1:end,i), "-", "Color", "k")
    end
end
for i_ = 1:n_edge
    if reaches_div_point(i_)
        color = [1 0.5 0];
    else
        color = "r";
    end
    plot3(x_edge(:,i_),y_edge(:,i_),z_edge(:,i_),"-","Color",color)
end

for i = 1:n_edge
    for j = 1:count_v(i)
        if intersect_index(j,i)
            plot3(x_edge(j,i),y_edge(j,i),z_edge(j,i),"o","Color","g")
        end
    end
end

%%
find_divergence_point(0.250105,-0.253822,0.246284,0,ds)

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
    epsilon_absolute = 1e-18;
    p1 = [x,y,z];
    p1 = gradiant_move(p1(1),p1(2),p1(3),const,ds);
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
        + nt_(2)*(p_(2)-p12(2)) + nt_(3)*(p_(3)-p12(3)));
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
            p1_ = [x(j_,i_),y(j_,i_),z(j_,i_)];
            p2_ = [x(j_+1,i_),y(j_+1,i_),z(j_+1,i_)];
            p12_ = (p1_+p2_)/2;
            if sqrt((p12_-p12)*transpose(p12_-p12)) < ds
                if mid_plane(p1_)*mid_plane(p2_) <= 0
                    n_ = [fx(p12_(1),p12_(2),p12_(3)),...
                        fy(p12_(1),p12_(2),p12_(3)),fz(p12_(1),p12_(2),p12_(3))];
                    if det([(p1-p12_);(p1_-p2_);n_])*det([(p2-p12_);(p1_-p2_);n_]) <= 0
                        value = [j_,i_];
                        return
                    end
                end
            end
        end
    end
    value = false;
end
                                         