clc;
clear;

Size = 1;
levels = 0;
n_points = 15000;
ds = 0.001;

f_init = levels(1);
P_start= zeros(3,3);
for k = 1:3
    normal = zeros(1,3);
    tangent = zeros(1,3);
    normal(1,mod(k-1,3)+1) = 1;
    tangent(1,mod(k,3)+1) = 1;
    P_start(k,:) = P_start(k,:)+(Size/2)*normal+ds*tangent;
    P_start(k,:) = find_level(P_start(k,1),P_start(k,2),P_start(k,3),...
        f_init,ds,tangent(1),tangent(2),tangent(3));
end

clear normal tangent

x_all = zeros(6,n_points);
y_all = zeros(6,n_points);
z_all = zeros(6,n_points);
loop_length = zeros(1,6);
face_indices = ones(1,7);
Px = @(x,y,z) (x-Size/2)*(x+Size/2);
Py = @(x,y,z) (y-Size/2)*(y+Size/2);
Pz = @(x,y,z) (z-Size/2)*(z+Size/2);

for k = 1:3
    normal = zeros(1,3);
    normal(1,mod(k-1,3)+1) = 1;
    x = zeros(1,n_points);
    x(1) = P_start(k,1);
    y = zeros(1,n_points);
    y(1) = P_start(k,2);
    z = zeros(1,n_points);
    z(1) = P_start(k,3);
    dir = 1;
    for i = 1:width(x)-1
        P_next = next_point_plane(x(i),y(i),z(i),f_init,ds,normal(1),normal(2),normal(3),dir);
        x(i+1) = P_next(1);
        y(i+1) = P_next(2);
        z(i+1) = P_next(3);
        if and(i > 2, (x(i+1)-x(1))^2+(y(i+1)-y(1))^2+(z(i+1)-z(1))^2 <= ds^2)
            x(i+1) = x(1);
            x_3D(k,:) = x;
            y(i+1) = y(1);
            y_3D(k,:) = y;
            z(i+1) = z(1);
            z_3D(k,:) = z;
            loop_length(1,k) = i+1;
            loop_length(1,k+3) = i+1;
            break
        end
        if abs(x(i+1)) > Size/2
            x(i+1) = Size/2*sign(x(i+1));
            break
        end
        if abs(y(i+1)) > Size/2
            y(i+1) = Size/2*sign(y(i+1));
            break
        end
        if abs(z(i+1)) > Size/2
            z(i+1) = Size/2*sign(z(i+1));
            break
        end
    end
end

for i = 2:7
    face_indices(1,i) = face_indices(1,i-1) + loop_length(1,i-1);
end

x = zeros(1,loop_length*ones(6,1));
y = zeros(1,loop_length*ones(6,1));
z = zeros(1,loop_length*ones(6,1));
for k = 1:3
    x(face_indices(k):face_indices(k+1)-1) = x_3D(k,1:loop_length(1,k));
    y(face_indices(k):face_indices(k+1)-1) = y_3D(k,1:loop_length(1,k));
    z(face_indices(k):face_indices(k+1)-1) = z_3D(k,1:loop_length(1,k));
    normal = zeros(1,3);
    normal(1,mod(k-1,3)+1) = 1;
    x(face_indices(k+3):face_indices(k+4)-1) = x_3D(k,1:loop_length(1,k))-Size*normal(1);
    y(face_indices(k+3):face_indices(k+4)-1) = y_3D(k,1:loop_length(1,k))-Size*normal(2);
    z(face_indices(k+3):face_indices(k+4)-1) = z_3D(k,1:loop_length(1,k))-Size*normal(3);
end

clear x_3D y_3D z_3D count loop_length normal tangent
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

is_div_line = zeros(1,width(x));
for i = 1:width(x)-1
    i_ = i+1;
    k = 0;
    for j = 1:6
        if i_ == face_indices(1,j)
            k = 1;
            break
        end
    end
    if k
        continue
    end
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

%% 
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
        x(:,div_indices(i)) = x_edge(:,i);
        y(:,div_indices(i)) = y_edge(:,i);
        z(:,div_indices(i)) = z_edge(:,i);
    end
end

for i = 2:n_edge-1
    if and(div_indices(i) == div_indices(i-1)+1, div_indices(i) == div_indices(i+1)-1)
        reaches_div_point(i) = reaches_div_point(i-1);
        [~,I] = min((x_edge(1:count_v(i),i)-p(reaches_div_point(i),1)).^2 ...
            +(y_edge(1:count_v(i),i)-p(reaches_div_point(i),2)).^2 ...
            +(z_edge(1:count_v(i),i)-p(reaches_div_point(i),3)).^2);
        curv_max(i) = I;
        count_v(i) = I+1;
        count_v_all(div_indices(i)) = count_v(i);
        x_edge(count_v(i):end,i) = p(reaches_div_point(i),1);
        y_edge(count_v(i):end,i) = p(reaches_div_point(i),2);
        z_edge(count_v(i):end,i) = p(reaches_div_point(i),3);
        x(:,div_indices(i)) = x_edge(:,i);
        y(:,div_indices(i)) = y_edge(:,i);
        z(:,div_indices(i)) = z_edge(:,i);
    end
end

clear count increment_count p i i_ j j_ i1 j1 k I x1 x2 y1 y2 z1 z2 t1 t2 dt dt_max
%%

intersect_index = zeros(max(count_v),n_edge);
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
            k_ = zeros(1,n_);
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
                            k_(1,i__) = i__;
                            break
                        end
                    end
                end
            end
            d_min = 0;
            for i__ = 1:n_
                if and(i__ ~= i, k_(1,i__))
                    d_dash = (x_(count_v(indices(i__))-1,i__) - x_(count_v(i_)-1,i))^2 ...
                        + (y_(count_v(indices(i__))-1,i__) - y_(count_v(i_)-1,i))^2 ...
                        + (z_(count_v(indices(i__))-1,i__) - z_(count_v(i_)-1,i))^2;
                    if d_min == 0
                        d_min = d_dash;
                        k = i__;
                    elseif d_dash < d_min
                        d_min = d_dash;
                        k = i__;
                    end
                end
            end
            if d_min
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
            end
            clear d_min d_dash k_
        end
    end
    clear indices n_
end
%%

frac = 20;
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
                elseif max(abs([x(count_v_all(I),I),y(count_v_all(I),I),z(count_v_all(I),I)])) == Size/2
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
clear lin_indices M M_ I a j_ frc old_count p x_new y_new z_new
%%
%Run from here after loading variables

for i = 1:frac:width(x)
    i_ = count_v_all(i);
    p = [x(i_,i),y(i_,i),z(i_,i)];
    M = max(abs(p));
    if is_div_line(i) || count_v_all(i) == 1
        continue
    elseif min(Size/2 - M) > ds || M == Size/2
        continue
    end
    M = Size;
    for j = [1:frac:i-frac,i+frac:frac:width(x)]
        j_ = count_v_all(j);
        if is_div_line(j) || j_ == 1
            continue
        end
        p1 = [x(1,j),y(1,j),z(1,j)];
        p2 = [x(j_,j),y(j_,j),z(j_,j)];
        if max(abs(p2)) ~= Size/2
            continue
        else
            [~,M_1] = max(abs(p1));
            [~,M_2] = max(abs(p2));
            [~,M_] = max(abs(p));
            if M_ ~= M_1 && M_ ~= M_2
                continue
            end
            [~,M_] = max(abs([x(1,i),y(1,i),z(1,i)]));
            if M_ ~= M_1 && M_ ~= M_2
                continue
            end
        end
        if (p1-p)*transpose(p1-p) < M^2
            M = sqrt((p1-p)*transpose(p1-p));
            dir = 0;
            I = j;
        end
        if (p2-p)*transpose(p2-p) < M^2
            M = sqrt((p2-p)*transpose(p2-p));
            dir = 1;
            I = j;
        end
    end
    if M > frac*ds
        continue
    end
    x_new = zeros(n_points+1,1);
    y_new = zeros(n_points+1,1);
    z_new = zeros(n_points+1,1);
    I_ = count_v_all(I);
    for j = 2:i_
        [M,J_] = min((x(2:I_-1,I)-x(j,i)).^2 + (y(2:I_-1,I)-y(j,i)).^2 ...
            + (z(2:I_-1,I)-z(j,i)).^2);
        J_ = J_+1; 
        if M < (frac*7*ds/8)^2
            J = j;
            break
        end
    end
    if dir == 0
        count_v_all(i) = max(J,I_-J_+1);
    else
        count_v_all(i) = max(J,J_);
    end
    for j = 1:count_v_all(i)
        frc = (j-1)/(count_v_all(i)-1);
        if dir == 0
            j_ = [floor(frc*(J-1))+1,floor(frc*(J-1))+2,...
                I_-floor(frc*(I_-J_)),I_-floor(frc*(I_-J_))-1];
            a = [j_(2)-(frc*(J-1)+1),(frc*(J-1)+1)-j_(1),...
                (I_-j_(4))-frc*(I_-J_),(frc*(I_-J_))-(I_-j_(3))];
        elseif dir == 1
            j_ = [floor(frc*(J-1))+1,floor(frc*(J-1))+2,...
                floor(frc*(J_-1))+1,floor(frc*(J_-1))+2];
            a = [j_(2)-(frc*(J-1)+1),(frc*(J-1)+1)-j_(1),...
                j_(4)-(frc*(J_-1)+1),(frc*(J_-1))+1-j_(3)];
        end
        if j_(2) > i_
            j_(2) = i_;
        elseif j_(2) == 0
            j_(2) = 1;
        end
        if j_(4) > I_
            j_(4) = I_;
        elseif j_(4) == 0
            j_(4) = 1;
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
end
clear a dir frc I I_ i_ J J_ j_ M M_ M_1 M_2 p p1 p2 x_new y_new z_new
%%

count_v_all(div_indices) = count_v;
x(:,div_indices) = x_edge;
y(:,div_indices) = y_edge;
z(:,div_indices) = z_edge;
x = x(:,1:frac:width(x));
y = y(:,1:frac:width(y));
z = z(:,1:frac:width(z));
count_v_all = count_v_all(:,1:frac:width(count_v_all));

indices = zeros(n_points+1,width(x));
indices(1:count_v_all(1),1) = transpose(1:count_v_all(1));
for i = 2:width(x)
    indices(1:count_v_all(i),i) = indices(count_v_all(i-1),i-1)+transpose(1:count_v_all(i));
end

for i = 1:width(x)
    I = count_v_all(i);
    if I == 1
        for j = [1:i-1,i+1:width(x)]
            J = count_v_all(j);
            if x(J,j) == x(1,i) && y(J,j) == y(1,i) && z(J,j) == z(1,i)
                count_v_all(i) = J;
                indices(1:J,i) = indices(J:-1:1,j);
                x(1:J,i) = x(J:-1:1,j);
                y(1:J,i) = y(J:-1:1,j);
                z(1:J,i) = z(J:-1:1,j);
                break
            end
        end
    elseif max(abs([x(I,i),y(I,i),z(I,i)])) ~= Size/2 && ~is_div_line(1+frac*(i-1))
        for j = [1:i-1,i+1:width(x)]
            J = count_v_all(j);    
            if ~((x(1,j)-x(1,i))^2+(y(1,j)-y(1,i))^2+(z(1,j)-z(1,i))^2 < (frac*ds)^2) &&...
                    ~((x(J,j)-x(1,i))^2+(y(J,j)-y(1,i))^2+(z(J,j)-z(1,i))^2 < (frac*ds)^2)
                continue
            end
            [M,J_] = min((x(1:J,j)-x(I,i)).^2+(y(1:J,j)-y(I,i)).^2+(z(1:J,j)-z(I,i)).^2);
            if M == 0
                dir = (x(1,j)-x(1,i))^2+(y(1,j)-y(1,i))^2+(z(1,j)-z(1,i))^2 ...
                    < (x(J,j)-x(1,i))^2+(y(J,j)-y(1,i))^2+(z(J,j)-z(1,i))^2;
                if dir
                    x(I+1:I+J-J_,i) = x(J_+1:J,j);
                    y(I+1:I+J-J_,i) = y(J_+1:J,j);
                    z(I+1:I+J-J_,i) = z(J_+1:J,j);
                    count_v_all(i) = count_v_all(i)+J-J_;
                    indices(I:count_v_all(i),i) = indices(J_:J,j);
                else
                    x(I+1:I+J_-1,i) = x(J_-1:-1:1,j);
                    y(I+1:I+J_-1,i) = y(J_-1:-1:1,j);
                    z(I+1:I+J_-1,i) = z(J_-1:-1:1,j);
                    count_v_all(i) = count_v_all(i)+J_-1;
                    indices(I:count_v_all(i),i) = indices(J_:-1:1,j);
                end
                break
            end
        end
    end
end
clear I J J_ M

p_array = zeros(max(indices,[],"all"),3);
for i = 1:width(indices)
    for j = 1:count_v_all(i)
        p_array(indices(j,i),:) = [x(j,i),y(j,i),z(j,i)];
    end
end

count = 1;
[~,I_old] = max(abs([x(1,1),y(1,1),z(1,1)]));
for i = 1:width(x)
    [~,I_new] = max(abs([x(1,i),y(1,i),z(1,i)]));
    if I_new ~= I_old
        count = count+1;
        face_indices(count) = i;
    end
    I_old = I_new;
end
face_indices(count+1) = width(x)+1;
clear count I_new I_old

%%
close all;
figure(1)
view(3)
hold on;
for i = 1:6
    m = [face_indices(i):face_indices(i+1)-1,face_indices(i)];
    for k = 1:width(m)-1
        if is_div_line((face_indices(i)+k-2)*frac+1)
            m = [m(1:k-1),m(k+1:end)];
        end
    end
    for j = 1:frac:min(h)
        array = indices(j,m);
        for k = 1:width(m)-1
            for j_ = max(j-frac,1):j+frac
                if array(k) == indices(j_,m(k+1))
                    array(k) = min(array(k),array(k+1));
                    array(k+1) = array(k);
                end
            end
        end
        plot3(p_array(array,1),p_array(array,2),p_array(array,3),"-o","MarkerSize",4,"Color","b")
    end
    array = zeros(1,width(m));
    for j = 1:width(m)
        array(1,j) = indices(floor(count_v_all(m(j))/2),m(j));
    end
    plot3(p_array(array,1),p_array(array,2),p_array(array,3),"-o","MarkerSize",4,"Color","b")
end
clear array
for i = 1:width(x)
    if ~is_div_line(1+frac*(i-1))
        plot3(x(1:count_v_all(i),i), y(1:count_v_all(i),i),...
            z(1:count_v_all(i),i), "-", "Color", "k")
        if count_v_all(i) > 1
            scatter3(x(1,i),y(1,i),z(1,i),36,[0 0 1])
        end
    end
end

%%
n_levels = width(levels);
x_3D = zeros([size(x),n_levels]);
y_3D = zeros([size(y),n_levels]);
z_3D = zeros([size(z),n_levels]); 
x_3D(:,:,1) = x;
y_3D(:,:,1) = y;
z_3D(:,:,1) = z;
for k = 2:n_levels
    level = levels(k);
    for i = 1:width(x)
        for j = 1:count_v_all(i)
            p = gradiant_move(x_3D(j,i,k-1),y_3D(j,i,k-1),z_3D(j,i,k-1),level,0.1);
            x_3D(j,i,k) = p(1);
            y_3D(j,i,k) = p(2);
            z_3D(j,i,k) = p(3);
        end
    end
end
%%

n = 100;
[X, Y, Z] = meshgrid(-Size/2:Size/n:Size/2);
OBJ = cos(2*pi*X/Size) + cos(2*pi*Y/Size) + cos(2*pi*Z/Size);

figure(2)
for i = 1:n_levels
    isosurface(X,Y,Z, (OBJ), levels(i)); 
    hold on;
end

xlim([-0.5*Size 0.5*Size])
ylim([-0.5*Size 0.5*Size])
zlim([-0.5*Size 0.5*Size])

for k = 1:n_levels
    for i = 1:6
        for j = 1:10
            plot3(x_3D(j,[face_indices(i):face_indices(i+1)-1,face_indices(i)],k),...
                y_3D(j,[face_indices(i):face_indices(i+1)-1,face_indices(i)],k),...
                z_3D(j,[face_indices(i):face_indices(i+1)-1,face_indices(i)],k),"-","Color","b")
        end
    end
    for i = 1:width(x)
        if ~is_div_line(1+frac*(i-1))
            plot3(x_3D(1:count_v_all(i),i,k), y_3D(1:count_v_all(i),i,k),...
                z_3D(1:count_v_all(i),i,k), "-", "Color", "k")
            if count_v_all(i) > 1
                scatter3(x_3D(1,i,k),y_3D(1,i,k),z_3D(1,i,k),36,[0 0 1])
                if max(abs([x(count_v_all(i),i),y(count_v_all(i),i),z(count_v_all(i),i)])) ~= Size/2 && ~is_div_line(i)
                    scatter3(x_3D(count_v_all(i),i,k),y_3D(count_v_all(i),i,k),z_3D(count_v_all(i),i,k),36,[1 0 0])
                end
            end
        end
    end
    % for i = 1:n_edge
    %     i_ = div_indices(i);
    %     plot3(x_3D(1:count_v(i),i_,k), y_3D(1:count_v(i),i_,k),z_3D(1:count_v(i),i_,k), "-", "Color", "r")
    %     j = curv_max(i);
    %     scatter3(x_3D(j,i_,k),y_3D(j,i_,k),z_3D(j,i_,k),30+reaches_div_point(i),[0 1 0])
    % end
end
%%

function value = f(x,y,z)
    value = cos(2*pi*x)+cos(2*pi*y)+cos(2*pi*z);
end

function value = fx(x,~,~) 
    value = -2*pi*sin(2*pi*x);
end

function value = fy(~,y,~) 
    value = -2*pi*sin(2*pi*y);
end

function value = fz(~,~,z) 
    value = -2*pi*sin(2*pi*z);
end

function value = fxx(x,~,~) 
    value = -4*pi^2*cos(2*pi*x);
end

function value = fyy(~,y,~) 
    value = -4*pi^2*cos(2*pi*y);
end

function value = fzz(~,~,z) 
    value = -4*pi^2*cos(2*pi*z);
end

function value = fxy(~,~,~) 
    value = 0;
end

function value = fyz(~,~,~) 
    value = 0;
end

function value = fzx(~,~,~) 
    value = 0;
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
    epsilon = 1e-9;
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
    epsilon_absolute = 1e-9;
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
        p1 = p1 - 0.1*H1*g_/(g_*transpose(g_));
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