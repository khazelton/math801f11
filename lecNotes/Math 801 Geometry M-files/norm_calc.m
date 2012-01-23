function [Np] = norm_calc(u,v,x1,x2,x3)

%% computing the differential
z1_x1 = diff(x1,u);
z2_x1 = diff(x1,v);

z1_x2 = diff(x2,u);
z2_x2 = diff(x2,v);

z1_x3 = diff(x3,u);
z2_x3 = diff(x3,v);

%% evaluating the values for differnt values.For this example I have taken a limited range of u and v to check for the correctness of the code.

k = 1;
for u = 1: 5
    for v = 1:5
        a_x = eval(z1_x1);
        b_x = eval(z1_x2);
        c_x = eval(z1_x3);
        
        final_x(k,:) = [a_x b_x c_x];

        a_y = eval(z2_x1);
        b_y = eval(z2_x2);
        c_y = eval(z2_x3);
        
        final_y(k,:) = [a_y b_y c_y];
        
        k = k+1;
    end
end

%%% computing the cross product.

 cross_val = cross(final_x,final_y);

 %% computing the dot product
 dot_val1 = dot(final_x,final_y);
 dot_val2 = dot(final_y,final_x);
 
 %% calculating the norm
 norm_cross_val = norm(cross_val);
 
 %% the final value
 Np = cross_val/norm_cross_val;
 
 return
 