function g = grav(GM, a, b, c, x, y, z)
    
    den_g = [a b c] + [x y z];
    g(1) = -(GM*x)/(norm(den_g)^3);
    g(2) = -(GM*y)/(norm(den_g)^3);
    g(3) = -(GM*z)/(norm(den_g)^3);
    
%     g = (GM*[x y z])./(vecnorm(den_g).^3);
end

