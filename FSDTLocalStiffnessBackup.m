function [stiffness] = FSDTLocalStiffnessBackup(coords)
    syms z xsi eta B3c B3f;
    
    %% MATERIAL PROPERTIES
    thickness = 0.01; % m
    E = 10920;
    nu = 0.3;
    G = E / (2 * (nu + 1));
    alpha = 5/6;
        
    N1 = -(1-eta)*(1-xsi)*(1+xsi+eta)/4;
    N2 = -(1-eta)*(1+xsi)*(1-xsi+eta)/4;
    N3 = -(1+eta)*(1+xsi)*(1-xsi-eta)/4;
    N4 = -(1+eta)*(1-xsi)*(1+xsi-eta)/4;

    N5 = (1-eta)*(1-xsi)*(1+xsi)/2;
    N6 = (1-eta)*(1+xsi)*(1+eta)/2;
    N7 = (1+eta)*(1-xsi)*(1+xsi)/2;
    N8 = (1-eta)*(1-xsi)*(1+eta)/2;
    
    Ni = [N1, N2, N3, N4, N5, N6, N7, N8];
    
    x = sum(Ni * coords(:, 1));
    y = sum(Ni * coords(:, 2));
    
    xdifxi = double(subs(subs(sum(diff(Ni, xsi) * coords(:, 1)), xsi, 1/sqrt(3)), eta, 1/sqrt(3)));
    xdifeta = double(subs(subs(sum(diff(Ni, eta) * coords(:, 1)), xsi, 1/sqrt(3)), eta, 1/sqrt(3)));
    ydifxi = double(subs(subs(sum(diff(Ni, xsi) * coords(:, 2)), xsi, 1/sqrt(3)), eta, 1/sqrt(3)));
    ydifeta = double(subs(subs(sum(diff(Ni, eta) * coords(:, 2)), xsi, 1/sqrt(3)), eta, 1/sqrt(3)));
    
    dNi = [diff(N1, xsi), diff(N2, xsi), diff(N3, xsi), diff(N4, xsi), diff(N5, xsi), diff(N6, xsi), diff(N7, xsi), diff(N8, xsi);
           diff(N1, eta), diff(N2, eta), diff(N3, eta), diff(N4, eta), diff(N5, eta), diff(N6, eta), diff(N7, eta), diff(N8, eta)];
    dNi = subs(subs(dNi, xsi, 1/sqrt(3)), eta, 1/sqrt(3));
    
    J = [diff(x, xsi), diff(x, eta); diff(y, xsi), diff(y, eta)];
    
    JTransInv = (inv(J)).';
    
    
    B1f = [1, 0, 0, 0; 0, 0, 0, 1; 0, 1, 1, 0];
    B2f = zeros(4, 4);
    
    B2f(1:2, 1:2) = JTransInv;
    B2f(3:4, 3:4) = JTransInv;
    
    
    B3f = sym('B3c', [4, 3*length(Ni)]);%zeros(3, 3 * length(Ni));
    for i=1:3:length(Ni) * 3
        index = (i - 1) / 3 + 1;
        B3f(1:4, i:i+2) = [0,    diff(Ni(index), xsi),   0;
                           0,    diff(Ni(index), eta),   0;
                           0,    0,                      diff(Ni(index), xsi);
                           0,    0,                      diff(Ni(index), eta)]; 
    end    
    B1c = [1, 0, 1, 0; 0, 1, 0, 1];
    
    B2c = zeros(4, 4);
    B2c(1:2, 1:2) = JTransInv;
    B2c(3:4, 3:4) = [1, 0;0, 1];
    
    B3c = sym('B3f', [4, 3*length(Ni)]);%zeros(3, 3 * length(Ni));
    for i=1:3:length(Ni) * 3
        index = (i - 1) / 3 + 1;
        B3c(1:4, i:i+2) = [diff(Ni(index), xsi),   0,         0;
                           diff(Ni(index), eta),   0,         0;
                           0,                      Ni(index), 0;
                           0,                      0,         Ni(index)]; 
    end
    
    Bf = B1f * B2f * B3f;
    Df = E/(1 - nu^2) * [1, nu, 0; nu, 1, 0;0, 0, (1-nu)/2];
    stiffness = thickness^3 / 12 * int(int(Bf.' * Df * Bf * det(J), xsi, -1, 1), eta, -1, 1);
    
    %stiffness = thickness^3 / 12 * 4 * subs(subs(Bf.' * Df * Bf * det(J), xsi, 0), 0);
    
    Bc = B1c * B2c * B3c;
    Dc = [G, 0; 0, G];
    
    %% INTEGRATION (CHANGE TO GAUSS QUADRATIC INTEGRATION AFTER WORKING)
    stiffness = stiffness + alpha * thickness * int(int(Bc.' * Dc * Bc * det(J), xsi, -1, 1), eta, -1, 1);
    %stiffness = stiffness + alpha * thickness * 4 * subs(subs(Bc.' * Dc * Bc * det(J), xsi, 0), eta, 0);
    %int(int(Bc.' * Dc * Bc * det(J), xsi, -1, 1), eta, -1, 1);
    
    
    stiffness = double(stiffness);
end