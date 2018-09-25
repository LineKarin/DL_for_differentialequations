function un = BDF2(U,deltaT,g,c)
un = [1+g -deltaT; c^2*deltaT 1+g]\((1+2*g)*U(:,2)-g*U(:,1));
end