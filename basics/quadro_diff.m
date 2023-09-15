function [dy] = quadro_diff(x,y);
% Дифференцирование функции, кусочно аппроксимируемой квадратичным
% полиномом Лагранжа (по три точки)
for i = 2:(length(x)-1)
    
    x1 = [ x(i-1) x(i) x(i+1) ];
    y1 = [ y(i-1) y(i) y(i+1) ];
    [aa] = lagrange_polynom_2(x1,y1,0);
    
    dy(i) = 2*aa(1)*x(i) + aa(2);
    
    
    if i == 2
%         dy(1) = 2*aa(1)*x(1) + aa(2);
        [a,b] = Line_coef(x(1),y(1),x(2),y(2));
        dy(1) = a;
    end
    if i == (length(x) - 1)
%         dyend = 2*aa(1)*x(end) + aa(2);
        [a,b] = Line_coef(x(end-1),y(end-1),x(end),y(end));
        dyend = a;
    end
end
dy = [ dy dyend ];
% dy = dy(2:end);
% xx = x(2:(end-1));
% dyx = interp1(xx,dy,x(1));
% dy = [ dyx dy ];
% dy(end+1) = interp1(x(1:(end-1)),dy,x(end));
