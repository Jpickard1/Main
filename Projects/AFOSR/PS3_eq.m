function dydt = PS5_eq(t,y,c)
%     Quadratic         Turing Connection
dy1 = -((y(1).^2+y(2).^2).*y(1))+(y(1)-y(2))+c*(y(3)-y(1));
dy2 = -((y(1).^2+y(2).^2).*y(2))+(y(2)+y(1))+c*(y(4)-y(2));

dy3 = -((y(3).^2+y(4).^2).*y(3))+(y(3)-y(4))+c*(y(1)-y(3));
dy4 = -((y(3).^2+y(4).^2).*y(4))+(y(4)+y(3))+c*(y(2)-y(4));

dydt = [dy1;dy2;dy3;dy4];
return;