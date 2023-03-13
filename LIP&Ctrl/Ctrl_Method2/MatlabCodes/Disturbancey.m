function d=Disturbancey(R)

t=R(1); step=R(2);

d=0;

F=100; m=60;
% F=335; m=60;
if step==3 && t<=0.1
    d=-F/m;
end
end