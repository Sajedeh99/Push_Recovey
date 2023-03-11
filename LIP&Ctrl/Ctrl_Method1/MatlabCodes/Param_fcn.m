function [Tnom,Lnom,Wnom,tau_nom]=Param_fcn(L_min,L_max,W_min,W_max,T_min,T_max,w0,Vx,Vy)

Bl1=L_min/abs(Vx);
Bu1=L_max/abs(Vx);
Bl2=W_min/abs(Vy);
Bu2=W_max/abs(Vy);
Bl3=T_min;
Bu3=T_max;
if Vx==0
    Bl=max(Bl2,Bl3);
    Bu=min(Bu2,Bu3);
elseif Vy==0
    Bl=max(Bl1,Bl3);
    Bu=min(Bu1,Bu3);
else
    Bl=max([Bl1,Bl2,Bl3]);
    Bu=min([Bu1,Bu2,Bu3]);
end

Tnom=(Bl+Bu)/2;
Lnom=Vx*(Bl+Bu)/2;
Wnom=Vy*(Bl+Bu)/2;
tau_nom=exp(w0*Tnom);

end