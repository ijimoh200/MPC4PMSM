function x_next = runge_kutta_pmsm(states_ini,Ts,u_c,di)
K1 = pmsm_plant([],states_ini,u_c,di);
K2 = pmsm_plant([],states_ini+K1.*(Ts/2),u_c,di);
K3 = pmsm_plant([],states_ini+K2.*(Ts/2),u_c,di);
K4 = pmsm_plant([],states_ini+K3.*Ts,u_c,di);

x_next = states_ini + (Ts/6).*(K1+2.*K2+2.*K3 + K4);
end