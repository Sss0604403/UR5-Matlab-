function Q = ObjectiveFunction_FOPID(params)
    Kp = params(1);
    Ki = params(2);
    Kd = params(3);
    lambda = params(4);
    mu = params(5);
    w1 = 0.999;
    w2 = 0.001;
    w3 = 100;

    model_name = 'ISSA_FOPID_Controller';
    sim_time = 4;
    try
        simOut = sim(model_name, 'SimulationMode', 'normal', 'StopTime', num2str(sim_time), 'SrcWorkspace', 'current');
        e_timeseries = simOut.e;
        u_timeseries = simOut.u;

        e_data = e_timeseries.Data;
        t_data = e_timeseries.Time;
        u_data = u_timeseries.Data;

        integrand1 = w1 * abs(e_data) + w2 * u_data.^2;
        integral_part = trapz(t_data, integrand1);

        penalty_part = 0;
        if any(e_data < 0)
            overshooot_indices = find(e_data < 0);
            penalty_integrand = w3 * abs(e_data(overshooot_indices));
            penalty_part = trapz(t_data(overshooot_indices), penalty_integrand);
        end
        Q = integral_part + penalty_part;
    catch ME
        fprintf('仿真失败。 参数：%s\n', num2str(params));
        disp(ME.message);
        Q = 1e10;
    end
end


