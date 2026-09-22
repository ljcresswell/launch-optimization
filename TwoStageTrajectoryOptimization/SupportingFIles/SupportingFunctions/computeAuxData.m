function aux = computeAuxData(T, X, params, stage_all)

    N = length(T);

    aux.T = zeros(N,1);
    aux.D  = zeros(N,1);
    aux.mdot = zeros(N,1);
    aux.gamma = zeros(N,1);
    aux.delta = zeros(N,1);
    aux.V = zeros(N,1);

    for i = 1:N

        t = T(i);
        Xi = X(i,:)';

        r = Xi(1);
        r_dot = Xi(3);
        v_theta = Xi(4);
        m = Xi(5);

        delta = params.deltaInterp(t);

        V = sqrt(r_dot^2 + v_theta^2);
        gamma = atan2(r_dot, v_theta);

        inAtmo = (r < params.rEndAtmo);
        
        stage= stage_all(i);

        switch stage

            case 1
                delta = 0;

                m_dot = m_dotFunc(params.specificHeatRatio1, ...
                                  params.At1, params.totalPressure1, ...
                                  params.totalTemp1, params.gasConstant1);

                Tforce = thrust(params.Ae1, params.At1, ...
                                params.specificHeatRatio1, ...
                                params.totalTemp1, params.totalPressure1, ...
                                params.gasConstant1, ...
                                params.radiusEarth, r);

            case 2
                m_dot = 0;
                Tforce = 0;

            case 3
                m_dot = m_dotFunc(params.specificHeatRatio2, ...
                                  params.At2, params.totalPressure2, ...
                                  params.totalTemp2, params.gasConstant2);

                Tforce = thrust(params.Ae2, params.At2, ...
                                params.specificHeatRatio2, ...
                                params.totalTemp2, params.totalPressure2, ...
                                params.gasConstant2, ...
                                params.radiusEarth, r);

            case 4
                m_dot = 0;
                Tforce = 0;
        end

        if inAtmo
            D = dragFunc(r, params.radiusEarth, V, ...
                         params.dragCoefficient, params.referenceArea);
        else
            D = 0;
        end

        % store
        aux.T(i)        = Tforce;
        aux.D(i)        = D;
        aux.mdot(i)     = m_dot;
        aux.gamma(i)    = gamma;
        aux.delta(i)    = delta;
        aux.V(i)        = V;
        aux.stage(i)    = stage;

    end
end