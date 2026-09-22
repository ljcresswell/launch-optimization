function m_dot = m_dotFunc(gamma, At, P0, T0, R)
    pwr = -(gamma+1)/(2*gamma-2);
    m_dot = At*P0/sqrt(T0)*sqrt(gamma/R)*((gamma+1)/2)^pwr;
end