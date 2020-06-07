function Ham=funz(m,vinf,w_frac_vinf)


m=m_add_unsteady_loads(m,[1 0 0]',w_frac_vinf);

Ham = m.Ham;

end