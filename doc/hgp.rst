****
Head-Gordon and Pople Method
****

=======
Primitive Function
=======

.....
Function Defined
.....

.. math::
    :nowrap:

    \begin{equation}
        [a]=\chi(l_a,m_a,n_a,\alpha,\mathbf{A})=x_A^{l_a}y_A^{m_a}z_A^{n_a} exp(-\alpha_a r_A^2)
    \end{equation}
    
    \begin{equation}
        [a+1_x]=\chi(l_a+1,m_a,n_a,\alpha,\mathbf{A})=x_A^{l_a+1}y_A^{m_a}z_A^{n_a} exp(-\alpha r_A^2)
    \end{equation}
    ......
    
    \begin{equation}
        0=\chi(0,0,0,\alpha,\mathbf{A})=exp(-\alpha r_A^2)
    \end{equation}


.....
Recurrence Relation
.....

.. math::

    \begin{equation}
    \begin{split}
        [(a+1_i)b, cd]^m &= \mathbf{PA}_i[ab,cd]^m + \mathbf{WP}_i[ab,cd]^{m+1} \\
                         &+ \frac{N_i(a)}{2\zeta}\left\{ [(a-1_i)b,cd]^m - \frac{\rho}{\zeta}[(a-1_i)b,cd]^{m+1}\right\} \\
                         &+ \frac{N_i(b)}{2\zeta}\left\{ [(a-1_i)b,cd]^m - \frac{\rho}{\eta}[(a-1_i)b,cd]^{m+1}\right\} \\
                         &+ \frac{N_i(c)}{\zeta+\eta}[ab,(c-1_i)d]^{m+1} \\
                         &+ \frac{N_i(d)}{\zeta+\eta}[ab,c(d-1_i)]^{m+1}
    \end{split}
    \end{equation}
    
    \begin{equation}
    \begin{split}
        [ab, (c+1_i)d]^m &= \mathbf{QC}_i[ab,cd]^m + \mathbf{WQ}_i[ab,cd]^{m+1} \\
                         &+ \frac{N_i(c)}{2\zeta}\left\{ [ab,(c-1_i)d]^m - \frac{\rho}{\zeta}[ab,(c-1_i)d]^{m+1}\right\} \\
                         &+ \frac{N_i(d)}{2\zeta}\left\{ [ab,c(d-1_i)]^m - \frac{\rho}{\eta}[ab,c(d-1_i)]^{m+1}\right\} \\
                         &+ \frac{N_i(a)}{\zeta+\eta}[(a-1_i)b,cd]^{m+1} \\
                         &+ \frac{N_i(b)}{\zeta+\eta}[a(b-1_i),cd]^{m+1}
    \end{split}
    \end{equation}

    where i=x, y, z and \( N_i(a)=l_a, m_a, n_a \)
    \[ \zeta=\alpha_a+\alpha_b \]
    \[ \eta =\alpha_c+\alpha_d \]
    \[ \rho = \frac{\zeta\eta}{\zeta+\eta} \]
    \[ \mathbf{P} = \frac{\alpha_a \mathbf{A} + \alpha_b\mathbf{B}}{\alpha_a + \alpha_b} \]
    \[ \mathbf{Q} = \frac{\alpha_c \mathbf{C} + \alpha_d\mathbf{D}}{\alpha_c + \alpha_d} \]
    \[ \mathbf{W} = \frac{\zeta \mathbf{P} + \eta\mathbf{Q}}{\zeta + \eta} \]


    Finally, we should obtain:
    \[ [ab,cd] = \sum_mC_m[00,00]^m \]
    and we can get the value of \( [00,00]^m \)

========
Basis Function
========

.. math::
    :nowrap:

    \begin{equation}
        (a)=\phi=\sum_u c_u\chi_u
    \end{equation}
    Where:  \(c_u\) is a constant coefficient.

    \begin{equation}
        (ab,cd)=(\chi_A\chi_B\mid\chi_C\chi_D)=\sum_u\sum_v\sum_s\sum_t[ab,cd]
    \end{equation}

.....
Recurrence Relation
.....

.. math::
    :nowrap:

    \begin{equation}
        (a(b+1_i),cd)=((a+1_i)b,cd)+\mathbf{AB}_i(ab,cd)
    \end{equation}
    
    \begin{equation}
        (ab,c(d+1_i))=(ab,(c+1_i)d)+\mathbf{CD}_i(ab,cd)
    \end{equation}
    
    
    Finally, we should obtain:
    \[(a(b+1_i),cd) = (e0,f0)+...\]

