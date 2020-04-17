
# Model description

## Fully-connected SIR

The population is composed of $N$ indistinguishable individuals,
each interacting equally with all the others.  An individual can be in
one of three states S (susceptible), I (infected) or R (recovered). 
The dynamics is defined by the transition rates (per individual)

\begin{align*}
  W^{(i)}_{S\to I} &= \beta \frac{N_I}{N-1}, \\
  W^{(i)}_{I\to R} &= \gamma,
\end{align*}

where $N_I$ is the number of individuals in state I.  The rate
constants $\beta$ and $\gamma$ are related to the infection time
$t_\text{inf}$ and $R_0$ through

$$ \gamma = \frac{1}{t_\text{inf}}, \qquad \beta = R_0 \gamma.$

The [parameter file](./sir_par.dat) gives $R_0$ and $t_\text{inf}$ as
well as the initial conditions $S_0$ and $I_0$, namely the initial
fraction of susceptible and infected.


## SIR with families

This is a SIR as above, but with the population divided in families.
The infection rate now has two terms, one which is proportional to the
global concentration of infected, and one proportional to the number
of infected in the same family.  This is a way to mimick social
confinement measures.  The rates for individual $i$ belonging to
family $f$ are

\begin{align*}
  W^{(i)}_{S\to I} &= \beta_\text{out} \frac{N_I}{N-1} + \beta_\text{in} N_{I,f}  , \\
  W^{(i)}_{I\to R} &= \gamma,
\end{align*}

where $N$ is the total population, $N_I$ the total number of infected,
and $N_{I,f}$ the number of infected in family $f$.  Families have $M$
members, with $M=1,\ldots,N_\text{max}$ and distributed with
probabilities $P_M$.

The [parameter file](./sir_f_par.dat) specifies the number of
families, $M_\text{max}$, the probabilites $P_M$ (actually weights, as
the need not be normalized), the rate constants $\beta_\text{in}$,
$\beta_\text{out}$ and $\gamma$, and the initial infected fraction $I_0$.
