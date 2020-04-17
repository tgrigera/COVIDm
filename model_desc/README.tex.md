
# Model description

## Fully-connected SIR

The population is composed of $N$ indistinguishable individuals,
each interacting equally with all the others.  An individual can be in
one of three states S (susceptible), I (infected) or R (recovered). 
The dynamics is defined by the transition rates (per individual)

\begin{align*}
  W^{(i)}_{S\to I} &= \beta \frac{N_I}{N-1}, \\
  W^{(i)}_{I\to R} &= \gamma.
\end{align*}

The rate constants $\beta$ and $\gamma$ are related to the infection
time $t_\text{inf}$ and $R_0$ through

$$ \gamma = \frac{1}{t_\text{inf}}, \qquad \beta = R_0 \gamma.$

The [./sir_par.dat][parameter file] gives $R_0$ and $t_\text{inf}$ as
well as the initial conditions $S_0$ and $I_0$, namely the initial
fraction of susceptible and infected.

