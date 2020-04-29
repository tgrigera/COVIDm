
# Model description

The models implemented in COVIDm are the following.

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

$$ \gamma = \frac{1}{t_\text{inf}}, \qquad \beta = R_0 \gamma. $$

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


## SEEIIR

This model has six states instead of the three of SIR: susceptible,
two exposed states, two infected states and recovered.  Exposed means
that the individual will turn into an infected, but does not yet
transmit the disease.  The two E and I states are identical, the
duplicated state serves to give a non-exponential distribution of the
duration of the infection.  The population is divided in families as
in the previous model.  Transition rates are

\begin{align*}
   W^{(i)}_{S\to E_1} & = \beta_\text{out} \frac{N_{I_1}+N_{I_2} } {N-1}  +
                          \beta_\text{in}  ( N_{I_1,f} + N_{I_2,f}), \\
   W^{(i)}_{E_1\to E_2} &= 2  \sigma, \\
   W^{(i)}_{E_2 \to I_1} &= 2  \sigma, \\
   W^{(i)}_{I_1 \to I_2} &= 2  \gamma, \\
   W^{(i)}_{I_2 \to R} &= 2  \gamma, \\
\end{align*}

with $N$ the total population, $N_X$ the number of individuals in
state $X$, and $N_{X,f}$ the number of members of family $f$ in state
$X$.  The [parameter file](./seeir_par.dat) specifies the number and
size distribuiton of families as in the previous case, plus the four
rate constants $\beta_\text{in}$, $\beta_{out}$, $\sigma$ and
$\gamma$.  Finally it gives the name of the file containing imported
infections.

The population is initialized to all S, so infected individuals must
be introduced in order to start the epidemic.  The external infections
are read from a two-column [file](./imported_infections.dat) giving
time and number of imported cases.  The second column is interpreted
as the cumulative total cases, not new cases.  At the specifed time a
number of S individuals are forced to state I$_1$ so that the total
count of imported cases matches the number given in the file.  Both
times and cases must be monotonically increasing.

The program outputs the total number of individuals in each
epidemiological state plus a cumulative count of infections according
to type: imported (this matches the number from the imported
infections file), close contact and communtiy.  Infections are counted
as close-contact if when the invidual transitions from $E_2$ to $I_1$,
the family has other members in states $I_1$, $I_2$ or $R$.  Otherwise
they are classified as community.

One can also give a time-dependent $\beta_\text{out}$ as a way to
model variaitons over time of social confinement measures.  To do
this, set a negative value $\beta_\text{out}$ in the parameter file
and give the name of a file with the
[time-dependent parameter](./beta_vs_time.dat): this is a two-column
file giving a time at which $\beta_\text{out}$ is changed, and its new value.
Be sure to give a reasonable value at $t=0$, otherwise the simulation will start
with a negative $\beta$.


## Hierarchical SEEIIR

The hierarchical SEEIIR model is a generalization of the SEEIIR with
families as described above, with several hierarchical levels of
aggregation: individuals are grouped in families, which are grouped in
neighborhoods, which are grouped in towns, etc.  Infections can result
through contacts with individuals of the same family, or neighborhood,
etc, and each of these levels carries a different $\beta$.

The individual that contracts the disease goes from state $S$ to $E_1$
to $E_2$ to $I_1$ to $I_2$ to $R$ as above, only that we now allow for
for the rate for $E_1\to E_2$ to be different from the rate for
$E_2\to I_1$ (and similarly for $I_1\to I_2$ and $I_2\to R$).  The $E$
states are not contagious.

The different aggregation categories (families, neighborhoods, etc.)
are called _levels_.  Level 0 is the individual, level 1 the family
and so on.  The [parameter file](./seeiir_h_par.dat) specifies the
number $L$ of levels, and gives at each level the number $M_l$ of
entities of level $l-1$ for each level $l$ group.  For example, if
$L=3$ there will be only one group at $l=3$ (the whole population),
$M_3$ groups at level 2, $M_2$ groups at level 1 (families) and each
level 1 group will have $M_1$ individuals.  It is possible to ask that
$M_l$ be randomly variable: if $M_l$ is negative, the program takes
$|M_l|$ as the maximum allowed $M_l$, then reads the probabilies for
$M_l$ taking the values $1,\ldots,|M_l|$ (actually the numbers are
weights rather than probabilities, since they do not need to be
normalized).

To write the transition rates, give all groups a unique number, and
let $\nu_{l,i}$ be the number of the group of level $l$ to which
individual $i$ belongs.  Calling $N_{\nu_{l,i}}(X)$ the number of
individuals of group $\nu_{i,l}$ in state $X$, the transition rates
are

\begin{align*}
   W^{(i)}_{S\to E_1} & = \sum_{l=2}^L \beta_l \frac{N_{\nu_{l,i}}(I_1)+N_{\nu_{i,l}}(I_2) }
                     {N_{\nu_{i,l}}-1}  +
                          \beta_1  \left( N_{\nu_{i,1}}(I_1) + N_{\nu_{i,1}}(I_2) \right), \\
   W^{(i)}_{E_1\to E_2} &=   \sigma_1, \\
   W^{(i)}_{E_2 \to I_1} &=  \sigma_2, \\
   W^{(i)}_{I_1 \to I_2} &=   \gamma_1, \\
   W^{(i)}_{I_2 \to R} &=   \gamma_2, \\
\end{align*}

The L+4 rate constants ($\beta_1,\ldots,\beta_L$, $\sigma_1$,
$\sigma_2$, $\gamma_1$, $\gamma_2$) are given at the end of the
parameter file, as a function of time: at the time given in the first
column, the constants are set to the given values.  It is necessary to
give a line with the values at time 0, additional lines with different
times are only given if changes are desired.

Finally, as in the SEEIIR model above, the population is initialized
to all S, so infected individuals must be introduced in order to start
the epidemic.  The external infections are read from a two-column
[file](./imported_infections.dat) giving time and number of imported
cases.  The second column is interpreted as the cumulative total
cases, not new cases.  At the specifed time a number of S individuals
are forced to state $I_1$ so that the total count of imported cases
matches the number given in the file.  Both times and cases must be
monotonically increasing.
