
# Model description

## Fully-connected SIR

The population is composed of <img src="/model_desc/tex/f9c4988898e7f532b9f826a75014ed3c.svg?invert_in_darkmode&sanitize=true" align=middle width=14.99998994999999pt height=22.465723500000017pt/> indistinguishable individuals,
each interacting equally with all the others.  An individual can be in
one of three states S (susceptible), I (infected) or R (recovered). 
The dynamics is defined by the transition rates (per individual)

<p align="center"><img src="/model_desc/tex/a97856c18c503c5ef63cc76d8c3e5775.svg?invert_in_darkmode&sanitize=true" align=middle width=130.0391466pt height=63.5745132pt/></p>

where <img src="/model_desc/tex/acba44151255e475354dab0ce025c585.svg?invert_in_darkmode&sanitize=true" align=middle width=19.92813569999999pt height=22.465723500000017pt/> is the number of individuals in state I.  The rate
constants <img src="/model_desc/tex/8217ed3c32a785f0b5aad4055f432ad8.svg?invert_in_darkmode&sanitize=true" align=middle width=10.16555099999999pt height=22.831056599999986pt/> and <img src="/model_desc/tex/11c596de17c342edeed29f489aa4b274.svg?invert_in_darkmode&sanitize=true" align=middle width=9.423880949999988pt height=14.15524440000002pt/> are related to the infection time
<img src="/model_desc/tex/7eb44e3c5f2ccb05ae52bafea039a0b4.svg?invert_in_darkmode&sanitize=true" align=middle width=21.91224089999999pt height=20.221802699999984pt/> and <img src="/model_desc/tex/12d208b4b5de7762e00b1b8fb5c66641.svg?invert_in_darkmode&sanitize=true" align=middle width=19.034022149999988pt height=22.465723500000017pt/> through

$<img src="/model_desc/tex/1773c59dba0725c6f84607dfa23930a5.svg?invert_in_darkmode&sanitize=true" align=middle width=160.2571509pt height=27.77565449999998pt/>

The [parameter file](./sir_par.dat) gives <img src="/model_desc/tex/12d208b4b5de7762e00b1b8fb5c66641.svg?invert_in_darkmode&sanitize=true" align=middle width=19.034022149999988pt height=22.465723500000017pt/> and <img src="/model_desc/tex/7eb44e3c5f2ccb05ae52bafea039a0b4.svg?invert_in_darkmode&sanitize=true" align=middle width=21.91224089999999pt height=20.221802699999984pt/> as
well as the initial conditions <img src="/model_desc/tex/7ea1d381135a7e0a3cd0eefafea7c973.svg?invert_in_darkmode&sanitize=true" align=middle width=16.632471899999988pt height=22.465723500000017pt/> and <img src="/model_desc/tex/88fbd05154e7d6a65883f20e1b18a817.svg?invert_in_darkmode&sanitize=true" align=middle width=13.77859724999999pt height=22.465723500000017pt/>, namely the initial
fraction of susceptible and infected.


## SIR with families

This is a SIR as above, but with the population divided in families.
The infection rate now has two terms, one which is proportional to the
global concentration of infected, and one proportional to the number
of infected in the same family.  This is a way to mimick social
confinement measures.  The rates for individual <img src="/model_desc/tex/77a3b857d53fb44e33b53e4c8b68351a.svg?invert_in_darkmode&sanitize=true" align=middle width=5.663225699999989pt height=21.68300969999999pt/> belonging to
family <img src="/model_desc/tex/190083ef7a1625fbc75f243cffb9c96d.svg?invert_in_darkmode&sanitize=true" align=middle width=9.81741584999999pt height=22.831056599999986pt/> are

<p align="center"><img src="/model_desc/tex/64ab1a073c0f7d5daf632abf7b9293a0.svg?invert_in_darkmode&sanitize=true" align=middle width=222.49021904999998pt height=63.5745132pt/></p>

where <img src="/model_desc/tex/f9c4988898e7f532b9f826a75014ed3c.svg?invert_in_darkmode&sanitize=true" align=middle width=14.99998994999999pt height=22.465723500000017pt/> is the total population, <img src="/model_desc/tex/acba44151255e475354dab0ce025c585.svg?invert_in_darkmode&sanitize=true" align=middle width=19.92813569999999pt height=22.465723500000017pt/> the total number of infected,
and <img src="/model_desc/tex/12d049577400c1be7fd05557a834032f.svg?invert_in_darkmode&sanitize=true" align=middle width=31.53213689999999pt height=22.465723500000017pt/> the number of infected in family <img src="/model_desc/tex/190083ef7a1625fbc75f243cffb9c96d.svg?invert_in_darkmode&sanitize=true" align=middle width=9.81741584999999pt height=22.831056599999986pt/>.  Families have <img src="/model_desc/tex/fb97d38bcc19230b0acd442e17db879c.svg?invert_in_darkmode&sanitize=true" align=middle width=17.73973739999999pt height=22.465723500000017pt/>
members, with <img src="/model_desc/tex/a96e0f223526e1719db57f23c8f3ce99.svg?invert_in_darkmode&sanitize=true" align=middle width=121.87186439999998pt height=22.465723500000017pt/> and distributed with
probabilities <img src="/model_desc/tex/c1c794a6d41ec36f4887107fb625fa05.svg?invert_in_darkmode&sanitize=true" align=middle width=24.323105399999992pt height=22.465723500000017pt/>.

The [parameter file](./sir_f_par.dat) specifies the number of
families, <img src="/model_desc/tex/0f94141a7db1d4304b554728f7b76be6.svg?invert_in_darkmode&sanitize=true" align=middle width=40.205640749999986pt height=22.465723500000017pt/>, the probabilites <img src="/model_desc/tex/c1c794a6d41ec36f4887107fb625fa05.svg?invert_in_darkmode&sanitize=true" align=middle width=24.323105399999992pt height=22.465723500000017pt/> (actually weights, as
the need not be normalized), the rate constants <img src="/model_desc/tex/f5c11fac6ffc2dccb58b3f98e6047614.svg?invert_in_darkmode&sanitize=true" align=middle width=20.279762249999987pt height=22.831056599999986pt/>,
<img src="/model_desc/tex/a549abb698d8c2c1db04d5502c7e297f.svg?invert_in_darkmode&sanitize=true" align=middle width=28.24783829999999pt height=22.831056599999986pt/> and <img src="/model_desc/tex/11c596de17c342edeed29f489aa4b274.svg?invert_in_darkmode&sanitize=true" align=middle width=9.423880949999988pt height=14.15524440000002pt/>, and the initial infected fraction <img src="/model_desc/tex/88fbd05154e7d6a65883f20e1b18a817.svg?invert_in_darkmode&sanitize=true" align=middle width=13.77859724999999pt height=22.465723500000017pt/>.
