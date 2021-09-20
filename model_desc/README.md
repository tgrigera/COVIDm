
# Model description

The models implemented in COVIDm are the following.

## Fully-connected SIR

The population is composed of <img src="svgs//f9c4988898e7f532b9f826a75014ed3c.svg?invert_in_darkmode" align=middle width=14.99998994999999pt height=22.465723500000017pt/> indistinguishable individuals,
each interacting equally with all the others.  An individual can be in
one of three states S (susceptible), I (infected) or R (recovered). 
The dynamics is defined by the transition rates (per individual)

<p align="center"><img src="svgs//a97856c18c503c5ef63cc76d8c3e5775.svg?invert_in_darkmode" align=middle width=130.0391466pt height=63.5745132pt/></p>

where <img src="svgs//acba44151255e475354dab0ce025c585.svg?invert_in_darkmode" align=middle width=19.92813569999999pt height=22.465723500000017pt/> is the number of individuals in state I.  The rate
constants <img src="svgs//8217ed3c32a785f0b5aad4055f432ad8.svg?invert_in_darkmode" align=middle width=10.16555099999999pt height=22.831056599999986pt/> and <img src="svgs//11c596de17c342edeed29f489aa4b274.svg?invert_in_darkmode" align=middle width=9.423880949999988pt height=14.15524440000002pt/> are related to the infection time
<img src="svgs//7eb44e3c5f2ccb05ae52bafea039a0b4.svg?invert_in_darkmode" align=middle width=21.91224089999999pt height=20.221802699999984pt/> and <img src="svgs//12d208b4b5de7762e00b1b8fb5c66641.svg?invert_in_darkmode" align=middle width=19.034022149999988pt height=22.465723500000017pt/> through

<p align="center"><img src="svgs//2497d8a0a73a80bfd4b178c2c4e07841.svg?invert_in_darkmode" align=middle width=163.213677pt height=35.45589465pt/></p>

The [parameter file](./sir_par.dat) gives <img src="svgs//12d208b4b5de7762e00b1b8fb5c66641.svg?invert_in_darkmode" align=middle width=19.034022149999988pt height=22.465723500000017pt/> and <img src="svgs//7eb44e3c5f2ccb05ae52bafea039a0b4.svg?invert_in_darkmode" align=middle width=21.91224089999999pt height=20.221802699999984pt/> as
well as the initial conditions <img src="svgs//7ea1d381135a7e0a3cd0eefafea7c973.svg?invert_in_darkmode" align=middle width=16.632471899999988pt height=22.465723500000017pt/> and <img src="svgs//88fbd05154e7d6a65883f20e1b18a817.svg?invert_in_darkmode" align=middle width=13.77859724999999pt height=22.465723500000017pt/>, namely the initial
fraction of susceptible and infected.


## SIR with families

This is a SIR as above, but with the population divided in families.
The infection rate now has two terms, one which is proportional to the
global concentration of infected, and one proportional to the number
of infected in the same family.  This is a way to mimick social
confinement measures.  The rates for individual <img src="svgs//77a3b857d53fb44e33b53e4c8b68351a.svg?invert_in_darkmode" align=middle width=5.663225699999989pt height=21.68300969999999pt/> belonging to
family <img src="svgs//190083ef7a1625fbc75f243cffb9c96d.svg?invert_in_darkmode" align=middle width=9.81741584999999pt height=22.831056599999986pt/> are

<p align="center"><img src="svgs//64ab1a073c0f7d5daf632abf7b9293a0.svg?invert_in_darkmode" align=middle width=222.49021904999998pt height=63.5745132pt/></p>

where <img src="svgs//f9c4988898e7f532b9f826a75014ed3c.svg?invert_in_darkmode" align=middle width=14.99998994999999pt height=22.465723500000017pt/> is the total population, <img src="svgs//acba44151255e475354dab0ce025c585.svg?invert_in_darkmode" align=middle width=19.92813569999999pt height=22.465723500000017pt/> the total number of infected,
and <img src="svgs//12d049577400c1be7fd05557a834032f.svg?invert_in_darkmode" align=middle width=31.53213689999999pt height=22.465723500000017pt/> the number of infected in family <img src="svgs//190083ef7a1625fbc75f243cffb9c96d.svg?invert_in_darkmode" align=middle width=9.81741584999999pt height=22.831056599999986pt/>.  Families have <img src="svgs//fb97d38bcc19230b0acd442e17db879c.svg?invert_in_darkmode" align=middle width=17.73973739999999pt height=22.465723500000017pt/>
members, with <img src="svgs//a96e0f223526e1719db57f23c8f3ce99.svg?invert_in_darkmode" align=middle width=121.87186439999998pt height=22.465723500000017pt/> and distributed with
probabilities <img src="svgs//c1c794a6d41ec36f4887107fb625fa05.svg?invert_in_darkmode" align=middle width=24.323105399999992pt height=22.465723500000017pt/>.

The [parameter file](./sir_f_par.dat) specifies the number of
families, <img src="svgs//0f94141a7db1d4304b554728f7b76be6.svg?invert_in_darkmode" align=middle width=40.205640749999986pt height=22.465723500000017pt/>, the probabilites <img src="svgs//c1c794a6d41ec36f4887107fb625fa05.svg?invert_in_darkmode" align=middle width=24.323105399999992pt height=22.465723500000017pt/> (actually weights, as
the need not be normalized), the rate constants <img src="svgs//f5c11fac6ffc2dccb58b3f98e6047614.svg?invert_in_darkmode" align=middle width=20.279762249999987pt height=22.831056599999986pt/>,
<img src="svgs//a549abb698d8c2c1db04d5502c7e297f.svg?invert_in_darkmode" align=middle width=28.24783829999999pt height=22.831056599999986pt/> and <img src="svgs//11c596de17c342edeed29f489aa4b274.svg?invert_in_darkmode" align=middle width=9.423880949999988pt height=14.15524440000002pt/>, and the initial infected fraction <img src="svgs//88fbd05154e7d6a65883f20e1b18a817.svg?invert_in_darkmode" align=middle width=13.77859724999999pt height=22.465723500000017pt/>.


## SEEIIR

This model has six states instead of the three of SIR: susceptible,
two exposed states, two infected states and recovered.  Exposed means
that the individual will turn into an infected, but does not yet
transmit the disease.  The two E and I states are identical, the
duplicated state serves to give a non-exponential distribution of the
duration of the infection.  The population is divided in families as
in the previous model.  Transition rates are

<p align="center"><img src="svgs//3f55e2889ffa9c70c1236a80b5cb8c2f.svg?invert_in_darkmode" align=middle width=343.537062pt height=155.87549669999999pt/></p>

with <img src="svgs//f9c4988898e7f532b9f826a75014ed3c.svg?invert_in_darkmode" align=middle width=14.99998994999999pt height=22.465723500000017pt/> the total population, <img src="svgs//4855e1a0b6043df68024a06e6eb9d355.svg?invert_in_darkmode" align=middle width=24.882493349999987pt height=22.465723500000017pt/> the number of individuals in
state <img src="svgs//cbfb1b2a33b28eab8a3e59464768e810.svg?invert_in_darkmode" align=middle width=14.908688849999992pt height=22.465723500000017pt/>, and <img src="svgs//f1e6937ee9361d7fcc96be9f2eb2ed5d.svg?invert_in_darkmode" align=middle width=35.778728699999995pt height=22.465723500000017pt/> the number of members of family <img src="svgs//190083ef7a1625fbc75f243cffb9c96d.svg?invert_in_darkmode" align=middle width=9.81741584999999pt height=22.831056599999986pt/> in state
<img src="svgs//cbfb1b2a33b28eab8a3e59464768e810.svg?invert_in_darkmode" align=middle width=14.908688849999992pt height=22.465723500000017pt/>.  The [parameter file](./seeir_par.dat) specifies the number and
size distribuiton of families as in the previous case, plus the four
rate constants <img src="svgs//f5c11fac6ffc2dccb58b3f98e6047614.svg?invert_in_darkmode" align=middle width=20.279762249999987pt height=22.831056599999986pt/>, <img src="svgs//b0adc8c70f0031d9ca1ea378d3a13c63.svg?invert_in_darkmode" align=middle width=28.52451029999999pt height=22.831056599999986pt/>, <img src="svgs//8cda31ed38c6d59d14ebefa440099572.svg?invert_in_darkmode" align=middle width=9.98290094999999pt height=14.15524440000002pt/> and
<img src="svgs//11c596de17c342edeed29f489aa4b274.svg?invert_in_darkmode" align=middle width=9.423880949999988pt height=14.15524440000002pt/>.  Finally it gives the name of the file containing imported
infections.

The population is initialized to all S, so infected individuals must
be introduced in order to start the epidemic.  The external infections
are read from a two-column [file](./imported_infections.dat) giving
time and number of imported cases.  The second column is interpreted
as the cumulative total cases, not new cases.  At the specifed time a
number of S individuals are forced to state I<img src="svgs//d0d302843978b7a7264dbb9cb2196dfd.svg?invert_in_darkmode" align=middle width=6.5525476499999895pt height=14.15524440000002pt/> so that the total
count of imported cases matches the number given in the file.  Both
times and cases must be monotonically increasing.

The program outputs the total number of individuals in each
epidemiological state plus a cumulative count of infections according
to type: imported (this matches the number from the imported
infections file), close contact and communtiy.  Infections are counted
as close-contact if when the invidual transitions from <img src="svgs//be655e6ff5e921809983a59b05ec05b4.svg?invert_in_darkmode" align=middle width=18.687266399999988pt height=22.465723500000017pt/> to <img src="svgs//d906cd9791e4b48a3b848558acda5899.svg?invert_in_darkmode" align=middle width=13.77859724999999pt height=22.465723500000017pt/>,
the family has other members in states <img src="svgs//d906cd9791e4b48a3b848558acda5899.svg?invert_in_darkmode" align=middle width=13.77859724999999pt height=22.465723500000017pt/>, <img src="svgs//9eff113852463b85a970d2d65d52280c.svg?invert_in_darkmode" align=middle width=13.77859724999999pt height=22.465723500000017pt/> or <img src="svgs//1e438235ef9ec72fc51ac5025516017c.svg?invert_in_darkmode" align=middle width=12.60847334999999pt height=22.465723500000017pt/>.  Otherwise
they are classified as community.

One can also give a time-dependent <img src="svgs//a549abb698d8c2c1db04d5502c7e297f.svg?invert_in_darkmode" align=middle width=28.24783829999999pt height=22.831056599999986pt/> as a way to
model variaitons over time of social confinement measures.  To do
this, set a negative value <img src="svgs//a549abb698d8c2c1db04d5502c7e297f.svg?invert_in_darkmode" align=middle width=28.24783829999999pt height=22.831056599999986pt/> in the parameter file
and give the name of a file with the
[time-dependent parameter](./beta_vs_time.dat): this is a two-column
file giving a time at which <img src="svgs//a549abb698d8c2c1db04d5502c7e297f.svg?invert_in_darkmode" align=middle width=28.24783829999999pt height=22.831056599999986pt/> is changed, and its new value.
Be sure to give a reasonable value at <img src="svgs//1c899e1c767eb4eac89facb5d1f2cb0d.svg?invert_in_darkmode" align=middle width=36.07293689999999pt height=21.18721440000001pt/>, otherwise the simulation will start
with a negative <img src="svgs//8217ed3c32a785f0b5aad4055f432ad8.svg?invert_in_darkmode" align=middle width=10.16555099999999pt height=22.831056599999986pt/>.


## Hierarchical SEEIIR

The hierarchical SEEIIR model is a generalization of the SEEIIR with
families as described above, with several hierarchical levels of
aggregation: individuals are grouped in families, which are grouped in
neighborhoods, which are grouped in towns, etc.  Infections can result
through contacts with individuals of the same family, or neighborhood,
etc, and each of these levels carries a different <img src="svgs//8217ed3c32a785f0b5aad4055f432ad8.svg?invert_in_darkmode" align=middle width=10.16555099999999pt height=22.831056599999986pt/>.

The individual that contracts the disease goes from state <img src="svgs//e257acd1ccbe7fcb654708f1a866bfe9.svg?invert_in_darkmode" align=middle width=11.027402099999989pt height=22.465723500000017pt/> to <img src="svgs//cc5d1b1ed1bb46a5b9c4cb510b29c8d8.svg?invert_in_darkmode" align=middle width=18.687266399999988pt height=22.465723500000017pt/>
to <img src="svgs//be655e6ff5e921809983a59b05ec05b4.svg?invert_in_darkmode" align=middle width=18.687266399999988pt height=22.465723500000017pt/> to <img src="svgs//d906cd9791e4b48a3b848558acda5899.svg?invert_in_darkmode" align=middle width=13.77859724999999pt height=22.465723500000017pt/> to <img src="svgs//9eff113852463b85a970d2d65d52280c.svg?invert_in_darkmode" align=middle width=13.77859724999999pt height=22.465723500000017pt/> to <img src="svgs//1e438235ef9ec72fc51ac5025516017c.svg?invert_in_darkmode" align=middle width=12.60847334999999pt height=22.465723500000017pt/> as above, only that we now allow for
for the rate for <img src="svgs//1129ab5181dbe8cda1795ec5fccdafec.svg?invert_in_darkmode" align=middle width=63.76704509999999pt height=22.465723500000017pt/> to be different from the rate for
<img src="svgs//4636946a03c94fd533ae50ebd4c2313e.svg?invert_in_darkmode" align=middle width=58.85837594999998pt height=22.465723500000017pt/> (and similarly for <img src="svgs//fa8ccd846b012105ff0c97ed5e7841be.svg?invert_in_darkmode" align=middle width=53.94970844999998pt height=22.465723500000017pt/> and <img src="svgs//1f444b2c14bf25cb4fe0d9355f156807.svg?invert_in_darkmode" align=middle width=52.77958454999999pt height=22.465723500000017pt/>).  The <img src="svgs//84df98c65d88c6adf15d4645ffa25e47.svg?invert_in_darkmode" align=middle width=13.08219659999999pt height=22.465723500000017pt/>
states are not contagious.

The different aggregation categories (families, neighborhoods, etc.)
are called _levels_.  Level 0 is the individual, level 1 the family
and so on.  The [parameter file](./seeiir_h_par.dat) specifies the
number <img src="svgs//ddcb483302ed36a59286424aa5e0be17.svg?invert_in_darkmode" align=middle width=11.18724254999999pt height=22.465723500000017pt/> of levels, and gives at each level the number <img src="svgs//a28763d9579bce99319cd8520393ee02.svg?invert_in_darkmode" align=middle width=20.171302799999992pt height=22.465723500000017pt/> of
entities of level <img src="svgs//b4ba0aa01606a3064bb00827756e5407.svg?invert_in_darkmode" align=middle width=33.53874479999999pt height=22.831056599999986pt/> for each level <img src="svgs//2f2322dff5bde89c37bcae4116fe20a8.svg?invert_in_darkmode" align=middle width=5.2283516999999895pt height=22.831056599999986pt/> group.  For example, if
<img src="svgs//0e312230f7e3bbd18f2262defcbf18b4.svg?invert_in_darkmode" align=middle width=41.32408334999999pt height=22.465723500000017pt/> there will be only one group at <img src="svgs//4f4dd296c2b8cc8358b3bcfa857581f7.svg?invert_in_darkmode" align=middle width=35.36518424999999pt height=22.831056599999986pt/> (the whole population),
<img src="svgs//7edb093218f7ba2b5c51886398fa9caf.svg?invert_in_darkmode" align=middle width=22.500061649999992pt height=22.465723500000017pt/> groups at level 2, <img src="svgs//dced8cd0d35e2af2d3499c10d7ee6289.svg?invert_in_darkmode" align=middle width=22.500061649999992pt height=22.465723500000017pt/> groups at level 1 (families) and each
level 1 group will have <img src="svgs//6f549764f2f97bec950c14de5352994a.svg?invert_in_darkmode" align=middle width=22.500061649999992pt height=22.465723500000017pt/> individuals.  It is possible to ask that
<img src="svgs//a28763d9579bce99319cd8520393ee02.svg?invert_in_darkmode" align=middle width=20.171302799999992pt height=22.465723500000017pt/> be randomly variable: if <img src="svgs//a28763d9579bce99319cd8520393ee02.svg?invert_in_darkmode" align=middle width=20.171302799999992pt height=22.465723500000017pt/> is negative, the program takes
<img src="svgs//346b9c00599fbecd7cf0c2b7ea91d5ef.svg?invert_in_darkmode" align=middle width=30.12564719999999pt height=24.65753399999998pt/> as the maximum allowed <img src="svgs//a28763d9579bce99319cd8520393ee02.svg?invert_in_darkmode" align=middle width=20.171302799999992pt height=22.465723500000017pt/>, then reads the probabilies for
<img src="svgs//a28763d9579bce99319cd8520393ee02.svg?invert_in_darkmode" align=middle width=20.171302799999992pt height=22.465723500000017pt/> taking the values <img src="svgs//0303c5a1ee56220c3888ff7d03329325.svg?invert_in_darkmode" align=middle width=74.87427254999999pt height=24.65753399999998pt/> (actually the numbers are
weights rather than probabilities, since they do not need to be
normalized).

To write the transition rates, give all groups a unique number, and
let <img src="svgs//47c8695a1a41562d3801560166cd7b65.svg?invert_in_darkmode" align=middle width=20.899093049999987pt height=14.15524440000002pt/> be the number of the group of level <img src="svgs//2f2322dff5bde89c37bcae4116fe20a8.svg?invert_in_darkmode" align=middle width=5.2283516999999895pt height=22.831056599999986pt/> to which
individual <img src="svgs//77a3b857d53fb44e33b53e4c8b68351a.svg?invert_in_darkmode" align=middle width=5.663225699999989pt height=21.68300969999999pt/> belongs.  Calling <img src="svgs//e26c6c2e6496281c577e8a53df97e30a.svg?invert_in_darkmode" align=middle width=61.236859199999984pt height=24.65753399999998pt/> the number of
individuals of group <img src="svgs//4539837d61243d7054004fb10a411d9a.svg?invert_in_darkmode" align=middle width=20.89909469999999pt height=14.15524440000002pt/> in state <img src="svgs//cbfb1b2a33b28eab8a3e59464768e810.svg?invert_in_darkmode" align=middle width=14.908688849999992pt height=22.465723500000017pt/>, the transition rates
are

<p align="center"><img src="svgs//c27f1ee56d6e9b1f453fdb232871394a.svg?invert_in_darkmode" align=middle width=274.78996379999995pt height=170.70283229999998pt/></p>

The L+4 rate constants (<img src="svgs//bf9b6c39addb984a496459e55f667729.svg?invert_in_darkmode" align=middle width=71.51812799999999pt height=22.831056599999986pt/>, <img src="svgs//9811b1e861c58c0f72de45f573e7eea4.svg?invert_in_darkmode" align=middle width=15.94565279999999pt height=14.15524440000002pt/>,
<img src="svgs//0318cc8a44e98dfa8db4cd5b6f731ed4.svg?invert_in_darkmode" align=middle width=15.94565279999999pt height=14.15524440000002pt/>, <img src="svgs//f4d2e02931a81ed37d695634192ebed9.svg?invert_in_darkmode" align=middle width=15.06318824999999pt height=14.15524440000002pt/>, <img src="svgs//b20edd7d59cbe61f6558fa4a05eaacea.svg?invert_in_darkmode" align=middle width=15.06318824999999pt height=14.15524440000002pt/>) are given at the end of the
parameter file, as a function of time: at the time given in the first
column, the constants are set to the given values.  It is necessary to
give a line with the values at time 0, additional lines with different
times are only given if changes are desired.

As in the SEEIIR model above, the population is initialized
to all S, so infected individuals must be introduced in order to start
the epidemic.  The external infections are read from a tree-column
[file](./imported_infections.dat) giving time and number of imported
cases.  The second column is interpreted as the cumulative total
cases, not new cases.  At the specifed time a number of S individuals
are forced to state <img src="svgs//d906cd9791e4b48a3b848558acda5899.svg?invert_in_darkmode" align=middle width=13.77859724999999pt height=22.465723500000017pt/> so that the total count of imported cases
matches the number given in the file.  Both times and cases must be
monotonically increasing.

The third column of this file gives the possibility to force a number
of susceptible individuals to state R, and to make them susceptible
again at a later time.  The number of forced R individuals is given as
a third column in the imported infections file.  This number need not
be monotonically increasing: a decrease indicates that some (picked at
random) _of the previouly forced to R_ individuals are to become
susceptible again.  If the program is `seeiir_h_force_recover_family`,
then the way this is done is to randomly choose a family where all
individuals are in state S or R, and then recover _all_ of the
susceptible individuals of the family.  This may lead to recovering
slightly more individuals than requested.  Similarly, when undoing the
recovery, all of the forcibly recovered members of a family are made
susceptible again.

Finally, there are two possible invocations.  The first one is

`seeiir_h parameterfile seed steps Nruns`

The second requests more details of individual's state at different
levels to be written to a given file:

`seeiir_h parameterfile seed steps Nruns detail_level detail_field detail_file`

`detail_field` is  a single letter (`S`, `I`, or `R`) which specifies what state one wishes to monitor in more detail.  `detail_file` is the name of an output file, which will hold, as a
function of time, the average number of  individuals in the requested state  at each
level (except level <img src="svgs//ddcb483302ed36a59286424aa5e0be17.svg?invert_in_darkmode" align=middle width=11.18724254999999pt height=22.465723500000017pt/> which contains the whole population, which data
is written to standard otput) and its variance.  In addition, for all
levels from <img src="svgs//ddcb483302ed36a59286424aa5e0be17.svg?invert_in_darkmode" align=middle width=11.18724254999999pt height=22.465723500000017pt/> down to `detail_level`, the number of individuals in the requested state at each
node will be written to the same file.


## Models on graphs

The SIR and SEEIIR model are also implemented on several graphs.  The
graphs are described by the connectivity mattrix <img src="svgs//3a91172ac0904d6988776a8ee5bd4140.svg?invert_in_darkmode" align=middle width=19.870705799999985pt height=22.465723500000017pt/>.  Each
element of the matrix is a real number, giving the _weight_ of the
bond between nodes <img src="svgs//77a3b857d53fb44e33b53e4c8b68351a.svg?invert_in_darkmode" align=middle width=5.663225699999989pt height=21.68300969999999pt/> and <img src="svgs//36b5afebdba34564d884d347484ac0c7.svg?invert_in_darkmode" align=middle width=7.710416999999989pt height=21.68300969999999pt/> (a weight of 0 means the nodes are not
connected).  Depending on the graph, <img src="svgs//3a91172ac0904d6988776a8ee5bd4140.svg?invert_in_darkmode" align=middle width=19.870705799999985pt height=22.465723500000017pt/> can be symmetric or not
(directed graph).

For the SIR model, one has

<p align="center"><img src="svgs//1f0f8859f3114b32bee8c732f801a952.svg?invert_in_darkmode" align=middle width=171.7176252pt height=69.11190165pt/></p>

where <img src="svgs//06db13d75f55ad4e807b0855961dc632.svg?invert_in_darkmode" align=middle width=82.63387769999999pt height=22.465723500000017pt/> is the state of individual <img src="svgs//77a3b857d53fb44e33b53e4c8b68351a.svg?invert_in_darkmode" align=middle width=5.663225699999989pt height=21.68300969999999pt/>, and
<img src="svgs//938dcbd8799e6313859df6a8d99ab862.svg?invert_in_darkmode" align=middle width=25.252984349999988pt height=22.831056599999986pt/> is Kronecker's delta, i. e. <img src="svgs//5463ae668d409c54ca02d37420c6e9b4.svg?invert_in_darkmode" align=middle width=30.75530369999999pt height=22.831056599999986pt/>
equals 1 if <img src="svgs//dcf5d53528d6bd547321757ca613e617.svg?invert_in_darkmode" align=middle width=45.29952239999999pt height=22.465723500000017pt/>, and 0 otherwise.

Similarly for SEEIIR,
<p align="center"><img src="svgs//4d372b7b9e672c7318e62562580c8ee4.svg?invert_in_darkmode" align=middle width=269.9906253pt height=161.41288515pt/></p>

The code is designed so that a model is programmed on a generic graph.
Graphs are implemented separately form the models, so that any of the
available models can run on any of the implemented graphs.

The following graphs have been implemented so far.

### Square lattice

<img src="svgs//c9b68e17641d236fb3019541ac2f8246.svg?invert_in_darkmode" align=middle width=50.82941984999999pt height=22.465723500000017pt/> if <img src="svgs//77a3b857d53fb44e33b53e4c8b68351a.svg?invert_in_darkmode" align=middle width=5.663225699999989pt height=21.68300969999999pt/> and <img src="svgs//36b5afebdba34564d884d347484ac0c7.svg?invert_in_darkmode" align=middle width=7.710416999999989pt height=21.68300969999999pt/> are neighbours on a square lattice, and 0
otherwise.  SIR and SEEIIR models run on this graph (`sir_sq` and `seeiir_sq`)

### Fully connected graph

The SIR model is implemented on a simple FC graph where <img src="svgs//c9b68e17641d236fb3019541ac2f8246.svg?invert_in_darkmode" align=middle width=50.82941984999999pt height=22.465723500000017pt/>
(`sir_fc`).  The SEEIIR model (`seeiir_fc`) is implemented on an FC
graph where the <img src="svgs//3a91172ac0904d6988776a8ee5bd4140.svg?invert_in_darkmode" align=middle width=19.870705799999985pt height=22.465723500000017pt/> are given by

<p align="center"><img src="svgs//6b64bb8a482d9288b7629712256ff178.svg?invert_in_darkmode" align=middle width=120.06764879999999pt height=37.099754999999995pt/></p>

and the <img src="svgs//3d13090ef3ed1448f3c4dc166d06ab4d.svg?invert_in_darkmode" align=middle width=13.948864049999989pt height=22.831056599999986pt/> are drawn from an exponential distribution of mean
<img src="svgs//c9b542784b9c76101bedbbaf0ec5da18.svg?invert_in_darkmode" align=middle width=22.950967049999992pt height=24.65753399999998pt/>.  See the [parameter file](./seeiir_fc_par.dat).
