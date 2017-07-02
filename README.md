# DBDA2Estan
This is a collection of models from Doing Bayesian Data Analysis, 2nd Edition implemented in Stan.  I have attempted to keep the implementation as close as possible to the JAGS versions developed by John Kruschke.  In most cases this was a straightforward translation of syntax.  Models that relied on discrete parameters were implemented by marginalizing over the discrete parameters and included in the code as a custom probability distribution (see the Stan reference manual for details).  

Please feel free to contact me with any questions or comments, either via github or at joseph.houpt@wright.edu
