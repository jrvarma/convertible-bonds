# Numerical valuation of convertible bonds and options on stocks. 

The software uses a binomial lattice with the stock price as the only state variable to value convertible bonds with call and put features.  The software does *not* use the warrant valuation approach which requires the volatility of equity (stocks plus warrants).  Instead, it ignores the dilution effect and uses stock price volatility which is more readily available. 

The software can also be used in valuing options on stocks. This is more useful for American puts, and to a lesser extent American calls on dividend paying stocks. European options can be valued more easily by the Black-Scholes formula. 

The software was written in 2003 and has not been actively used and maintained for many years now. It was last revised in 2013 to remove compiler warnings and errors while using a more recent version of `gcc`.
