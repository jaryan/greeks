.onLoad <- function(libname, pkgname) {
#  library.dynam("greeks", pkgname, libname)
}

# temp names for .Call BSput and BScall
GreeksNames <- c("value", "delta", "gamma", "vega", "theta", "rho", "vanna", 
"charm", "zomma", "speed", "colour", "DvegaDtime", "vomma", "dualdelta", 
"dualgamma")
