welcome<- function(){
msg <- c(paste0(
"Welcome to
 _____                                                                            ______________________
|     |                                                                           |                    |
|     |                                                                           |                    |
|     |                                                                           |      _________     |
|     |                                                                           |     |        |     |
|     |                                                                           |     |        |     |
|     |                     ____________________       ____________________       |     |________|     |
|     |                    |                    |     |                    |      |                    |
|     |______________      |                    |     |                    |      |                    |
|                    |     |      ______________|     |      ______________|      |      ______________|
|                    |     |     |                    |     |                     |     |
|      _________     |     |     |                    |     |                     |     |
|     |        |     |     |     |                    |     |                     |     |
|     |        |     |     |     |                    |     |                     |     |
|     |________|     |     |     |______________      |     |______________       |     |
|                    |     |                    |     |                    |      |     |
|                    |     |                    |     |                    |      |     |
|____________________|     |____________________|     |____________________|      |_____|                 version ",

packageVersion("bccp")),"\nType 'citation(\"bccp\")' for citing this R package in publications.")
 return(msg)
}
.onAttach <- function(libname, pkgname) {
  mess <- welcome()
  if(!interactive())
    mess[1] <- paste("Package 'bccp' version", packageVersion("bccp"))
    packageStartupMessage(mess)
  invisible()
  }

