.onAttach <- function(libname, pkgname) {
  # Define pieces and empty chessboard
  pieces <- c("\u265C", "\u265E", "\u265D", "\u265B", "\u265A", "\u265D", "\u265E", "\u265C",
              "\u265F", "\u265F", "\u265F", "\u265F", "\u265F", "\u265F", "\u265F", "\u265F",
              "\u2656", "\u2658", "\u2657", "\u2655", "\u2654", "\u2657", "\u2658", "\u2656",
              "\u2659", "\u2659", "\u2659", "\u2659", "\u2659", "\u2659", "\u2659", "\u2659")
  chessboard <- matrix("  ", nrow = 8, ncol = 8)
  
  # Assign pieces to random positions
  positions <- sample(1:64)
  chessboard[positions] <- pieces
  
  # Create the string that represents the chessboard
  row <- paste(rep("|---------------------", 1), "|", sep = "")
  chessboard_str <- ""
  for (i in 1:8) {
    chessboard_str <- paste0(chessboard_str, row, "\n")
    chessboard_str <- paste0(chessboard_str, "|", paste(chessboard[i, ], collapse = "|"), "|\n")
  }
  chessboard_str <- paste0(chessboard_str, row)
  
  packageStartupMessage(paste("
                                                
", chessboard_str, "

Flavio Gioia welcomes you to cmahalanobis!
Type 'citation(\"cmahalanobis\")' for citing this R package in publications.
"))
}
