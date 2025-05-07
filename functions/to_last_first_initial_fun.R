
# FUNCTION TO EXTRACT LAST NAME AND FIRST INITIAL OF AUTHORS ###################

to_last_first_initial_fun <- function(author_vec) {
  
  sapply(author_vec, function(name) {
    
    parts <- unlist(strsplit(trimws(name), ",\\s*"))
    
    if (length(parts) >= 2) {
      
      paste0(tolower(parts[1]), ",", substr(tolower(parts[2]), 1, 1))
      
    } else {
      
      tolower(name)
      
    }
  })
}