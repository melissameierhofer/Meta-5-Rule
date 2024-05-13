## function for 5% rule

semi_colon_splitter <- function(input1, input2, names = c("input1","input2")){
  
  df        <- data.frame(input1,input2)  
  df$input1 <- as.factor(df$input1)
  df$input2 <- as.factor(df$input2)
  
  to_separate <- levels(df$input1)[grepl(";", levels(df$input1))]
  
  df_all <- df[df$input1 %in% to_separate ,]
  df     <- df[!df$input1 %in% to_separate,]
  df$input1 <- droplevels(df$input1)
  
  df_all$input1 <- as.character(df_all$input1)
  
  for(i in nrow(df_all)) {
    
    df_i <- df_all[i,]
    split   <- strsplit(df_all$input1, ";")[[1]]
    split   <- trimws(split, which = c("both"))
    
    df <- rbind(df,data.frame(input1  = split,
                              input2  = rep(df_i$input2, length(split))))
    
  }
  
  colnames(df) <- names
  return(df)
}


