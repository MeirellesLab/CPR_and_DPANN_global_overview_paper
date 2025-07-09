library(tidyverse)

####################### Clean up the pairwise permanova results for ecosystem ########################################

# Load the file as fixed-width (since spacing is inconsistent)
txt_eco <- readLines("data/statistics/pairwise_permanova_ecosystem.txt")

# Extract headers
header1 <- txt_eco[1]
header2 <- txt_eco[38]  # second header, can be ignored or used for validation

# Extract numbered lines for main and continuation parts
main_part <- txt_eco[2:37]
tail_part <- txt_eco[39:length(txt_eco)]

# Parse main part (lines 2–37)
main_df <- lapply(main_part, function(line) {
  # Extract row number
  row_id <- as.integer(sub("^\\s*(\\d+).*", "\\1", line))
  # Remove the number prefix
  clean_line <- sub("^\\s*\\d+\\s+", "", line)
  # Match pattern: extract text before the first numeric block (Df)
  split_pos <- regexpr("\\s+\\d+\\s+", clean_line, perl = TRUE)
  pair <- substr(clean_line, 1, split_pos - 1)
  rest <- substr(clean_line, split_pos, nchar(clean_line))
  rest <- gsub(" +", "\t", trimws(rest))
  fields <- unlist(strsplit(rest, "\t"))
  return(c(id = row_id, pair, fields))
})

main_df <- as.data.frame(do.call(rbind, main_df), stringsAsFactors = FALSE)

# Parse tail part (lines 39–74)
tail_df <- lapply(tail_part, function(line) {
  row_id <- as.integer(sub("^\\s*(\\d+).*", "\\1", line))
  clean_line <- sub("^\\s*\\d+\\s+", "", line)
  fields <- unlist(strsplit(gsub(" +", "\t", trimws(clean_line)), "\t"))
  return(c(id = row_id, fields))
})

tail_df <- as.data.frame(do.call(rbind, tail_df), stringsAsFactors = FALSE)

# Merge both parts by row number
colnames(main_df) <- c("id", "Comparison", "Df", "SumsOfSqs", "F.Model", "R2", "p.value")
colnames(tail_df) <- c("id", "p.adjusted", "sig")

# Merge
final_df_eco <- merge(main_df, tail_df, by = "id")
final_df_eco <- final_df_eco[order(as.numeric(final_df_eco$id)), ]  # sort by row number
final_df_eco$id <- NULL  # remove id column

# Convert numeric columns
num_cols <- c("Df", "SumsOfSqs", "F.Model", "R2", "p.value", "p.adjusted")
final_df_eco[num_cols] <- lapply(final_df_eco[num_cols], as.numeric)

# Export to CSV
write.csv(final_df_eco, "data/statistics/pairwise_permanova_ecosystem_clean.txt", row.names = FALSE)

#Exlude objects exept final_df_eco
rm(list = setdiff(ls(), "final_df_eco"))
