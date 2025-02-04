#'@title Permnova analysis
#'@description This script performs the peranova analysis and get letters for
#'the pairwise comparisons.
#'@param data Dataframe containing the data to be analyzed.
#'@param response Response variable.
#'@param full_model Full model to be used in the analysis.
#'@param predictor_letters Predictor variable to be used in the pairwise
#'comparisons.
#'@param file_name Name of the file to be saved.
#'@param method Method to be used in the analysis.
#'@author Bright Mage

library(vegan)
library(multcomp)
library(lmPerm)

do_peranova <-
  function(
      data,
      response,
      full_model,
      predictor_letters,
      file_name,
      method = "jaccard",
      parallel = 1,
      permutations = 10,
      posterior = TRUE) {

    formula <- as.formula(paste0("data[[response]] ~ ", full_model))

    if (!file.exists(file_name)) {
      print(paste0(file_name, " permanova not found!, running..."))
      permanova <- adonis2(
        formula,                               data = data,
        permutations = permutations,   parallel = parallel,
        method = method
      )
      capture.output(
        permanova, file = file_name
      )
      if (posterior == TRUE) {
        # Getting pairwise comparisons
        model_aov <- aovp(
          data[[response]] ~ data[[predictor_letters]],
          perm = "",  data = data
        )

        model_pairwise <- glht(
          model_aov, linfct = mcp("data[[predictor_letters]]" = "Tukey")
        )
        tuk_cld <- cld(model_pairwise)
        capture.output(
          tuk_cld$mcletters$Letters,
          file = paste0(file_name, "_letters")
        )
      }


   } else {
     print(paste0("Output for:", file_name,  "already exists!"))
   }


}



# backup
# community_distancematrix <- vegdist(
#   x = subset(microgroups_richness, microGroup == "Bonafide")$richness,
#   method = "euclidian"
# )
# model_aov <- pairwise.adonis(
#   x = community_distancematrix,
#   factors = subset(microgroups_richness, microGroup == "Bonafide")$ecosystem,
#   perm = 4999,
#   p.adjust.m = "fdr"
# )
