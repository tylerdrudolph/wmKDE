#' Subset a GPS collar telemetry database according to user preferences and/or analytical objectives.
#'
#' @param df data.frame
#' @param pop factor
#' @param saison factor
#' @param sex factor
#' @param annee integer
#' @param id character
#' @param nmin vector
#'
#' @return original df with specified filters applied
#' @export
#'
subsetDB <- function(df, pop=NULL, saison=NULL, sex=NULL, annee=NULL,
                     id=NULL, nmin=50) {

  ## Population(s) of interest
  if(shiny::isTruthy(pop)) {
    pop <- match.arg(pop, choices=as.character(unique(df$Pop)), several.ok=TRUE)
    df <- df[df$Pop==pop,]
    if(is.factor(df$Pop)) df$Pop <- droplevels(df$Pop)
  }

  ## Season(s) of interest
  if(shiny::isTruthy(saison)) {
    saison <- match.arg(saison, choices=as.character(unique(df$saison)), several.ok = T)
    df <- df[df$saison %in% saison,]
    if(is.factor(df$saison)) df$saison <- droplevels(df$saison)
  }

  ## Sex(s) of interest
  if(shiny::isTruthy(sex)) {
    sex <- match.arg(sex, choices=c('F','M'), several.ok = T)
    df <- df[df$Sex %in% sex,]
    if(is.factor(df$Sex)) df$Sex <- droplevels(df$Sex)
  }

  ## Year(s) of interest
  if(shiny::isTruthy(annee)) {
    annee <- as.integer(match.arg(as.character(annee), choices=as.character(unique(df$An)), several.ok = T))
    df <- df[df$An %in% annee,]
    if(is.factor(df$An)) df$An <- droplevels(df$An)
  }

  ## Individual(s) of interest
  if(shiny::isTruthy(id)) {
    id <- match.arg(id, choices=as.character(unique(df$IDAnimal)), several.ok = T)
    df <- df[df$IDAnimal %in% id,]
    if(is.factor(df$IDAnimal)) df$IDAnimal <- droplevels(df$IDAnimal)
  }

  ## Exclude individuals with insufficient sample sizes
  ntab <- table(df$IDAnimal)
  df <- df[df$IDAnimal %in% names(ntab)[ntab >= nmin], ]
  if(is.factor(df$IDAnimal)) df$IDAnimal <- droplevels(df$IDAnimal)

  methods::slot(df, 'data') <- droplevels(methods::slot(df, 'data'))

  return(df)

}

