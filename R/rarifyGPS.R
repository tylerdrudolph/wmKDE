##########################################################################################################################################
## Rar?faction d'une base de donnees t?l?m?triques afin de produire une observation par individu par jour (r?f?rence = midi)

rarifyGPS = function(df, date.heure=NULL, id=NULL, reftime = chron::times(c("12:00:00"))) {

  if(is.null(date.heure)) stop("Il faut fournir un champs 'date-heure'")
  if(!is.null(date.heure) & !inherits(date.heure, "POSIXct")) stop("date.heure doit relever de la classe 'POSIXct'")

  df$date.heure <- date.heure

  if(is.null(id)) stop("Il faut fournir un champs 'id'")
  df$id <- id

  ## S'assurer d'avoir un champs intitul? "time" ayant le format 'times' (package::chron)
  df$time = chron::times(format(df$date.heure, format="%H:%M:%S"))

  ## Soustraire toutes les dates-heure de midi afin d'avoir les valeurs m?dianes
  df$diff <- abs(reftime - df$time)
  df$diff <- as.numeric(df$diff)

  ## Trier les valeurs par id, ensuite jour, et finalement diff
  df <- df[order(df$id, df$date.heure, df$diff), ]

  ## Extraire les observations correspondant au critère de base
  df.min <- stats::aggregate(df$diff, list(id = df$id, date = as.Date(df$date.heure)), min, na.rm = TRUE)
    colnames(df.min) <- c("id", "date", "diff")  # renommer les colonnes
  df.subset <- subset(df, paste(df$id, as.Date(df$date.heure), df$diff) %in% paste(df.min$id, df.min$date, df.min$diff))

  ## Identifier et éliminer les exceptions s'agissant de plus d'une observation par jour/id
  df.subset.tab <- data.frame(table(df.subset$id, as.Date(df.subset$date.heure)))
  colnames(df.subset.tab) <- c("id", "date", "Freq")
  df.subset.tab <- subset(df.subset.tab, df.subset.tab$Freq > 1)

  ## Éliminer les doublons
  if(nrow(df.subset.tab) > 0) {
    df.subset <- df.subset[!duplicated(paste(df.subset$id, df.subset$date.heure)),]
    # system.time(x <- lapply(1:50, function(m) which(df.subset$id==df.subset.tab$id[m] & as.Date(df.subset$date.heure)==as.character(df.subset.tab$date[m]))))
    #
    #
    # ## séparer les doublons
    # #dupes <- df.subset[paste(df.subset$id, as.Date(df.subset$date.heure)) %in% paste(df.subset.tab$id, df.subset.tab$date),] #duplicated(paste(df.subset$id, as.Date(df.subset$date.heure)), fromLast=TRUE), ]
    # #df.subset <- df.subset[!paste(df.subset$id, as.Date(df.subset$date.heure)) %in% paste(df.subset.tab$id, df.subset.tab$date),]
    # # tirer une localisation aléatoire parmi celles disponibles
    # for(m in 1:nrow(df.subset.tab)) {
    #   cat(m, '/', nrow(df.subset.tab), str_c('(', round(m/nrow(df.subset.tab), digits=3), '%)'), '\n')
    #   df.subset <- filter(df.subset, id==df.subset.tab$id[m] & as.Date(date.heure)==as.character(df.subset.tab$date[m]))
    #
    #                      slice(filter(dupes, id==df.subset.tab$id[m] & as.Date(date.heure)==as.character(df.subset.tab$date[m])), sample.int(n=df.subset.tab$Freq[m], size=1)))
    #
    # }
  }

  return(df.subset)

}

