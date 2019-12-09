#' Wrapper for coxph
#'
#' Perform cox regression using coxph function from survival package
#' @param response The response variable
#' @param time Time variable, e.g. time to treatment, time to last known alive...
#' @param endpoint Event variable, e.g. treatment, died...
#' @param scale Whether to center and scale the response variable
#' @export
#' @import survival
#'
#'
com <- function(response, time, endpoint, scale =FALSE) {

  if (scale) {
    #calculate z-score
    response <- (response - mean(response, na.rm = TRUE))/sd(response, na.rm=TRUE)
  }
  surv <- coxph(Surv(time, endpoint) ~ response)


  tibble(p = summary(surv)[[7]][,5],
         HR = summary(surv)[[7]][,2],
         lower = summary(surv)[[8]][,3],
         higher = summary(surv)[[8]][,4])
}


#Internal function to format floats
formatNum <- function(i, limit = 0.01, digits =1, format="e") {
  r <- sapply(i, function(n) {
    if (n < limit) {
      formatC(n, digits = digits, format = format)
    } else {
      format(n, digits = digits)
    }
  })
  return(r)
}

#' Wrapper for plotting Kaplan-Meier plot
#'
#' Using ggsurvplot function from survminer to plot Kaplan-Meier plot
#' @param response The response variable
#' @param time Time variable, e.g. time to treatment, time to last known alive...
#' @param endpoint Event variable, e.g. treatment, died...
#' @param titlePlot Plot title.
#' @param pval The p value that needs to be shown in the plot. If NULL, the p value will be calculated automatically.
#' @param stat The method to stratify samples. If the reponse variable is catagorical, only "binary" should be used. If the response variable is continous, either "median" or "maxstat" can be used.
#' @param maxTime Truncate the maximal time duration.
#' @param showP Whether to show P value in the plot.
#' @param showTable Whether to show the "number at risks" table under the KM plot.
#' @param ylab Label for y axis.
#' @param xlab Label for x axis.
#' @param table_ration The ration of height between the KM plot and the table below.
#' @export
#' @import survminer cowplot maxstat
#'
#'
km <- function(response, time, endpoint, titlePlot = "KM plot", pval = NULL,
               stat = "median", maxTime =NULL, showP = TRUE, showTable = FALSE,
               ylab = "Fraction", xlab = "Time (years)",
               table_ratio = c(0.7,0.3)) {


  #function for km plot
  survS <- tibble(time = time,
                  endpoint = endpoint)

  if (!is.null(maxTime))
    survS <- mutate(survS, endpoint = ifelse(time > maxTime, FALSE, endpoint),
                    time = ifelse(time > maxTime, maxTime, time))

  if (stat == "maxstat") {
    ms <- maxstat.test(Surv(time, endpoint)  ~ response,
                       data = survS,
                       smethod = "LogRank",
                       minprop = 0.2,
                       maxprop = 0.8,
                       alpha = NULL)

    survS$group <- factor(ifelse(response >= ms$estimate, "high", "low"))
    p <- com(survS$group, survS$time, survS$endpoint)$p

  } else if (stat == "median") {
    med <- median(response, na.rm = TRUE)
    survS$group <- factor(ifelse(response >= med, "high", "low"))
    p <- com(survS$group, survS$time, survS$endpoint)$p

  } else if (stat == "binary") {
    survS$group <- factor(response)
    if (nlevels(survS$group) > 2) {
      sdf <- survdiff(Surv(survS$time,survS$endpoint) ~ survS$group)
      p <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
    } else {
      p <- com(survS$group, survS$time, survS$endpoint)$p
    }
  }

  if (is.null(pval)) {
    if(p< 1e-16) {
      pAnno <- bquote(italic("P")~"< 1e-16")
    } else {
      pval <- formatNum(p, digits = 1)
      pAnno <- bquote(italic("P")~"="~.(pval))
    }

  } else {
    pval <- formatNum(pval, digits = 1)
    pAnno <- bquote(italic("P")~"="~.(pval))
  }

  if (!showP) pAnno <- ""

  p <- ggsurvplot(survfit(Surv(time, endpoint) ~ group, data = survS),
                  data = survS, pval = FALSE,  conf.int = FALSE,
                  legend = ifelse(showTable, "none","top"),
                  ylab = "Fraction", xlab = "Time (years)", title = titlePlot,
                  pval.coord = c(0,0.1), risk.table = showTable, legend.labs = sort(unique(survS$group)),
                  ggtheme = theme_classic() + theme(plot.title = element_text(hjust =0.5), panel.border = element_blank()))
  if (!showTable) {
    p <- p$plot + annotate("text",label=pAnno, x = 0.1, y=0.1, hjust =0, size =5)
    return(p)
  } else {
    #construct a gtable
    pp <- p$plot + annotate("text",label=pAnno, x = 0.1, y=0.1, hjust =0, size=5)
    pt <- p$table + ylab("") + xlab("") + theme(plot.title = element_text(hjust=0, size =10))
    p <- plot_grid(pp,pt, rel_heights = table_ratio, nrow =2, align = "v")
    return(p)
  }
}

#' Multivariate Cox regression
#'
#' Wrapper function to run multivariate Cox regression using coxph
#' @param survTab A tidy table contain the survival data in three columns: identifier, time, event
#' @param riskTab A tidy table that contain the risk factors that need to be tested. One column must contain the sample identifiers. All other columns will be treated as risk factors and tested in the multivariate Cox model.
#' @param time The column name of the time variable in the survTab.
#' @param endpoint The column name of the event variable in the survTab.
#' @param id The column name of the sample identifier, which is used to join the two table.
#' @export
#' @import survival
runCox <- function(survTab, riskTab, time, endpoint, id = "patientID") {
  survTab <- select(survTab, !!id, !!time, !!endpoint) %>%
    dplyr::rename(time = !!time, endpoint = !!endpoint) %>%
    filter(!is.na(time), !is.na(endpoint))
  testTab <- right_join(survTab, riskTab, by = id) %>%
    select(-!!id)
  surv1 <- coxph(
    Surv(time, endpoint) ~
      .,
    data = testTab)
  return(surv1)
}

#' Forest plot for multivariate Cox regression result.
#' Plot the hazard ration and confidence interval for each risk factors from a multivariate Cox regression results.
#' @param survRes The output of runCox function
#' @param title The plot for the title.
#' @export
#' @import survival
plotHazard <- function(survRes, title = "") {
  sumTab <- summary(survRes)$coefficients
  confTab <- summary(survRes)$conf.int
  #correct feature name
  nameOri <- rownames(sumTab)
  nameMod <- substr(nameOri, 1, nchar(nameOri) -1)
  plotTab <- tibble(feature = rownames(sumTab),
                    nameMod = substr(nameOri, 1, nchar(nameOri) -1),
                    HR = sumTab[,2],
                    p = sumTab[,5],
                    Upper = confTab[,4],
                    Lower = confTab[,3]) %>%
    mutate(feature = ifelse(nameMod %in% names(survRes$xlevels), nameMod, feature)) %>%
    #mutate(feature = str_replace(feature, "[.]","/")) %>%
    #mutate(feature = str_replace(feature, "[_]","-")) %>%
    arrange(desc(abs(p))) %>% mutate(feature = factor(feature, levels = feature)) %>%
    mutate(type = ifelse(HR >1 ,"up","down")) %>%
    mutate(Upper = ifelse(Upper > 10, 10, Upper))

  ggplot(plotTab, aes(x=feature, y = HR, color = type)) +
    geom_hline(yintercept = 1, linetype = "dotted", color = "grey50") +
    geom_point(position = position_dodge(width=0.8), size=3, color = "black") +
    geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.3, size=1,color = "grey20") +
    geom_text(position = position_nudge(x = 0.3),
              aes(y = HR, label =  sprintf("italic(P)~'='~'%s'",
                                           formatNum(p, digits = 1))),
              color = "black", size =5, parse = TRUE) +
    expand_limits(y=c(-0.5,0))+
    ggtitle(title) + scale_y_log10() +
    ylab("Hazard ratio") +
    coord_flip() +
    theme_bw() +
    theme(legend.position = "none", axis.title.y = element_blank(), panel.grid = element_blank())
}
