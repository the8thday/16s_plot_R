#! /urs/bin/env Rscript
# 生存曲线

library(survival)
library(survminer)

fit<- survfit(Surv(time, status) ~ sex, data = lung)
res.sum <- surv_summary(fit)
ggsurv <- ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  data = lung,             # data used to fit survival curves.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimates of survival curves.
  palette = c("#E7B800", "#2E9FDF"),
  xlim = c(0,500),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in days",   # customize X axis label.
  break.time.by = 100,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T,# colour risk table text annotations.
  risk.table.height = 0.25, # the height of the risk table
  risk.table.y.text = FALSE,# show bars instead of names in text annotations
  # in legend of risk table.
  ncensor.plot = TRUE,      # plot the number of censored subjects at time t
  ncensor.plot.height = 0.25,
  conf.int.style = "step",  # customize style of confidence intervals
  surv.median.line = "hv",  # add the median survival pointer.
  legend.labs = 
    c("Male", "Female")    # change legend labels.
)
ggsurv
res$table <- res$table + theme(axis.line = element_blank())
res$plot <- res$plot + labs(title = "Survival Curves")

surv_diff <- survdiff(Surv(time, status) ~ sex, data = lung)
# cox比例风险模型 ---------------------------------------------------------------
# 具体某一因素对生存的影响程度,exp(coef)为HR

reg <- coxph(Surv(time, status) ~ sex, data = lung)
survminer::ggforest(reg, data = lung)

fit2 <- coxph(Surv(time, status)~sex+age+ph.ecog+ph.karno+pat.karno+meal.cal+wt.loss, data=lung)
survminer::ggforest(fit2, data = lung)


# 处理年龄 --------------------------------------------------------------------
surv_cutpoint(lung$age)
lung$agecat <- cut(lung$age, breaks=c(0, 70, Inf), labels=c("young", "old"))
lung <- lung %>% 
  mutate(agecat=cut(age, breaks=c(0, 70, Inf), labels=c("young", "old")))
ggsurvplot(survfit(Surv(time, status)~agecat, data=lung), pval=TRUE)
