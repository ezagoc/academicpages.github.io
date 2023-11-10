
library(didimputation)
library(tidyverse)
library(lfe)
library(DIDmultiplegt)
library(did)

# Functions:

# Create a tidier for "multiplegt" objects
tidy.event_ch_dh = function(x, level = 0.95) {
  ests = x[grepl("^placebo_|^effect|^dynamic_", names(x))]
  ret = data.frame(
    term      = names(ests),
    estimate  = as.numeric(ests),
    std.error = as.numeric(x[grepl("^se_placebo|^se_effect|^se_dynamic", names(x))]),
    N         = as.numeric(x[grepl("^N_placebo|^N_effect|^N_dynamic", names(x))])
  ) |>
    # For CIs we'll assume standard normal distribution
    within({
      conf.low  = estimate - std.error*(qnorm(1-(1-level)/2))
      conf.high = estimate + std.error*(qnorm(1-(1-level)/2))
    })
  return(ret)
}

# Data from Callaway and Sant'Anna (2021)

data(mpdta)

mpdta <- mpdta |> mutate(post = ifelse(year >= first.treat, 1, 0),
                         post = ifelse(first.treat == 0, 0, post),
                         post_treat = post*treat)

twfe <- felm(lemp ~ post_treat | countyreal + year | 0 | countyreal,
             data = mpdta)

coef <- c(twfe$coefficients)
se_twfe <- c(twfe$se)

df_twfe <- tibble(type = "TWFE", coef = coef, se = se_twfe, 
                  conf.low = coef - 1.96*se_twfe, 
                  conf.high = coef + 1.96*se_twfe)

borus <- did_imputation(data = mpdta, yname = "lemp", gname = "first.treat",
                        first_stage = ~ 0 | countyreal + year,
                        tname = "year", idname = "countyreal") |> 
  mutate(type = "Borusyak") |> 
  select(type, coef = estimate, se = std.error, conf.low, conf.high)

ch_dh <- did_multiplegt(df = mpdta, Y = "lemp", G = "countyreal", T = "year", 
                        D = "post_treat", brep = 40, 
                        cluster   = 'countyreal')

coef_ch_dh <- ch_dh$effect

se_ch_dh <- ch_dh$se_effect

df_ch_dh <- tibble(type = "Ch & D'H", 
                   coef = coef_ch_dh, 
                   se = se_ch_dh, 
                   conf.low = coef - 1.96*se, 
                   conf.high = coef + 1.96*se)

results_full_panel <- rbind(df_twfe, borus, df_ch_dh) |> mutate(panel_type = 'full')


### Event Study:

mpdta <- mpdta %>% group_by(countyreal) %>%
  mutate(time_to_event = year - first.treat) |> 
  mutate(time_to_event2 = ifelse(treat == 0, 0, time_to_event))

dd_plot_reg <- feols(lemp ~ i(time_to_event2, treat, ref = -1) | 
                     year + countyreal, cluster = "countyreal", data = mpdta)

time <- c(-4:3)
coefs_twfe <- c(dd_plot_reg$coefficients[1:3], 0, dd_plot_reg$coefficients[4:7])
se_twfe <- c(dd_plot_reg$se[1:3], 0, dd_plot_reg$se[4:7])


iplot(dd_plot_reg)

df_event <- tibble(type = 'Canonical', time = time, 
                   coef = coefs_twfe, se = se_twfe, 
                   conf.low = coef - 1.96*se, 
                   conf.high = coef + 1.96*se)


## Borus:

event_borus <- did_imputation(data = mpdta, yname = "lemp", gname = "first.treat",
                              first_stage = ~ 0 | year + countyreal,
                              tname = "year", idname = "countyreal",
                              pretrends = TRUE, horizon = TRUE) |> 
  mutate(time = as.numeric(term), type = 'Borusyak') |> 
  rename(coef = estimate, se = std.error) |> select(-c(lhs, term))

# Un poco de manipulación para generar un df que podamos graficar:

aux_1 <- event_borus |> filter(time < 0)
aux_2 <- event_borus |> filter(time >= 0)
aux_3 <- tibble(time = -1, se = 0, coef = 0,
                conf.low = 0, conf.high = 0, type = "Borusyak")

event_borus1 <- rbind(aux_1, aux_3, aux_2)

## CHDH

library(broom)

event_ch_dh <- did_multiplegt(df = mpdta, Y = "lemp", G = "countyreal", 
                              T = "year", 
                              D = "post_treat", brep = 21,
                              cluster = 'countyreal', dynamic = 3, placebo = 3)
# Tidy frame:

tidy_ch_dh = tidy.event_ch_dh(event_ch_dh) |> within({
  term = gsub("^placebo_", "-", term)
  term = gsub("^effect", "0", term)
  term = gsub("^dynamic_", "", term)
  term = as.integer(term)
}) |> select(time = term, coef = estimate, se = std.error, conf.low, conf.high) |> 
  mutate(type = "Ch & D'H")


### Graphs:

event_can <- rbind(df_event |> mutate(time = time - .15,
                                      type_model = 'TWFE'), 
                   event_borus1 |> mutate(type_model = 'Borusyak'),
                   tidy_ch_dh |> mutate(time = time + .15,
                                         type_model = "Ch & D'H"))

colors <- c("#000000", "#0072B2","#D55E00")
borus_plot <- ggplot(data = event_can, mapping = aes(y = coef, x = time)) +
  geom_point(aes(colour = factor(type_model)), size = 2) + 
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high, colour = factor(type_model)), width=0.2) +
  geom_hline(yintercept = 0, linetype="solid", color ="grey", 2) +
  geom_vline(xintercept = 0, linetype="dashed", color ="red", 2) +
  theme_bw() +
  ylab("Estimated Coefficient (95% IC)") + 
  xlab("Time since Treatment") + 
  scale_color_manual(name = "Model", values= colors) +
  theme(legend.position = "right")
borus_plot

ggsave(borus_plot, filename = 'G_Event_Canonical.pdf', device = cairo_pdf,
       dpi = 300, width = 18, height = 13, units = 'cm')