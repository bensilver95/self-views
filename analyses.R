library(tidyverse)
library(lmerTest)
library(brms)
library(corrplot)
setwd("/Users/bensilver/Google Drive/My Drive/Grad School/Projects/Self_Memory/")

# load raw data
pre <- read_csv("data/pre_new.csv")
post1 <- read_csv("data/post1_new.csv")
post2 <- read_csv("data/post2_new.csv")
post3 <- read_csv("data/post3_new.csv")

##### Q1 - Traits #####
#### Clean/prep #####
pre_traits <- pre %>% 
  select(c("prolificid","age","gender",
           "lifesat_2","po_spectrum_scale_1","po_engaged_1","change_traits_5",
           starts_with("traits"))) %>% 
  rename(change_traits_pre = change_traits_5) %>% 
  pivot_longer(cols = starts_with("traits"),
               names_to = "trait",
               names_prefix = ("traits_"),
               values_to = "pre") %>% 
  rename(lifesat = lifesat_2,
         pospectrum = po_spectrum_scale_1,
         engaged = po_engaged_1)

post1_traits <- post1 %>% 
  select(c("prolificid", "change_traits_1",
           starts_with("traits"))) %>% 
  rename(change_traits_post1 = change_traits_1) %>% 
  pivot_longer(cols = starts_with("traits"),
               names_to = "trait",
               names_prefix = ("traits_"),
               values_to = "post1")

post2_traits <- post2 %>% 
  select(c("prolificid","change_traits_1",
           starts_with("traits"))) %>% 
  rename(change_traits_post2 = change_traits_1) %>% 
  pivot_longer(cols = starts_with("traits"),
               names_to = "trait",
               names_prefix = ("traits_"),
               values_to = "post2")

post3_traits <- post3 %>% 
  select(c("prolificid","change_traits_1",
           starts_with("traits"))) %>% 
  rename(change_traits_post3 = change_traits_1) %>% 
  pivot_longer(cols = starts_with("traits"),
               names_to = "trait",
               names_prefix = ("traits_"),
               values_to = "post3")

traits <- pre_traits %>% 
  full_join(post1_traits) %>% 
  full_join(post2_traits) %>% 
  full_join(post3_traits) %>% 
  mutate(post1_pre = abs(post1 - pre),
         post2_pre = abs(post2 - pre),
         post3_pre = abs(post3 - pre),
         post2_post1 = abs(post2 - post1),
         post3_post1 = abs(post3 - post1),
         post3_post2 = abs(post3 - post2),
         trait = as.integer(trait),
         positive = if_else(trait < 11,1,0)) %>% 
  select(-contains("change"))

#### models ########
traits_model <- traits %>% 
  select(-c("pre","post1","post2","post3")) %>% 
  pivot_longer(cols = contains("_"),
               names_to = "comparison",
               values_to = "difference") %>% 
  mutate(pre_comp = if_else(str_detect(comparison,"pre"),1,0))

# full model
summary(lmer(difference ~ pre_comp + 
               (1 + pre_comp | prolificid),
             data = traits_model))

# single comparisons
summary(lmer(difference ~ pre_comp + 
               (1 + pre_comp | prolificid),
             data = traits_model %>% 
               filter(comparison == "post3_pre" | 
                        comparison == "post3_post1")
))

# political spectrum
summary(lmer(difference ~ pre_comp*pospectrum + 
               (1 + pre_comp | prolificid),
             data = traits_model %>% 
               filter(comparison == "post1_pre" | 
                        comparison == "post2_post1")))

## non-absolute value to test direction of movement
extremity_traits <- traits %>% 
  mutate(post1_pre_nonabs = post1 - pre,
         post2_pre_nonabs = post2 - pre,
         post3_pre_nonabs = post3 - pre,
         post2_post1_nonabs = post2 - post1,
         post3_post1_nonabs = post3 - post1,
         post3_post2_nonabs = post3 - post2)

extremity_traits_model <- extremity_traits %>% 
  select(-c("pre","post1","post2","post3")) %>% 
  pivot_longer(cols = contains("nonabs"),
               names_to = "comparison",
               values_to = "difference") %>% 
  mutate(pre_comp = if_else(str_detect(comparison,"pre"),1,0))

# basic model
summary(lmer(difference ~ pre_comp*positive + 
               (1 + pre_comp | prolificid),
             data = extremity_traits_model))

# single comparisons
summary(lmer(difference ~ pre_comp*positive + 
               (1 + pre_comp | prolificid),
             data = extremity_traits_model %>% 
               filter(comparison == "post1_pre_nonabs" | 
                        comparison == "post2_post1_nonabs")))

#### figures #######

traits_model$comparison <- factor(traits_model$comparison,
                                  levels = c("post1_pre",
                                             "post2_pre",
                                             "post3_pre",
                                             "post2_post1",
                                             "post3_post1",
                                             "post3_post2"))

traits_figure <- traits_model %>% 
  group_by(prolificid,comparison,pre_comp,pospectrum) %>% 
  summarize(difference = mean(difference)) %>% 
  mutate(comparison = recode(comparison,
                             "post1_pre" = "T1-T0","post2_pre" = "T2-T0",
                             "post3_pre" = "T3-T0","post2_post1" = "T2-T1",
                             "post3_post1" = "T3-T1","post3_post2" = "T3-T2"),
         pre_comp = factor(pre_comp,levels = c(1,0)))

ggplot(traits_figure, 
       aes(x = comparison,y = difference,
           color = pre_comp)) +
  geom_jitter(alpha = .1, size = 2) +
  stat_summary(fun = mean, fun.min = mean, fun.max = mean, 
               width = .9, geom = "crossbar") +
  theme_classic() +
  labs(x = "Comparison Type", y = "Self-View Change\n(Absolute Difference)",
       title = "Changes in self-views of traits", color = "Timepoint\nComparison") +
  scale_color_manual(values = c("steelblue","olivedrab"),
                     labels = c("Across election","After election")) +
  theme(axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 12),
        plot.title = element_text(size = 18, hjust = .5))
ggsave("figs/paper/traits_q1.jpg", scale = .7)
##### Q1 - Values ######
#### Clean/prep #####
pre_values <- pre %>% 
  select(c("prolificid","age","gender",
           "lifesat_2","po_spectrum_scale_1","po_engaged_1","change_values_4",
           starts_with("values"))) %>% 
  rename(change_values_pre = change_values_4) %>% 
  pivot_longer(cols = starts_with("values"),
               names_to = "value",
               names_prefix = ("values_"),
               values_to = "pre") %>% 
  rename(lifesat = lifesat_2,
         pospectrum = po_spectrum_scale_1,
         engaged = po_engaged_1)

post1_values <- post1 %>% 
  select(c("prolificid","change_values_1",
           starts_with("values"))) %>% 
  rename(change_values_post1 = change_values_1) %>% 
  pivot_longer(cols = starts_with("values"),
               names_to = "value",
               names_prefix = ("values_"),
               values_to = "post1") %>% 
  mutate(post1 = scales::rescale(post1, to = c(1,7), from = c(0,7)))

post2_values <- post2 %>% 
  select(c("prolificid","change_values_1",
           starts_with("values"))) %>% 
  rename(change_values_post2 = change_values_1) %>% 
  pivot_longer(cols = starts_with("values"),
               names_to = "value",
               names_prefix = ("values_"),
               values_to = "post2") %>% 
  mutate(post2 = scales::rescale(post2, to = c(1,7), from = c(0,7)))

post3_values <- post3 %>% 
  select(c("prolificid","change_values_1",
           starts_with("values"))) %>% 
  rename(change_values_post3 = change_values_1) %>% 
  pivot_longer(cols = starts_with("values"),
               names_to = "value",
               names_prefix = ("values_"),
               values_to = "post3")

values <- pre_values %>% 
  full_join(post1_values) %>% 
  full_join(post2_values) %>% 
  full_join(post3_values) %>% 
  mutate(post1_pre = abs(post1 - pre),
         post2_pre = abs(post2 - pre),
         post3_pre = abs(post3 - pre),
         post2_post1 = abs(post2 - post1),
         post3_post1 = abs(post3 - post1),
         post3_post2 = abs(post3 - post2),
         value = as.integer(value)) %>% 
  select(-contains("change"))

#### models ######
values_model <- values %>% 
  select(-c("pre","post1","post2","post3")) %>% 
  pivot_longer(cols = contains("_"),
               names_to = "comparison",
               values_to = "difference") %>% 
  mutate(pre_comp = if_else(str_detect(comparison,"pre"),1,0))

# simple model
summary(lmer(difference ~ pre_comp + 
               (1 + pre_comp | prolificid),
             data = values_model))

# single comparisons
summary(lmer(difference ~ pre_comp + 
               (1 + pre_comp | prolificid),
             data = values_model %>% 
               filter(comparison == "post1_pre" | 
                        comparison == "post2_post1")))

# political ideology
summary(lmer(difference ~ pre_comp*pospectrum + 
               (1 + pre_comp | prolificid),
             data = values_model %>% 
               filter(comparison == "post1_pre" | 
                        comparison == "post2_post1")))

## testing correlation between traits and values ####
traits_md <- traits %>% 
  group_by(prolificid) %>% 
  summarize(post1_pre = sum(post1_pre, na.rm = T),
            post2_pre = sum(post2_pre, na.rm = T),
            post3_pre = sum(post3_pre, na.rm = T),
            post2_post1 = sum(post2_post1, na.rm = T),
            post3_post1 = sum(post3_post1, na.rm = T),
            post3_post2 = sum(post3_post1, na.rm = T))

values_md <- values %>% 
  group_by(prolificid) %>% 
  summarize(post1_pre = sum(post1_pre, na.rm = T),
            post2_pre = sum(post2_pre, na.rm = T),
            post3_pre = sum(post3_pre, na.rm = T),
            post2_post1 = sum(post2_post1, na.rm = T),
            post3_post1 = sum(post3_post1, na.rm = T),
            post3_post2 = sum(post3_post1, na.rm = T))

cor.test(traits_md$post1_pre,values_md$post1_pre)
cor.test(traits_md$post2_pre,values_md$post2_pre)
cor.test(traits_md$post3_pre,values_md$post3_pre)

traits_model_md <- traits_model %>% 
  group_by(prolificid) %>% 
  summarize(difference = sum(difference, na.rm = T))

values_model_md <- values_model %>% 
  group_by(prolificid) %>% 
  summarize(difference = sum(difference, na.rm = T))

cor.test(traits_model_md$difference,values_model_md$difference)
# significantly different, strongest for post1_pre and post2_pre

cor.test(traits_model_md$difference,events_sims_model$similarity)


df_cors <- data.frame(traits_md$post1_pre,traits_md$post2_pre,
                      traits_md$post3_pre,traits_md$post2_post1,
                      traits_md$post3_post1,traits_md$post3_post2,
                      values_md$post1_pre,values_md$post2_pre,
                      values_md$post3_pre,values_md$post2_post1,
                      values_md$post3_post1,values_md$post3_post2,
                      events_sims$post1_pre,events_sims$post2_pre,
                      events_sims$post3_pre,events_sims$post2_post1,
                      events_sims$post3_post1,events_sims$post3_post2)

M <- cor(df_cors)

testRes = cor.mtest(df_cors, conf.level = 0.95)
corrplot(M, p.mat = testRes$p, sig.level = .1,  addrect = 2)

#### figures #####

values_model$comparison <- factor(values_model$comparison,
                                  levels = c("post1_pre",
                                             "post2_pre",
                                             "post3_pre",
                                             "post2_post1",
                                             "post3_post1",
                                             "post3_post2"))

values_figure <- values_model %>% 
  group_by(prolificid,comparison,pre_comp,pospectrum) %>% 
  summarize(difference = mean(difference)) %>% 
  mutate(comparison = recode(comparison,
                             "post1_pre" = "T1-T0","post2_pre" = "T2-T0",
                             "post3_pre" = "T3-T0","post2_post1" = "T2-T1",
                             "post3_post1" = "T3-T1","post3_post2" = "T3-T2"),
         pre_comp = factor(pre_comp,levels = c(1,0)))


ggplot(values_figure, 
       aes(x = comparison,y = difference,
           color = as.factor(pre_comp))) +
  geom_jitter(alpha = .1, size = 2) +
  stat_summary(fun = mean, fun.min = mean, fun.max = mean, 
               width = .9, geom = "crossbar") +
  theme_classic() +
  labs(x = "Comparison Type", y = "Self-View Change\n(Absolute Difference)",
       title = "Changes in self-views of values",color = "Timepoint\nComparison") +
  scale_color_manual(values = c("steelblue","olivedrab"),
                     labels = c("Across election","After election")) +
  theme(axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 12),
        plot.title = element_text(size = 18, hjust = .5))
ggsave("figs/paper/values_q1.jpg", scale = .7)
##### Q1 - Events ######
#### Clean/prep #####
## personal, valence #####
pre_events <- pre %>% 
  select(c("prolificid","age","po_spectrum_scale_1",
           contains("unrelated_valence"))) %>% 
  rename("event1" = "unrelated_valence_1",
         "event2" = "unrelated_valence_2",
         "event3" = "unrelated_valence_3") %>% 
  pivot_longer(cols = starts_with("event"),
               names_to = "event",
               values_to = "pre_event_valence")

pre_events_memory <- post2 %>% 
  select(c("prolificid","prev_unr_1_valence_c_1",
           "prev_unr_2_valence_c_1","prev_unr_3_valence_c_1")) %>% 
  rename("event1" = "prev_unr_1_valence_c_1",
         "event2" = "prev_unr_2_valence_c_1",
         "event3" = "prev_unr_3_valence_c_1") %>% 
  pivot_longer(cols = starts_with("event"),
               names_to = "event",
               values_to = "pre_event_memory_valence")

post2_events <- post2 %>% 
  select(c("prolificid",
           contains("unrelated_valence"))) %>% 
  rename("event1" = "unrelated_valence_1",
         "event2" = "unrelated_valence_2",
         "event3" = "unrelated_valence_3") %>% 
  pivot_longer(cols = starts_with("event"),
               names_to = "event",
               values_to = "post2_event_valence")

post2_events_memory <- post3 %>% 
  select(c("prolificid","t4_unrelated_1_val2_1","t4_unrelated_2_val2_1",
           "t4_unrelated_3_val2_1")) %>% 
  rename("event1" = "t4_unrelated_1_val2_1",
         "event2" = "t4_unrelated_2_val2_1",
         "event3" = "t4_unrelated_3_val2_1") %>% 
  pivot_longer(cols = starts_with("event"),
               names_to = "event",
               values_to = "post2_event_memory_valence")

events_personal_valence <- pre_events %>% 
  full_join(pre_events_memory) %>% 
  full_join(post2_events) %>% 
  full_join(post2_events_memory)

## political, valence ##########
pre_events <- pre %>% 
  select(c("prolificid","age","po_spectrum_scale_1",
           starts_with("related_valence"))) %>% 
  rename("event1" = "related_valence_1",
         "event2" = "related_valence_2",
         "event3" = "related_valence_3") %>% 
  pivot_longer(cols = starts_with("event"),
               names_to = "event",
               values_to = "pre_event_valence")

pre_events_memory <- post2 %>% 
  rename("prev_rel_3_valence_p_1" = "prev_rel_event_3_1",
         "prev_rel_3_valence_c_1" = "prev_rel_3_valence_p_1",
         "prev_rel_3_imp_p_1" = "prev_rel_3_valence_c_1",
         "prev_rel_3_imp_c_1" = "prev_rel_3_imp_p_1") %>% 
  select(c("prolificid","prev_rel_1_valence_c_1",
           "prev_rel_2_valence_c_1","prev_rel_3_valence_c_1")) %>% 
  rename("event1" = "prev_rel_1_valence_c_1",
         "event2" = "prev_rel_2_valence_c_1",
         "event3" = "prev_rel_3_valence_c_1") %>% 
  pivot_longer(cols = starts_with("event"),
               names_to = "event",
               values_to = "pre_event_memory_valence")

post2_events <- post2 %>% 
  select(c("prolificid",
           starts_with("related_valence"))) %>% 
  rename("event1" = "related_valence_1",
         "event2" = "related_valence_2",
         "event3" = "related_valence_3") %>% 
  pivot_longer(cols = starts_with("event"),
               names_to = "event",
               values_to = "post2_event_valence")

post2_events_memory <- post3 %>% 
  select(c("prolificid","t4_related_1_val2_1","t4_related_2_val2_1",
           "t4_related_3_val2_1")) %>% 
  rename("event1" = "t4_related_1_val2_1",
         "event2" = "t4_related_2_val2_1",
         "event3" = "t4_related_3_val2_1") %>% 
  pivot_longer(cols = starts_with("event"),
               names_to = "event",
               values_to = "post2_event_memory_valence")

events_political_valence <- pre_events %>% 
  full_join(pre_events_memory) %>% 
  full_join(post2_events) %>% 
  full_join(post2_events_memory)

## personal, importance ########
pre_events <- pre %>% 
  select(c("prolificid","age","po_spectrum_scale_1",
           contains("unrelated_importance"))) %>% 
  rename("event1" = "unrelated_importance_1",
         "event2" = "unrelated_importance_2",
         "event3" = "unrelated_importance_3") %>% 
  pivot_longer(cols = starts_with("event"),
               names_to = "event",
               values_to = "pre_event_importance")

pre_events_memory <- post2 %>% 
  select(c("prolificid","prev_unr_1_imp_c_1",
           "prev_unr_2_imp_c_1","prev_unr_3_imp_c_1")) %>% 
  rename("event1" = "prev_unr_1_imp_c_1",
         "event2" = "prev_unr_2_imp_c_1",
         "event3" = "prev_unr_3_imp_c_1") %>% 
  pivot_longer(cols = starts_with("event"),
               names_to = "event",
               values_to = "pre_event_memory_importance") %>% 
  mutate(pre_event_memory_importance = scales::rescale(pre_event_memory_importance, 
                                                       to = c(1,7), from = c(0,7)))

post2_events <- post2 %>% 
  select(c("prolificid",
           contains("unrelated_importance"))) %>% 
  rename("event1" = "unrelated_importance_1",
         "event2" = "unrelated_importance_2",
         "event3" = "unrelated_importance_3") %>% 
  pivot_longer(cols = starts_with("event"),
               names_to = "event",
               values_to = "post2_event_importance")

post2_events_memory <- post3 %>% 
  select(c("prolificid","t4_unrelated_1_imp2_1","t4_unrelated_2_imp2_1",
           "t4_unrelated_3_imp2_1")) %>% 
  rename("event1" = "t4_unrelated_1_imp2_1",
         "event2" = "t4_unrelated_2_imp2_1",
         "event3" = "t4_unrelated_3_imp2_1") %>% 
  pivot_longer(cols = starts_with("event"),
               names_to = "event",
               values_to = "post2_event_memory_importance") %>% 
  mutate(post2_event_memory_importance = scales::rescale(post2_event_memory_importance, 
                                                         to = c(1,7), from = c(0,7)))

events_personal_importance <- pre_events %>% 
  full_join(pre_events_memory) %>% 
  full_join(post2_events) %>% 
  full_join(post2_events_memory)

## political, importance ##########
pre_events <- pre %>% 
  select(c("prolificid","age","po_spectrum_scale_1",
           starts_with("related_importance"))) %>% 
  rename("event1" = "related_importance_1",
         "event2" = "related_importance_2",
         "event3" = "related_importance_3") %>% 
  pivot_longer(cols = starts_with("event"),
               names_to = "event",
               values_to = "pre_event_importance")

# fixes survey issue
pre_events_memory <- post2 %>% 
  rename("prev_rel_3_valence_p_1" = "prev_rel_event_3_1",
         "prev_rel_3_valence_c_1" = "prev_rel_3_valence_p_1",
         "prev_rel_3_imp_p_1" = "prev_rel_3_valence_c_1",
         "prev_rel_3_imp_c_1" = "prev_rel_3_imp_p_1") %>% 
  select(c("prolificid","prev_rel_1_imp_c_1",
           "prev_rel_2_imp_c_1","prev_rel_3_imp_c_1")) %>% 
  rename("event1" = "prev_rel_1_imp_c_1",
         "event2" = "prev_rel_2_imp_c_1",
         "event3" = "prev_rel_3_imp_c_1") %>% 
  pivot_longer(cols = starts_with("event"),
               names_to = "event",
               values_to = "pre_event_memory_importance") %>% 
  mutate(pre_event_memory_importance = scales::rescale(pre_event_memory_importance, 
                                                       to = c(1,7), from = c(0,7)))

post2_events <- post2 %>% 
  select(c("prolificid",
           starts_with("related_importance"))) %>% 
  rename("event1" = "related_importance_1",
         "event2" = "related_importance_2",
         "event3" = "related_importance_3") %>% 
  pivot_longer(cols = starts_with("event"),
               names_to = "event",
               values_to = "post2_event_importance")

post2_events_memory <- post3 %>% 
  select(c("prolificid","t4_related_1_imp2_1","t4_related_2_imp2_1",
           "t4_related_3_imp2_1")) %>% 
  rename("event1" = "t4_related_1_imp2_1",
         "event2" = "t4_related_2_imp2_1",
         "event3" = "t4_related_3_imp2_1") %>% 
  pivot_longer(cols = starts_with("event"),
               names_to = "event",
               values_to = "post2_event_memory_importance") %>% 
  mutate(post2_event_memory_importance = scales::rescale(post2_event_memory_importance, 
                                                         to = c(1,7), from = c(0,7)))

events_political_importance <- pre_events %>% 
  full_join(pre_events_memory) %>% 
  full_join(post2_events) %>% 
  full_join(post2_events_memory)

#### models ######
## initial prep for running models ####

# absolute value valence
events_pev_model_abs <- events_personal_valence %>% 
  mutate(pre_events_change = abs(pre_event_memory_valence - pre_event_valence),
         post2_events_change = abs(post2_event_memory_valence - post2_event_valence)) %>% 
  select(-contains("valence")) %>% 
  pivot_longer(cols = contains("change"),
               names_to = "timepoint",
               values_to = "valence_change") %>% 
  mutate(pre_comp = if_else(str_detect(timepoint,"pre"),1,0))

events_pov_model_abs <- events_political_valence %>% 
  mutate(pre_events_change = abs(pre_event_memory_valence - pre_event_valence),
         post2_events_change = abs(post2_event_memory_valence - post2_event_valence)) %>% 
  select(-contains("valence")) %>% 
  pivot_longer(cols = contains("change"),
               names_to = "timepoint",
               values_to = "valence_change") %>% 
  mutate(pre_comp = if_else(str_detect(timepoint,"pre"),1,0))

events_pev_model_abs <- events_pev_model_abs %>% 
  mutate(type = 0)
events_pov_model_abs <- events_pov_model_abs %>% 
  mutate(type = 1)

events_valence_model_abs <- bind_rows(events_pev_model_abs,
                                      events_pov_model_abs)

# not absolute value valence
events_pev_model <- events_personal_valence %>% 
  mutate(pre_events_change = pre_event_memory_valence - pre_event_valence,
         post2_events_change = post2_event_memory_valence - post2_event_valence) %>% 
  select(-contains("valence")) %>% 
  pivot_longer(cols = contains("change"),
               names_to = "timepoint",
               values_to = "valence_change") %>% 
  mutate(pre_comp = if_else(str_detect(timepoint,"pre"),1,0))

events_pov_model <- events_political_valence %>% 
  mutate(pre_events_change = pre_event_memory_valence - pre_event_valence,
         post2_events_change = post2_event_memory_valence - post2_event_valence) %>% 
  select(-contains("valence")) %>% 
  pivot_longer(cols = contains("change"),
               names_to = "timepoint",
               values_to = "valence_change") %>% 
  mutate(pre_comp = if_else(str_detect(timepoint,"pre"),1,0))

events_pev_model <- events_pev_model %>% 
  mutate(type = 0)
events_pov_model <- events_pov_model %>% 
  mutate(type = 1)

events_valence_model <- bind_rows(events_pev_model,
                                  events_pov_model)
events_valence_model <- events_valence_model %>% 
  mutate(pospectrum.b = as.factor(if_else(po_spectrum_scale_1 < 50, 
                                          "Liberal", "Conservative")),
         pospectrum.d = as.factor(if_else(po_spectrum_scale_1 < 50, 0, 1)))

# not absolute value importance
events_pei_model <- events_personal_importance %>% 
  mutate(pre_events_change = pre_event_memory_importance - pre_event_importance,
         post2_events_change = post2_event_memory_importance - post2_event_importance) %>% 
  select(-contains("importance")) %>% 
  pivot_longer(cols = contains("change"),
               names_to = "timepoint",
               values_to = "importance_change") %>% 
  mutate(pre_comp = if_else(str_detect(timepoint,"pre"),1,0))

events_poi_model <- events_political_importance %>% 
  mutate(pre_events_change = pre_event_memory_importance - pre_event_importance,
         post2_events_change = post2_event_memory_importance - post2_event_importance) %>% 
  select(-contains("importance")) %>% 
  pivot_longer(cols = contains("change"),
               names_to = "timepoint",
               values_to = "importance_change") %>% 
  mutate(pre_comp = if_else(str_detect(timepoint,"pre"),1,0))

events_pei_model <- events_pei_model %>% 
  mutate(type = 0)
events_poi_model <- events_poi_model %>% 
  mutate(type = 1)

events_importance_model <- bind_rows(events_pei_model,
                                     events_poi_model)
events_importance_model <- events_importance_model %>% 
  mutate(pospectrum.b = as.factor(if_else(po_spectrum_scale_1 < 50, 
                                          "Liberal", "Conservative")),
         pospectrum.d = as.factor(if_else(po_spectrum_scale_1 < 50, 0, 1)))

# absolute value importance 
events_pei_model_abs <- events_personal_importance %>% 
  mutate(pre_events_change = abs(pre_event_memory_importance - pre_event_importance),
         post2_events_change = abs(post2_event_memory_importance - post2_event_importance)) %>% 
  select(-contains("importance")) %>% 
  pivot_longer(cols = contains("change"),
               names_to = "timepoint",
               values_to = "importance_change") %>% 
  mutate(pre_comp = if_else(str_detect(timepoint,"pre"),1,0))


events_poi_model_abs <- events_political_importance %>% 
  mutate(pre_events_change = abs(pre_event_memory_importance - pre_event_importance),
         post2_events_change = abs(post2_event_memory_importance - post2_event_importance)) %>% 
  select(-contains("importance")) %>% 
  pivot_longer(cols = contains("change"),
               names_to = "timepoint",
               values_to = "importance_change") %>% 
  mutate(pre_comp = if_else(str_detect(timepoint,"pre"),1,0))

events_pei_model_abs <- events_pei_model_abs %>% 
  mutate(type = 0)
events_poi_model_abs <- events_poi_model_abs %>% 
  mutate(type = 1)

events_importance_model_abs <- bind_rows(events_pei_model_abs,
                                         events_poi_model_abs)

# semantic similarity
events_sims <- read_csv('data/events.csv')
events_sims_model <- events_sims %>% 
  pivot_longer(cols = c("post1_pre","post2_pre","post3_pre",
                        "post2_post1","post3_post1","post3_post2"),
               names_to = "comparison",
               values_to = "similarity") %>% 
  mutate(pre_comp = if_else(str_detect(comparison,"pre"),1,0))
## running models ######
# not absolute value valence
summary(lmer(valence_change ~ pre_comp*type +
               (1 + pre_comp*type | prolificid),
             data = events_valence_model))

vm1 <- brm(valence_change ~ pre_comp*type +
             (1 + pre_comp*type | prolificid),
           data = events_valence_model,
           chains = 2, iter = 6000)

vm1a <- brm(valence_change ~ pre_comp*pospectrum.d +
              (1 + pre_comp*pospectrum.d | prolificid),
            data = events_valence_model %>% 
              filter(type == 1),
            chains = 2, iter = 6000)

vm1b <- brm(valence_change ~ pre_comp*pospectrum.d +
              (1 + pre_comp*pospectrum.d | prolificid),
            data = events_valence_model %>% 
              filter(type == 0),
            chains = 2, iter = 6000)

# absolute value valence
summary(lmer(valence_change ~ pre_comp*type +
               (1 + pre_comp*type | prolificid),
             data = events_valence_model_abs))

vm2 <- brm(valence_change ~ pre_comp*type +
             (1 + pre_comp*type | prolificid),
           data = events_valence_model_abs, cores = 4,
           chains = 2, iter = 6000)

# not absolute value importance
summary(lmer(importance_change ~ pre_comp*type +
               (1 + pre_comp*type | prolificid),
             data = events_importance_model))

im1 <- brm(importance_change ~ pre_comp*type +
             (1 + pre_comp*type | prolificid),
           data = events_importance_model,
           cores = 4, chains = 2, iter = 6000)

im1a <- brm(importance_change ~ pre_comp*pospectrum.d +
              (1 + pre_comp*pospectrum.d | prolificid),
            data = events_importance_model %>% 
              filter(type == 1),
            chains = 2, iter = 6000)

im1b <- brm(importance_change ~ pre_comp*pospectrum.d +
              (1 + pre_comp*pospectrum.d | prolificid),
            data = events_importance_model %>% 
              filter(type == 0),
            chains = 2, iter = 6000)

# absolute value importance
summary(lmer(importance_change ~ pre_comp*type +
               (1 + pre_comp*type | prolificid),
             data = events_importance_model_abs))

im2 <- brm(importance_change ~ pre_comp*type +
             (1 + pre_comp*type | prolificid),
           data = events_importance_model_abs,
           cores = 4, chains = 2, iter = 6000)

# semantic similarity
summary(lmer(similarity ~ pre_comp + 
               (1 + pre_comp | prolificid),
             data = events_sims_model))
#### figures ######
ggplot(events_valence_model_abs_fig,
       aes(x = type, y = valence_change,
           fill = pre_comp, group = pre_comp)) +
  geom_bar(stat = "summary", fun = "mean",
           position = "dodge", color = "black") +
  geom_errorbar(aes(ymin = mean_vc - sem_vc, ymax = mean_vc + sem_vc),
                position = position_dodge(width = .9), width = .25) +
  theme_classic() +
  scale_x_discrete(labels = c("Personal events","Political events")) +
  scale_fill_manual(values = c("steelblue","olivedrab"), 
                    labels = c("Across election","After election")) +
  labs(x = "Event type",y = "Event valence change\n(Absolute difference)",
       title = "Changes in perceived valence\nof specific events",
       fill = "Timepoint\ncomparison") +
  coord_cartesian(ylim = c(0,50)) +
  geom_segment(x = 1, xend = 2, y = 42, yend = 42, size = .75) +
  geom_segment(x = .75, xend = 1.25, y = 35, yend = 35, size = .75) +
  geom_segment(x = 1.75, xend = 2.25, y = 35, yend = 35, size = .75) +
  annotate("text",label = "*",x = 1.5, y = 42.2, size = 10) +
  annotate("text",label = "*",x = 1, y = 35.2, size = 10) +
  annotate("text",label = "*",x = 2, y = 35.2, size = 10) +
  theme(axis.title = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 18, hjust = .5))
ggsave("figs/paper/event_valence_q1.jpg", scale = .9)

ggplot(events_importance_model_abs_fig,
       aes(x = type, y = importance_change,
           fill = pre_comp, group = pre_comp)) +
  geom_bar(stat = "summary", fun = "mean",
           position = "dodge", color = "black") +
  geom_errorbar(aes(ymin = mean_vc - sem_vc, ymax = mean_vc + sem_vc),
                position = position_dodge(width = .9), width = .25) +
  theme_classic() +
  scale_x_discrete(labels = c("Personal events","Political events")) +
  scale_fill_manual(values = c("steelblue","olivedrab"), 
                    labels = c("Across election","After election")) +
  labs(x = "Event type",y = "Event importance change\n(Absolute difference)",
       title = "Changes in perceived importance\nof specific events",
       fill = "Timepoint\ncomparison") +
  coord_cartesian(ylim = c(0,7)) +
  geom_segment(x = 1, xend = 2, y = 2.5, yend = 2.5, size = .75) +
  annotate("text",label = "*",x = 1.5, y = 2.7, size = 10) +
  theme(axis.title = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 18, hjust = .5))
ggsave("figs/paper/event_importance_q1.jpg", scale = .9)
##### Q2 - Traits ######
#### Clean/prep ####
pre_traits_mem <- post1 %>% 
  select(c("prolificid",
           starts_with("imagine_traits"))) %>% 
  pivot_longer(cols = starts_with("imagine_traits"),
               names_to = "trait",
               names_prefix = ("imagine_traits_"),
               values_to = "pre_imagine")

post1_traits_mem <- post2 %>% 
  select(c("prolificid",
           starts_with("imagine_traits"))) %>% 
  pivot_longer(cols = starts_with("imagine_traits"),
               names_to = "trait",
               names_prefix = ("imagine_traits_"),
               values_to = "post1_imagine")

post2_traits_mem <- post3 %>% 
  select(c("prolificid",
           starts_with("imagine_traits"))) %>% 
  pivot_longer(cols = starts_with("imagine_traits"),
               names_to = "trait",
               names_prefix = ("imagine_traits_"),
               values_to = "post2_imagine")

traits_mem <- pre_traits %>% 
  full_join(post1_traits) %>% 
  full_join(post2_traits) %>% 
  full_join(pre_traits_mem) %>% 
  full_join(post1_traits_mem) %>% 
  full_join(post2_traits_mem) %>% 
  mutate(pre_mem = abs(pre_imagine - pre),
         post1_mem = abs(post1_imagine - post1),
         post2_mem = abs(post2_imagine - post2),
         trait = as.integer(trait),
         positive = if_else(trait < 11,1,0))

# direction of rating change
extremity_traits_mem <- traits_mem %>% 
  mutate(pre_mem_nonabs = pre_imagine - pre,
         post1_mem_nonabs = post1_imagine - post1,
         post2_mem_nonabs = post2_imagine - post2)
#### models ####
traits_mem_model <- traits_mem %>% 
  select(-c("pre","post1","post2",
            contains("imagine"))) %>% 
  pivot_longer(cols = contains("mem"),
               names_to = "comparison",
               values_to = "difference") %>% 
  mutate(pre_comp = if_else(str_detect(comparison,"pre"),1,0))

summary(lmer(difference ~ pre_comp + 
               (1 + pre_comp | prolificid),
             data = traits_mem_model %>% 
               filter(comparison == "pre_mem" | 
                        comparison == "post1_mem")))

# direction of rating change
extremity_traits_mem_model <- extremity_traits_mem %>% 
  select(-c("pre","pre_imagine","post1","post1_imagine",
            "post2","post2_imagine")) %>% 
  pivot_longer(cols = contains("nonabs"),
               names_to = "comparison",
               values_to = "difference") %>% 
  mutate(pre_comp = if_else(str_detect(comparison,"pre"),1,0))

summary(lmer(difference ~ pre_comp*positive + 
               (1 + pre_comp | prolificid),
             data = extremity_traits_mem_model %>% 
               filter(comparison == "pre_mem_nonabs" | 
                        comparison == "post1_mem_nonabs")))

# political ideology
summary(lmer(difference ~ pre_comp*pospectrum + 
               (1 + pre_comp | prolificid),
             data = extremity_traits_mem_model %>% 
               filter(comparison == "pre_mem_nonabs" | 
                        comparison == "post1_mem_nonabs")))

#### figure ######
traits_mem_figure <- traits_mem_model %>% 
  group_by(prolificid,comparison,pre_comp,pospectrum, positive) %>% 
  summarize(difference = mean(difference)) %>% 
  mutate(pre_comp_fig = if_else(pre_comp == 1,0,1))

ggplot(traits_mem_figure %>% 
         filter(comparison == "pre_mem" | 
                  comparison == "post1_mem"),
       aes(x = pre_comp_fig, y = difference, color = pre_comp_fig)) +
  geom_jitter(alpha = .3, width = .2) +
  geom_smooth(method = "lm", color = "black", linewidth = 1.5) +
  theme_classic() +
  scale_color_gradient(low = "steelblue",high = "olivedrab") +
  scale_x_continuous(breaks = c(0,1),
                     labels = c("Memory\nfor T0","Memory\nfor T1")) +
  labs(x = "Comparison Type",title = "Self-view memory for traits",
       y = "Self-view memory\n(Absolute Difference)") +
  theme(axis.title = element_text(size = 14),
        legend.position = "none",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.title = element_text(size = 16, hjust = .5))
ggsave("figs/paper/traits_q2.jpg")
##### Q2 - Values ######
#### Clean/prep #####
pre_values_mem <- post1 %>% 
  select(c("prolificid",
           starts_with("imagine_values"))) %>% 
  pivot_longer(cols = starts_with("imagine_values"),
               names_to = "value",
               names_prefix = ("imagine_values_"),
               values_to = "pre_imagine")

post1_values_mem <- post2 %>% 
  select(c("prolificid",
           starts_with("imagine_values"))) %>% 
  pivot_longer(cols = starts_with("imagine_values"),
               names_to = "value",
               names_prefix = ("imagine_values_"),
               values_to = "post1_imagine")

post2_values_mem <- post3 %>% 
  select(c("prolificid",
           starts_with("imagine_values"))) %>% 
  pivot_longer(cols = starts_with("imagine_values"),
               names_to = "value",
               names_prefix = ("imagine_values_"),
               values_to = "post2_imagine")

values_mem <- pre_values %>% 
  full_join(post1_values) %>% 
  full_join(post2_values) %>% 
  full_join(pre_values_mem) %>% 
  full_join(post1_values_mem) %>% 
  full_join(post2_values_mem) %>% 
  mutate(pre_mem = abs(pre_imagine - pre),
         post1_mem = abs(post1_imagine - post1),
         post2_mem = abs(post2_imagine - post2),
         value = as.integer(value))

# direction of rating change
extremity_values_mem <- values_mem %>% 
  mutate(pre_mem_nonabs = pre_imagine - pre,
         post1_mem_nonabs = post1_imagine - post1,
         post2_mem_nonabs = post2_imagine - post2,
         pospectrum.d = as.factor(if_else(pospectrum < 50,0, 1)))

#### models #####
values_mem_model <- values_mem %>% 
  select(-c("pre","post1","post2",
            contains("imagine"))) %>% 
  pivot_longer(cols = contains("mem"),
               names_to = "comparison",
               values_to = "difference") %>% 
  mutate(pre_comp = if_else(str_detect(comparison,"pre"),1,0))

summary(lmer(difference ~ pre_comp + 
               (1 + pre_comp | prolificid),
             data = values_mem_model %>% 
               filter(comparison == "pre_mem" | 
                        comparison == "post1_mem")))

# direction of rating change
extremity_values_mem_model <- extremity_values_mem %>% 
  select(-c("pre","pre_imagine","post1","post1_imagine",
            "post2","post2_imagine")) %>% 
  pivot_longer(cols = contains("nonabs"),
               names_to = "comparison",
               values_to = "difference") %>% 
  mutate(pre_comp = if_else(str_detect(comparison,"pre"),1,0))

summary(lmer(difference ~ pre_comp*pospectrum + 
               (1 + pre_comp | prolificid),
             data = extremity_values_mem_model %>% 
               filter(comparison == "pre_mem_nonabs" | 
                        comparison == "post1_mem_nonabs")))
#### figure #####
values_mem_figure <- values_mem_model %>% 
  group_by(prolificid,comparison,pre_comp,pospectrum) %>% 
  summarize(difference = mean(difference)) %>% 
  mutate(pre_comp_fig = if_else(pre_comp == 1,0,1))

# basic fig with 1 comp
ggplot(values_mem_figure %>% 
         filter(comparison == "pre_mem" | 
                  comparison == "post1_mem"),
       aes(x = pre_comp_fig, y = difference, color = pre_comp_fig)) +
  geom_jitter(alpha = .3, width = .2) +
  geom_smooth(method = "lm", color = "black", linewidth = 1.5) +
  theme_classic() +
  scale_color_gradient(low = "steelblue",high = "olivedrab") +
  scale_x_continuous(breaks = c(0,1),
                     labels = c("Memory\nfor T0","Memory\nfor T1")) +
  labs(x = "Comparison Type",title = "Self-view memory for values",
       y = "Self-view memory\n(Absolute Difference)") +
  annotate("text", label = "*", x = .5,y = 2, size = 10) +
  theme(axis.title = element_text(size = 14),
        legend.position = "none",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.title = element_text(size = 16, hjust = .5))
ggsave("figs/paper/values_q2.jpg")
##### Q2 - Events #####
#### Clean/prep #####
## personal, valence ####
pre_events <- pre %>% 
  select(c("prolificid","age","po_spectrum_scale_1",
           contains("unrelated_valence"))) %>% 
  rename("event1" = "unrelated_valence_1",
         "event2" = "unrelated_valence_2",
         "event3" = "unrelated_valence_3") %>% 
  pivot_longer(cols = starts_with("event"),
               names_to = "event",
               values_to = "pre_event_valence")

pre_events_memory <- post2 %>% 
  select(c("prolificid","prev_unr_1_valence_p_1",
           "prev_unr_2_valence_p_1","prev_unr_3_valence_p_1")) %>% 
  rename("event1" = "prev_unr_1_valence_p_1",
         "event2" = "prev_unr_2_valence_p_1",
         "event3" = "prev_unr_3_valence_p_1") %>% 
  pivot_longer(cols = starts_with("event"),
               names_to = "event",
               values_to = "pre_event_memory_valence")

post2_events <- post2 %>% 
  select(c("prolificid",
           contains("unrelated_valence"))) %>% 
  rename("event1" = "unrelated_valence_1",
         "event2" = "unrelated_valence_2",
         "event3" = "unrelated_valence_3") %>% 
  pivot_longer(cols = starts_with("event"),
               names_to = "event",
               values_to = "post2_event_valence")

post2_events_memory <- post3 %>% 
  select(c("prolificid","t4_unrelated_1_val1_1","t4_unrelated_2_val1_1",
           "t4_unrelated_3_val1_1")) %>% 
  rename("event1" = "t4_unrelated_1_val1_1",
         "event2" = "t4_unrelated_2_val1_1",
         "event3" = "t4_unrelated_3_val1_1") %>% 
  pivot_longer(cols = starts_with("event"),
               names_to = "event",
               values_to = "post2_event_memory_valence")

events_personal_valence_2 <- pre_events %>% 
  full_join(pre_events_memory) %>% 
  full_join(post2_events) %>% 
  full_join(post2_events_memory)

## political, valence ##########
pre_events <- pre %>% 
  select(c("prolificid","age","po_spectrum_scale_1",
           starts_with("related_valence"))) %>% 
  rename("event1" = "related_valence_1",
         "event2" = "related_valence_2",
         "event3" = "related_valence_3") %>% 
  pivot_longer(cols = starts_with("event"),
               names_to = "event",
               values_to = "pre_event_valence")

pre_events_memory <- post2 %>% 
  rename("prev_rel_3_valence_p_1" = "prev_rel_event_3_1",
         "prev_rel_3_valence_c_1" = "prev_rel_3_valence_p_1",
         "prev_rel_3_imp_p_1" = "prev_rel_3_valence_c_1",
         "prev_rel_3_imp_c_1" = "prev_rel_3_imp_p_1") %>% 
  select(c("prolificid","prev_rel_1_valence_p_1",
           "prev_rel_2_valence_p_1","prev_rel_3_valence_p_1")) %>% 
  rename("event1" = "prev_rel_1_valence_p_1",
         "event2" = "prev_rel_2_valence_p_1",
         "event3" = "prev_rel_3_valence_p_1") %>% 
  pivot_longer(cols = starts_with("event"),
               names_to = "event",
               values_to = "pre_event_memory_valence")

post2_events <- post2 %>% 
  select(c("prolificid",
           starts_with("related_valence"))) %>% 
  rename("event1" = "related_valence_1",
         "event2" = "related_valence_2",
         "event3" = "related_valence_3") %>% 
  pivot_longer(cols = starts_with("event"),
               names_to = "event",
               values_to = "post2_event_valence")

post2_events_memory <- post3 %>% 
  select(c("prolificid","t4_related_1_val1_1","t4_related_2_val1_1",
           "t4_related_3_val1_1")) %>% 
  rename("event1" = "t4_related_1_val1_1",
         "event2" = "t4_related_2_val1_1",
         "event3" = "t4_related_3_val1_1") %>% 
  pivot_longer(cols = starts_with("event"),
               names_to = "event",
               values_to = "post2_event_memory_valence")

events_political_valence_2 <- pre_events %>% 
  full_join(pre_events_memory) %>% 
  full_join(post2_events) %>% 
  full_join(post2_events_memory)

## personal, importance ########
pre_events <- pre %>% 
  select(c("prolificid","age","po_spectrum_scale_1",
           contains("unrelated_importance"))) %>% 
  rename("event1" = "unrelated_importance_1",
         "event2" = "unrelated_importance_2",
         "event3" = "unrelated_importance_3") %>% 
  pivot_longer(cols = starts_with("event"),
               names_to = "event",
               values_to = "pre_event_importance")

pre_events_memory <- post2 %>% 
  select(c("prolificid","prev_unr_1_imp_p_1",
           "prev_unr_2_imp_p_1","prev_unr_3_imp_p_1")) %>% 
  rename("event1" = "prev_unr_1_imp_p_1",
         "event2" = "prev_unr_2_imp_p_1",
         "event3" = "prev_unr_3_imp_p_1") %>% 
  pivot_longer(cols = starts_with("event"),
               names_to = "event",
               values_to = "pre_event_memory_importance") %>% 
  mutate(pre_event_memory_importance = scales::rescale(pre_event_memory_importance, 
                                                       to = c(1,7), from = c(0,7)))

post2_events <- post2 %>% 
  select(c("prolificid",
           contains("unrelated_importance"))) %>% 
  rename("event1" = "unrelated_importance_1",
         "event2" = "unrelated_importance_2",
         "event3" = "unrelated_importance_3") %>% 
  pivot_longer(cols = starts_with("event"),
               names_to = "event",
               values_to = "post2_event_importance")

post2_events_memory <- post3 %>% 
  select(c("prolificid","t4_unrelated_1_imp1_1","t4_unrelated_2_imp1_1",
           "t4_unrelated_3_imp1_1")) %>% 
  rename("event1" = "t4_unrelated_1_imp1_1",
         "event2" = "t4_unrelated_2_imp1_1",
         "event3" = "t4_unrelated_3_imp1_1") %>% 
  pivot_longer(cols = starts_with("event"),
               names_to = "event",
               values_to = "post2_event_memory_importance") %>% 
  mutate(post2_event_memory_importance = scales::rescale(post2_event_memory_importance, 
                                                         to = c(1,7), from = c(0,7)))

events_personal_importance_2 <- pre_events %>% 
  full_join(pre_events_memory) %>% 
  full_join(post2_events) %>% 
  full_join(post2_events_memory)

## political, importance ##########
pre_events <- pre %>% 
  select(c("prolificid","age","po_spectrum_scale_1",
           starts_with("related_importance"))) %>% 
  rename("event1" = "related_importance_1",
         "event2" = "related_importance_2",
         "event3" = "related_importance_3") %>% 
  pivot_longer(cols = starts_with("event"),
               names_to = "event",
               values_to = "pre_event_importance")

pre_events_memory <- post2 %>% 
  rename("prev_rel_3_valence_p_1" = "prev_rel_event_3_1",
         "prev_rel_3_valence_c_1" = "prev_rel_3_valence_p_1",
         "prev_rel_3_imp_p_1" = "prev_rel_3_valence_c_1",
         "prev_rel_3_imp_c_1" = "prev_rel_3_imp_p_1") %>% 
  select(c("prolificid","prev_rel_1_imp_p_1",
           "prev_rel_2_imp_p_1","prev_rel_3_imp_p_1")) %>% 
  rename("event1" = "prev_rel_1_imp_p_1",
         "event2" = "prev_rel_2_imp_p_1",
         "event3" = "prev_rel_3_imp_p_1") %>% 
  pivot_longer(cols = starts_with("event"),
               names_to = "event",
               values_to = "pre_event_memory_importance") %>% 
  mutate(pre_event_memory_importance = scales::rescale(pre_event_memory_importance, 
                                                       to = c(1,7), from = c(0,7)))

post2_events <- post2 %>% 
  select(c("prolificid",
           starts_with("related_importance"))) %>% 
  rename("event1" = "related_importance_1",
         "event2" = "related_importance_2",
         "event3" = "related_importance_3") %>% 
  pivot_longer(cols = starts_with("event"),
               names_to = "event",
               values_to = "post2_event_importance")

post2_events_memory <- post3 %>% 
  select(c("prolificid","t4_related_1_imp2_1","t4_related_2_imp2_1",
           "t4_related_3_imp2_1")) %>% 
  rename("event1" = "t4_related_1_imp2_1",
         "event2" = "t4_related_2_imp2_1",
         "event3" = "t4_related_3_imp2_1") %>% 
  pivot_longer(cols = starts_with("event"),
               names_to = "event",
               values_to = "post2_event_memory_importance")%>% 
  mutate(post2_event_memory_importance = scales::rescale(post2_event_memory_importance, 
                                                         to = c(1,7), from = c(0,7)))

events_political_importance_2 <- pre_events %>% 
  full_join(pre_events_memory) %>% 
  full_join(post2_events) %>% 
  full_join(post2_events_memory)

#### models #####

## initial prep for running models #####
# valence absolute value 
events_pev_model_2_abs <- events_personal_valence_2 %>% 
  mutate(pre_events_change = abs(pre_event_memory_valence - pre_event_valence),
         post2_events_change = abs(post2_event_memory_valence - post2_event_valence)) %>% 
  select(-contains("valence")) %>% 
  pivot_longer(cols = contains("change"),
               names_to = "timepoint",
               values_to = "valence_change") %>% 
  mutate(pre_comp = if_else(str_detect(timepoint,"pre"),1,0))

events_pov_model_2_abs <- events_political_valence_2 %>% 
  mutate(pre_events_change = abs(pre_event_memory_valence - pre_event_valence),
         post2_events_change = abs(post2_event_memory_valence - post2_event_valence)) %>% 
  select(-contains("valence")) %>% 
  pivot_longer(cols = contains("change"),
               names_to = "timepoint",
               values_to = "valence_change") %>% 
  mutate(pre_comp = if_else(str_detect(timepoint,"pre"),1,0))

events_pev_model_2_abs <- events_pev_model_2_abs %>% 
  mutate(type = 0)
events_pov_model_2_abs <- events_pov_model_2_abs %>% 
  mutate(type = 1)

events_valence_model_2_abs <- bind_rows(events_pev_model_2_abs,
                                        events_pov_model_2_abs)

# valence not absolute value
events_pev_model_2 <- events_personal_valence_2 %>% 
  mutate(pre_events_change = pre_event_memory_valence - pre_event_valence,
         post2_events_change = post2_event_memory_valence - post2_event_valence) %>% 
  select(-contains("valence")) %>% 
  pivot_longer(cols = contains("change"),
               names_to = "timepoint",
               values_to = "valence_change") %>% 
  mutate(pre_comp = if_else(str_detect(timepoint,"pre"),1,0))

events_pov_model_2 <- events_political_valence_2 %>% 
  mutate(pre_events_change = pre_event_memory_valence - pre_event_valence,
         post2_events_change = post2_event_memory_valence - post2_event_valence) %>% 
  select(-contains("valence")) %>% 
  pivot_longer(cols = contains("change"),
               names_to = "timepoint",
               values_to = "valence_change") %>% 
  mutate(pre_comp = if_else(str_detect(timepoint,"pre"),1,0))

events_pev_model_2 <- events_pev_model_2 %>% 
  mutate(type = 0)
events_pov_model_2 <- events_pov_model_2 %>% 
  mutate(type = 1)

events_valence_model_2 <- bind_rows(events_pev_model_2,
                                    events_pov_model_2)

# importance absolute value
events_pei_model_2_abs <- events_personal_importance_2 %>% 
  mutate(pre_events_change = abs(pre_event_memory_importance - pre_event_importance),
         post2_events_change = abs(post2_event_memory_importance - post2_event_importance)) %>% 
  select(-contains("importance")) %>% 
  pivot_longer(cols = contains("change"),
               names_to = "timepoint",
               values_to = "importance_change") %>% 
  mutate(pre_comp = if_else(str_detect(timepoint,"pre"),1,0))

events_poi_model_2_abs <- events_political_importance_2 %>% 
  mutate(pre_events_change = abs(pre_event_memory_importance - pre_event_importance),
         post2_events_change = abs(post2_event_memory_importance - post2_event_importance)) %>% 
  select(-contains("importance")) %>% 
  pivot_longer(cols = contains("change"),
               names_to = "timepoint",
               values_to = "importance_change") %>% 
  mutate(pre_comp = if_else(str_detect(timepoint,"pre"),1,0))

events_pei_model_2_abs <- events_pei_model_2_abs %>% 
  mutate(type = 0)
events_poi_model_2_abs <- events_poi_model_2_abs %>% 
  mutate(type = 1)

events_importance_model_2_abs <- bind_rows(events_pei_model_2_abs,
                                           events_poi_model_2_abs)

# importance not absolute value 
events_pei_model_2 <- events_personal_importance_2 %>% 
  mutate(pre_events_change = pre_event_memory_importance - pre_event_importance,
         post2_events_change = post2_event_memory_importance - post2_event_importance) %>% 
  select(-contains("importance")) %>% 
  pivot_longer(cols = contains("change"),
               names_to = "timepoint",
               values_to = "importance_change") %>% 
  mutate(pre_comp = if_else(str_detect(timepoint,"pre"),1,0))

events_poi_model_2 <- events_political_importance_2 %>% 
  mutate(pre_events_change = pre_event_memory_importance - pre_event_importance,
         post2_events_change = post2_event_memory_importance - post2_event_importance) %>% 
  select(-contains("importance")) %>% 
  pivot_longer(cols = contains("change"),
               names_to = "timepoint",
               values_to = "importance_change") %>% 
  mutate(pre_comp = if_else(str_detect(timepoint,"pre"),1,0))

events_pei_model_2 <- events_pei_model_2 %>% 
  mutate(type = 0)
events_poi_model_2 <- events_poi_model_2 %>% 
  mutate(type = 1)

events_importance_model_2 <- bind_rows(events_pei_model_2,
                                       events_poi_model_2)

## running models ######
# valence not absolute value
summary(lmer(valence_change ~ pre_comp*type +
               (1 + pre_comp*type | prolificid),
             data = events_valence_model_2))

vm3 <- brm(valence_change ~ pre_comp*type +
             (1 + pre_comp*type | prolificid),
           data = events_valence_model_2,
           chains = 2, iter = 6000)

# valence absolute value
summary(lmer(valence_change ~ pre_comp*type +
               (1 + pre_comp*type | prolificid),
             data = events_valence_model_2_abs))

vm4 <- brm(valence_change ~ pre_comp*type +
             (1 + pre_comp*type | prolificid),
           data = events_valence_model_2_abs,
           chains = 2, iter = 6000)

# importance not absolute value
summary(lmer(importance_change ~ pre_comp*type +
               (1 + pre_comp*type | prolificid),
             data = events_importance_model_2))

im3 <- brm(importance_change ~ pre_comp*type +
             (1 + pre_comp*type | prolificid),
           data = events_importance_model_2,
           chains = 2, iter = 6000)

# importance absolute value
summary(lmer(importance_change ~ pre_comp*type +
               (1 + pre_comp*type | prolificid),
             data = events_importance_model_2_abs))

im4 <- brm(importance_change ~ pre_comp*type +
             (1 + pre_comp*type | prolificid),
           data = events_importance_model_2_abs,
           chains = 2, iter = 6000)

#### figures #######
events_valence_model_2_abs_fig <- events_valence_model_2_abs %>% 
  mutate(pre_comp = if_else(pre_comp == 0, 1,0),
         pre_comp = factor(pre_comp,levels = c(0,1)),
         type = factor(type)) %>% 
  group_by(pre_comp,type) %>% 
  mutate(mean_vc = mean(valence_change, na.rm = T),
         sem_vc = sd(valence_change, na.rm = T)/ sqrt(n()))

ggplot(events_valence_model_2_abs_fig,
       aes(x = type, y = valence_change,
           fill = pre_comp, group = pre_comp)) +
  geom_bar(stat = "summary", fun = "mean",
           position = "dodge", color = "black") +
  geom_errorbar(aes(ymin = mean_vc - sem_vc, ymax = mean_vc + sem_vc),
                position = position_dodge(width = .9), width = .25) +
  theme_classic() +
  scale_x_discrete(labels = c("Personal events","Political events")) +
  scale_fill_manual(values = c("steelblue","olivedrab"), 
                    labels = c("Across election","After election")) +
  labs(x = "Event type",y = "Event importance memory\n(Absolute difference)",
       title = "Memory for perceived valence\nof specific events",
       fill = "Timepoint\ncomparison") +
  coord_cartesian(ylim = c(0,50)) +
  geom_segment(x = 1, xend = 2, y = 42, yend = 42, size = .75) +
  geom_segment(x = .75, xend = 1.25, y = 35, yend = 35, size = .75) +
  geom_segment(x = 1.75, xend = 2.25, y = 35, yend = 35, size = .75) +
  annotate("text",label = "*",x = 1.5, y = 42.2, size = 10) +
  annotate("text",label = "*",x = 1, y = 35.2, size = 10) +
  annotate("text",label = "*",x = 2, y = 35.2, size = 10) +
  theme(axis.title = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 18, hjust = .5))
ggsave("figs/paper/event_valence_q2.jpg", scale = .9)

events_importance_model_2_abs_fig <- events_importance_model_2_abs %>% 
  mutate(pre_comp = if_else(pre_comp == 0, 1,0),
         pre_comp = factor(pre_comp,levels = c(0,1)),
         type = factor(type)) %>% 
  group_by(pre_comp,type) %>% 
  mutate(mean_vc = mean(importance_change, na.rm = T),
         sem_vc = sd(importance_change, na.rm = T)/ sqrt(n()))

ggplot(events_importance_model_2_abs_fig,
       aes(x = type, y = importance_change,
           fill = pre_comp, group = pre_comp)) +
  geom_bar(stat = "summary", fun = "mean",
           position = "dodge", color = "black") +
  geom_errorbar(aes(ymin = mean_vc - sem_vc, ymax = mean_vc + sem_vc),
                position = position_dodge(width = .9), width = .25) +
  theme_classic() +
  scale_x_discrete(labels = c("Personal events","Political events")) +
  scale_fill_manual(values = c("steelblue","olivedrab"), 
                    labels = c("Across election","After election")) +
  labs(x = "Event type",y = "Event importance memory\n(Absolute difference)",
       title = "Memory for perceived importance\nof specific events",
       fill = "Timepoint\ncomparison") +
  coord_cartesian(ylim = c(0,7)) +
  geom_segment(x = 1, xend = 2, y = 2.5, yend = 2.5, size = .75) +
  annotate("text",label = "*",x = 1.5, y = 2.7, size = 10) +
  theme(axis.title = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 18, hjust = .5))
ggsave("figs/paper/event_importance_q2.jpg", scale = .9)

