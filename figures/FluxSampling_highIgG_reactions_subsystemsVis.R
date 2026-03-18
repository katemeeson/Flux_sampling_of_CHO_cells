library(tidyverse)
library(patchwork)
phases = c("[4]", "[5]", "[6, 7, 8]", "[11, 12, 14]")
phases_flt = c("[4]", "[6, 7, 8]", "[11, 12, 14]")
phase_names = c('Early exponential','Late exponential','Stationary/death')

# Toggle line 7 and 8 to plot absolute number and proportion :)
#dat <- read_csv("fluxSampling_subsystems_for_ggplot2.csv") %>% 
dat <- read_csv("fluxSampling_subsystems_proportions_for_ggplot2.csv") %>% 
  mutate(down = -down) %>% 
  pivot_longer(c(down, up), names_to = "direction") %>% 
  filter(abs(value) > 1) %>% 
  filter(phase != "[5]") %>%
  pivot_wider(names_from = "phase", values_fill = NA) %>% 
  #interchange 'phase_names'/'phases' to label with name of phase or time points
  pivot_longer(any_of(phase_names), names_to = "phase") %>%
  pivot_wider(names_from = direction, values_fill = NA) %>% 
  arrange(sumValues) %>% 
  mutate(order = 1:nrow(.),
         phase = factor(phase, levels = phase_names, ordered = T))

g_up <- ggplot(dat, aes(x = up, y = reorder(ss, order), color = phase)) +
  geom_segment(aes(x=up, xend=0, y=reorder(ss, order), 
                   yend=reorder(ss, order)),
               alpha= 0.6) +
  geom_point( size=5, alpha = 0.6) +
  #  facet_wrap(~direction, scales = "free_x") +
  theme_minimal() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(),
        panel.grid.minor.x = element_line(),
        strip.text = element_text(vjust = 1),
        axis.title = element_blank(),
        axis.text.y = element_blank()) +
  scale_x_continuous(labels = abs,
                     limits=c(0, max(c(abs(dat$down), abs(dat$up)), na.rm = T))
  ) +
 # scale_y_discrete(position = "right") +
 # geom_vline(xintercept = 0, color = "white", linewidth = 4) +
  scale_color_brewer(palette = "Pastel2", name = "") +
  annotate("segment", 
           xend = max(c(abs(dat$down), abs(dat$up)), na.rm = T)*0.75, 
           x = max(c(abs(dat$down), abs(dat$up)), na.rm = T)*0.25,
           y = 2, yend = 2,
           colour = "darkgrey", alpha=0.6, size= 0.5, 
           arrow=arrow(length=unit(0.2,"cm"), type = "closed")) +
  geom_text(x= max(c(abs(dat$down), abs(dat$up)), na.rm = T)*0.45,
            y = 2,
            label = "Up in high\nIgG solutions",
            size= 3,
            color = "darkgrey")
g_down <- ggplot(dat, aes(x = down, y = reorder(ss, order), color = phase)) +
  geom_segment(aes(x=down, xend=0, y=reorder(ss, order), 
                   yend=reorder(ss, order)),
               alpha = 0.6) +
  geom_point( size=5, alpha = 0.6) +
  #  facet_wrap(~direction, scales = "free_x") +
  theme_minimal() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(),
        panel.grid.minor.x = element_line(),
        strip.text = element_text(vjust = 1),
        axis.title = element_blank(), 
        axis.text.y.right = element_text(hjust = 0.5, margin= margin(r = 0, l = 0))) +
  scale_x_continuous(labels = abs,
                     limits=c(-max(c(abs(dat$down), abs(dat$up)), na.rm = T), 0)
  ) +
  scale_y_discrete(position = "right") +
 # geom_vline(xintercept = 0, color = "white", linewidth = 4) +
  scale_color_brewer(palette = "Pastel2", name = "") +
  annotate("segment", 
           xend = -max(c(abs(dat$down), abs(dat$up)), na.rm = T)*0.75, 
           x = -max(c(abs(dat$down), abs(dat$up)), na.rm = T)*0.25,
           y = 2, yend = 2,
           colour = "darkgrey", alpha=0.6, size= 0.5, 
           arrow=arrow(length=unit(0.2,"cm"), type = "closed")) +
  geom_text(x= -max(c(abs(dat$down), abs(dat$up)), na.rm = T)*0.45,
            y = 2,
            label = "Down in high\nIgG solutions",
            size= 3,
            color = "darkgrey")
g_down + g_up + plot_layout(guides = "collect")

#ggsave("results/figs/fluxSampling_subsystems_nReactions.tiff")


