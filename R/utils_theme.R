## store theme element for consistent plotting. 

make_base_theme <- function(base_size = 5, base_family = "sans") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      text = element_text(size = base_size, family = base_family),
      axis.title = element_text(size = base_size, family = base_family),
      axis.text = element_text(size = base_size, family = base_family),
      legend.title = element_text(size = base_size, family = base_family),
      legend.text = element_text(size = base_size, family = base_family),
      strip.text = element_text(size = base_size, family = base_family),
      legend.key.height = unit(3, "mm"),
      legend.key.width  = unit(3, "mm")
    )
}