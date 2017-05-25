#'@rdname setmytheme
#'@title Set My Theme
#'@description set theme for ggplot, to create consistent figures
#'@export
setmytheme <- function(){
  theme_set(theme_bw())
  theme_update(axis.text.x = element_text(size = 20),
               axis.text.y = element_text(size = 20),
               axis.title.x = element_text(size = 25, margin=margin(20,0,0,0)),
               axis.title.y = element_text(size = 25, angle = 90, margin = margin(0,20,0,0)),
               legend.text = element_text(size = 20),
               legend.title = element_text(size = 20),
               title = element_text(size = 30),
               strip.text = element_text(size = 25),
               strip.background = element_rect(fill="white"),
               panel.spacing = unit(2, "lines"),
               legend.position = "bottom")
}
