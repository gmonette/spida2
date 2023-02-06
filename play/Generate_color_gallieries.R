#'
#' Generate color galleries
#'

spida2::cols()$my25 %>% pal  

raster_pdf('colors_palette.pdf')
for(p in palette.pals()){
  for (i in 2:12) {
    palette.colors(i,palette = p) %>% pal(main = paste0('palette.colors(',i,', ',p,')'))
  }
}
dev.off()
raster_pdf('colors_hcl_qualitative.pdf')
for(p in hcl.pals(type='qualitative')){
  for (i in 2:12) {
    hcl.colors(i,palette = p) %>% pal(main = paste0('hcl.colors(',i,', ',p,')'))
  }
}
dev.off()
raster_pdf('colors_hcl_sequential.pdf')
for(p in hcl.pals(type='sequential')){
  for (i in 2:12) {
    hcl.colors(i,palette = p) %>% pal(main = paste0('hcl.colors(',i,', ',p,')'))
  }
}
dev.off()
raster_pdf('colors_hcl_diverging.pdf')
for(p in hcl.pals(type='diverging')){
  for (i in 2:12) {
    hcl.colors(i,palette = p) %>% pal(main = paste0('hcl.colors(',i,', ',p,')'))
  }
}
dev.off()
