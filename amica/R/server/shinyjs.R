observe({
  hide(selector = "#navbar li a[data-value=qctab]")
  hide(selector = "#navbar li a[data-value=quanttab]")
  hide(selector = "#navbar li a[data-value=comparemicatab]")
})

shinyjs::onclick("toggleAdvancedHeatmap",
                 shinyjs::toggle(id = "advancedHeatmap"))

shinyjs::onclick("heatmapParams",
                 shinyjs::toggle(id = "toggle_heatmap_params"))

shinyjs::onclick("toggleCompareAmicaInput",
                 shinyjs::toggle(id = "compareAmicaInput"))


### PLOT PARAMS ONCLICK
shinyjs::onclick("pcaParams",
                 shinyjs::toggle(id = "toggle_pca_params"))

shinyjs::onclick("boxplotParams",
                 shinyjs::toggle(id = "toggle_boxplot_params"))

shinyjs::onclick("densityParams",
                 shinyjs::toggle(id = "toggle_density_params"))

shinyjs::onclick("corParams",
                 shinyjs::toggle(id = "toggle_cor_params"))

shinyjs::onclick("cvParams",
                 shinyjs::toggle(id = "toggle_cv_params"))

shinyjs::onclick("contaminantsParams",
                 shinyjs::toggle(id = "toggle_contaminants_params"))

shinyjs::onclick("abundantParams",
                 shinyjs::toggle(id = "toggle_abundant_params"))

shinyjs::onclick("barplotIdParams",
                 shinyjs::toggle(id = "toggle_barplotId_params"))

shinyjs::onclick("barplotMvParams",
                 shinyjs::toggle(id = "toggle_barplotMv_params"))

shinyjs::onclick("overlapHeatmapParams",
                 shinyjs::toggle(id = "toggle_overlapHeatmap_params"))

shinyjs::onclick("scatterParams",
                 shinyjs::toggle(id = "toggle_scatter_params"))

shinyjs::onclick("volcanoParams",
                 shinyjs::toggle(id = "toggle_volcano_params"))

shinyjs::onclick("maParams",
                 shinyjs::toggle(id = "toggle_ma_params"))

shinyjs::onclick("fcParams",
                 shinyjs::toggle(id = "toggle_fc_params"))

shinyjs::onclick("profileParams",
                 shinyjs::toggle(id = "toggle_profile_params"))

shinyjs::onclick("showNodeTable",
                 shinyjs::toggle(id = "toggle_node_table"))

shinyjs::onclick("oraBarParams",
                 shinyjs::toggle(id = "toggle_oraBar_params"))

shinyjs::onclick("fcamicaParams",
                 shinyjs::toggle(id = "toggle_fcamica_params"))

shinyjs::onclick("scatteramicaParams",
                 shinyjs::toggle(id = "toggle_scatteramica_params"))

shinyjs::onclick("multiNameChange",
                 shinyjs::toggle(id = "toggle_multi_name_change"))

shinyjs::onclick("upsetParams",
                 shinyjs::toggle(id = "toggle_upset_params"))

shinyjs::onclick("eulerParams",
                 shinyjs::toggle(id = "toggle_euler_params"))

shinyjs::onclick("corAmicasParams",
                 shinyjs::toggle(id = "toggle_corAmicas_params"))