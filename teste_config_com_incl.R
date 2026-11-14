## Configuracao de TESTE: criterio normal (caminho "com inclinacao por UC"
## para especies que convergem), em pastas separadas.
## Uso: CONFIG_ANALISE=teste_config_com_incl.R Rscript roda_especie.R "Cebus unicolor"
source("00_config.R")
versao <- "teste_com_incl"
dir_ajustes    <- paste0("ajustes_", versao)
dir_resultados <- paste0("resultados_", versao)
