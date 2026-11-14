## Configuracao de TESTE: forca o caminho "sem inclinacao por UC"
## (conv_max_rhat irrealmente baixo), em pastas separadas.
## Uso: CONFIG_ANALISE=teste_config_sem_incl.R Rscript roda_especie.R "Cebus unicolor"
source("00_config.R")
versao <- "teste_sem_incl"
dir_ajustes    <- paste0("ajustes_", versao)
dir_resultados <- paste0("resultados_", versao)
conv_max_rhat  <- 1.00001   ## nenhum ajuste passa: testa o reajuste sem inclinacao
