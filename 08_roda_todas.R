## Roda o fluxo (05 -> 06 -> 07) para todas as especies, uma por vez, cada uma
## num processo R separado (evita conflitos entre ubms e brms e libera memoria).
##
## Configuracoes (pastas, faixas, criterios): 00_config.R
## Retomavel: especies ja concluidas (<dir_resultados>/<Especie>_concluido.txt)
## sao puladas; dentro de cada especie, ajustes ja salvos em <dir_ajustes>/ sao
## reaproveitados. Para refazer uma especie, apague o arquivo _concluido.txt
## (e os .rds em <dir_ajustes>/, se quiser reajustar os modelos).
##
## Uso no terminal (a partir da pasta do projeto):
##   mkdir -p logs
##   nohup Rscript 08_roda_todas.R > logs/08_geral.log 2>&1 &
## Log de cada especie: logs/<versao>_<Especie>.log

library(tidyverse)
source(Sys.getenv("CONFIG_ANALISE", "00_config.R"))  ## configuracoes (00_config.R por padrao)

## Especies: as da tabela de ocorrencia conferida
oc <- read.csv("data/conferencia_ocorrencia_especies_nas_ucs.csv")
especies <- gsub(".", " ", names(oc)[-(1:2)], fixed = TRUE)

dir.create("logs", showWarnings = FALSE)
dir.create(dir_resultados, showWarnings = FALSE)
concluido <- function(sp) file.exists(file.path(dir_resultados,
                                                paste0(gsub(" ", "_", sp), "_concluido.txt")))

status <- tibble(especie = especies, situacao = NA_character_, minutos = NA_real_)
for (i in seq_along(especies)) {
    sp <- especies[i]
    if (concluido(sp)) { status$situacao[i] <- "ja concluida"; next }
    message(sprintf("[%d/%d] %s ... (%s)", i, length(especies), sp, format(Sys.time(), "%H:%M")))
    log <- file.path("logs", paste0(versao, "_", gsub(" ", "_", sp), ".log"))
    t0 <- Sys.time()
    saida <- system2(file.path(R.home("bin"), "Rscript"),
                     c("roda_especie.R", shQuote(sp)), stdout = log, stderr = log)
    status$minutos[i]  <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    status$situacao[i] <- if (saida == 0 && concluido(sp)) "ok" else "ERRO (ver log)"
    message(sprintf("      %s em %.1f min", status$situacao[i], status$minutos[i]))
}
print(as.data.frame(status), digits = 3)

################################################################################
## Junta os resultados de todas as especies
################################################################################
escolha <- read_csv(file.path(dir_ajustes, "escolha_modelos.csv"), show_col_types = FALSE) |>
    select(especie, modelo_principal)
junta <- function(sufixo) {
    arqs <- list.files(dir_resultados, pattern = paste0("_", sufixo, "\\.csv$"), full.names = TRUE)
    ## exclui os arquivos combinados de rodadas anteriores (todas_*.csv)
    arqs <- arqs[!startsWith(basename(arqs), "todas_")]
    map_dfr(arqs, function(a) mutate(read_csv(a, show_col_types = FALSE),
                                     across(any_of("cnuc"), as.character))) |>
        left_join(escolha, by = "especie") |>
        mutate(principal = modelo == modelo_principal) |>
        select(-modelo_principal) |>
        arrange(especie, desc(principal))
}
todas_tend_media <- junta("tend_media")
todas_dens_ind   <- junta("dens_ind")
todas_tend_uc    <- junta("tend_uc")
write_csv(todas_tend_media, file.path(dir_resultados, "todas_tend_media.csv"))
write_csv(todas_dens_ind,   file.path(dir_resultados, "todas_dens_ind.csv"))
write_csv(todas_tend_uc,    file.path(dir_resultados, "todas_tend_uc.csv"))
file.copy(file.path(dir_ajustes, "comparacao_modelos.csv"),
          file.path(dir_resultados, "todas_comparacao_modelos.csv"), overwrite = TRUE)

## Resumo (modelo principal): tendencia de GRUPOS (resultado principal), de
## individuos com e sem censura no tamanho, e UCs por categoria (grupos)
fmt <- function(m, a, b) sprintf("%.1f (%.1f a %.1f)", m, a, b)
resumo <-
    todas_tend_media |>
    filter(principal, componente %in% c("grupos", "individuos", "individuos_sc")) |>
    mutate(valor = fmt(media, ic_inf, ic_sup)) |>
    select(especie, modelo, inclinacao_uc, componente, valor) |>
    pivot_wider(names_from = componente, values_from = valor,
                names_prefix = "pct_ano_") |>
    left_join(todas_tend_media |>
                  filter(principal, componente == "grupos") |>
                  select(especie, P_declinio_grupos = P_declinio),
              by = "especie") |>
    left_join(todas_tend_uc |>
                  filter(principal, componente == "grupos") |>
                  count(especie, evidencia) |>
                  pivot_wider(names_from = evidencia, values_from = n, values_fill = 0),
              by = "especie")
print(as.data.frame(resumo))
write_csv(resumo, file.path(dir_resultados, "resumo_tendencias.csv"))
