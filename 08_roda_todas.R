## Roda o fluxo (05 -> 06 -> 07) para todas as especies, uma por vez, cada uma
## num processo R separado (evita conflitos entre ubms e brms e libera memoria).
##
## Retomavel: especies ja concluidas (resultados/<Especie>_concluido.txt) sao
## puladas; dentro de cada especie, ajustes ja salvos em ajustes_grupos/ sao
## reaproveitados. Para refazer uma especie, apague o arquivo _concluido.txt
## (e os .rds em ajustes_grupos/, se quiser reajustar os modelos).
##
## Uso (no R, a partir da pasta do projeto):
##   source("08_roda_todas.R")
## ou no terminal:
##   Rscript 08_roda_todas.R
## Log de cada especie: logs/<Especie>.log

library(tidyverse)

## Criterio de UC classificavel (grupos detectados, anos com deteccao)
crit_min_grupos <- 10
crit_min_anos   <- 3

## Especies: as da tabela de ocorrencia conferida
oc <- read.csv("data/conferencia_ocorrencia_especies_nas_ucs.csv")
especies <- gsub(".", " ", names(oc)[-(1:2)], fixed = TRUE)

dir.create("logs", showWarnings = FALSE)
dir.create("resultados", showWarnings = FALSE)
concluido <- function(sp) file.exists(file.path("resultados",
                                                paste0(gsub(" ", "_", sp), "_concluido.txt")))

status <- tibble(especie = especies, situacao = NA_character_, minutos = NA_real_)
for (i in seq_along(especies)) {
    sp <- especies[i]
    if (concluido(sp)) { status$situacao[i] <- "ja concluida"; next }
    message(sprintf("[%d/%d] %s ... (%s)", i, length(especies), sp, format(Sys.time(), "%H:%M")))
    log <- file.path("logs", paste0(gsub(" ", "_", sp), ".log"))
    t0 <- Sys.time()
    saida <- system2(file.path(R.home("bin"), "Rscript"),
                     c("roda_especie.R", shQuote(sp), crit_min_grupos, crit_min_anos),
                     stdout = log, stderr = log)
    status$minutos[i]  <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    status$situacao[i] <- if (saida == 0 && concluido(sp)) "ok" else "ERRO (ver log)"
    message(sprintf("      %s em %.1f min", status$situacao[i], status$minutos[i]))
}
print(as.data.frame(status), digits = 3)

################################################################################
## Junta os resultados de todas as especies
################################################################################
escolha <- read_csv("ajustes_grupos/escolha_modelos.csv", show_col_types = FALSE) |>
    select(especie, modelo_principal)
junta <- function(sufixo) {
    arqs <- list.files("resultados", pattern = paste0("_", sufixo, "\\.csv$"), full.names = TRUE)
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
write_csv(todas_tend_media, "resultados/todas_tend_media.csv")
write_csv(todas_dens_ind,   "resultados/todas_dens_ind.csv")
write_csv(todas_tend_uc,    "resultados/todas_tend_uc.csv")
file.copy("ajustes_grupos/comparacao_modelos.csv", "resultados/todas_comparacao_modelos.csv",
          overwrite = TRUE)

## Resumo: tendencia de individuos (modelo principal) e UCs por categoria
resumo <-
    todas_tend_media |>
    filter(principal, componente == "individuos") |>
    select(especie, modelo, pct_ano = media, ic_inf, ic_sup, P_declinio) |>
    left_join(todas_tend_uc |>
                  filter(principal) |>
                  count(especie, evidencia) |>
                  pivot_wider(names_from = evidencia, values_from = n, values_fill = 0),
              by = "especie")
print(as.data.frame(resumo), digits = 3)
write_csv(resumo, "resultados/resumo_tendencias.csv")
