## Prepara os dados agregados por UC-ano para o modelo de GRUPOS
## e os dados de tamanho de grupo para o modelo censurado.
## Requer, antes de executar:
##   especie : nome valido da especie, ex. "Tayassu pecari"
## Produz:
##   umf_ua    : unmarkedFrameDS, uma linha por UC-ano, contagem de grupos
##   dados_tam : um avistamento por linha, para o modelo de tamanho de grupo

## Contagem de grupos (numero de avistamentos) por visita e classe de distancia
contagem <- "grupos"
source("03_prepara_dados_da_especie.R")

################################################################################
## Agregacao por UC-ano
################################################################################
## Sob Poisson, somar as visitas de um UC-ano nao muda a informacao sobre a
## tendencia; a pseudorreplicacao e tratada no modelo pelo efeito aleatorio
## (1 | cnuc_ano), que aqui e um efeito por observacao.
ua <-
    bind_cols(site_covs, as.data.frame(y_grp)) |>
    group_by(cnuc, ano) |>
    summarise(across(all_of(bin_cols), sum),
              tlength = sum(tlength),       ## km percorridos no UC-ano
              n_visitas = n(),
              .groups = "drop") |>
    mutate(ano_c    = (ano - ano_min) / (ano_max - ano_min),
           cnuc_ano = factor(paste(cnuc, ano, sep = "_")),
           cnuc     = factor(cnuc)) |>
    arrange(cnuc, ano)

y_ua <- as.matrix(ua[, bin_cols]); storage.mode(y_ua) <- "integer"
site_covs_ua <- as.data.frame(select(ua, cnuc, ano, ano_c, cnuc_ano, tlength, n_visitas))

stopifnot(sum(y_ua) == sum(y_grp),
          isTRUE(all.equal(sum(site_covs_ua$tlength), sum(site_covs$tlength))),
          !anyDuplicated(site_covs_ua$cnuc_ano))

umf_ua <- unmarkedFrameDS(
    y           = y_ua,
    siteCovs    = site_covs_ua,
    dist.breaks = dist_breaks,
    survey      = "line",
    tlength     = site_covs_ua$tlength,
    unitsIn     = "km")   ## distancias e comprimentos em km

################################################################################
## Dados de tamanho de grupo (censurado)
################################################################################
## Usa os mesmos avistamentos do modelo de grupos (UCs de ocorrencia, distancias
## truncadas). Exclui n_animais ausente/zero e os anomalos (acima do maximo
## plausivel). Contagens parciais entram como censuradas a direita.
##
## Variavel resposta: s = n_animais - 1 (todo grupo tem ao menos 1 animal).
## Censura no brms: cens = "right" usa log P(S > y). Uma contagem parcial de n
## animais informa que o grupo tinha n ou mais, ou seja S >= n - 1, isto e,
## S > n - 2. Por isso, para as parciais, s_obs = n - 2.
## Parciais com n = 1 nao informam nada (todo grupo tem >= 1) e sao excluidas.
dados_tam <-
    dados |>
    filter(!is.na(n_animais), n_animais > 0, !anomalo,
           !(cont_tot_par %in% "parcial" & n_animais == 1)) |>
    mutate(parcial  = cont_tot_par %in% "parcial",
           ## registros sem informacao de contagem (NA) sao tratados como totais
           censura  = if_else(parcial, "right", "none"),
           s_obs    = if_else(parcial, n_animais - 2, n_animais - 1),
           ano_c    = (ano - ano_min) / (ano_max - ano_min),
           dist_c   = (dist.m - mean(dist.m)) / sd(dist.m),
           cnuc     = factor(cnuc)) |>
    select(cnuc, ano, ano_c, data, dist.m, dist_c, n_animais, parcial,
           censura, s_obs)

message(sprintf(
    "04 (%s): %d UC-anos, %d UCs, %d grupos, %.0f km; tamanho: %d avistamentos (%d parciais).",
    especie, nrow(ua), nlevels(ua$cnuc), sum(y_ua), sum(ua$tlength),
    nrow(dados_tam), sum(dados_tam$parcial)))
