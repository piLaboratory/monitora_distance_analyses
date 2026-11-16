## Gera as tabelas (csv) usadas no relatorio de metodos e resultados
## (relatorio_<versao>.Rmd). Le os modelos ja ajustados (05 e 06) e os
## resultados combinados (07 e 08); NAO ajusta nenhum modelo.
##
## Uso no terminal (a partir da pasta do projeto, depois do 08_roda_todas.R):
##   Rscript 09_tabelas_relatorio.R > logs/09_tabelas.log 2>&1
## Depois:
##   Rscript -e 'rmarkdown::render("relatorio_v2.Rmd")'
##
## Saidas em relatorio_<versao>/tabelas/ (separador decimal: ponto):
##   T1_modelos_selecionados.csv   modelo de grupos e de tamanho por especie
##   T2_deteccao.csv               sigma, largura efetiva (ESW) e p de deteccao
##   T3_densidade_media.csv        densidade media no periodo (grupos e individuos)
##   T4_tendencia_media.csv        tendencia media (%/ano) por especie
##   T5_tendencia_uc.csv           tendencia por UC (especies com inclinacao por UC)
##   S1_amostragem.csv             dados de amostragem por especie
##   S2_convergencia.csv           diagnosticos de convergencia (grupos e tamanho)
##   S3_parametros_grupos.csv      parametros dos modelos de grupos (escala log)
##   S4_deteccao_parametros.csv    parametros de deteccao (escala log e em metros)
##   S5_qualidade_ajuste.csv       qui-quadrado e valor-p bayesiano; dispersao
##   S6_parametros_tamanho.csv     parametros dos modelos de tamanho de grupo
##   S7_comparacao_g_sem_g_com.csv os dois modelos de grupos lado a lado
##   S8_contagens_parciais.csv     fracao de contagens parciais por ano
##   S9_densidade_por_ano.csv      densidade prevista em cada ano (UC tipica)
##
## Retomavel: os calculos de cada especie ficam em relatorio_<versao>/cache/;
## para refazer uma especie, apague o arquivo dela no cache.
## O mais demorado e o gof() (alguns minutos por especie).

library(ubms); library(brms); library(tidyverse)
source(Sys.getenv("CONFIG_ANALISE", "00_config.R"))

dir_rel   <- paste0("relatorio_", versao)
dir_tab   <- file.path(dir_rel, "tabelas")
dir_cache <- file.path(dir_rel, "cache")
dir.create(dir_tab, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_cache, recursive = TRUE, showWarnings = FALSE)

## Configuracao do amostrador usada no 05 (para a Tabela S2)
cfg_grupos <- list(adapt_delta = 0.99, max_treedepth = 12)
cfg_tamanho <- list(adapt_delta = 0.95, max_treedepth = 10)

escolha <- read_csv(file.path(dir_ajustes, "escolha_modelos.csv"), show_col_types = FALSE)
comparacao <- read_csv(file.path(dir_ajustes, "comparacao_modelos.csv"), show_col_types = FALSE)
especies <- escolha$especie

## Resumo de um vetor de sorteios
q3 <- function(x) c(media = mean(x), ic_inf = unname(quantile(x, 0.025)),
                    ic_sup = unname(quantile(x, 0.975)))

## Convergencia de um ajuste do ubms, por grupo de parametros
conv_ubms <- function(fit) {
    s <- as.data.frame(rstan::summary(fit@stanfit)$summary)
    grupo <- case_when(grepl("^beta_state", rownames(s)) ~ "densidade: efeitos fixos",
                       grepl("^sigma_state", rownames(s)) ~ "densidade: desvios dos efeitos aleatorios",
                       grepl("^beta_det", rownames(s))   ~ "deteccao: efeitos fixos",
                       TRUE ~ NA_character_)
    sim <- fit@stanfit@sim
    s |>
        mutate(grupo = grupo) |>
        filter(!is.na(grupo)) |>
        summarise(n_parametros = n(), Rhat_max = max(Rhat, na.rm = TRUE),
                  n_eff_min = min(n_eff, na.rm = TRUE), .by = grupo) |>
        mutate(cadeias = sim$chains, iter = sim$iter, aquecimento = sim$warmup,
               adapt_delta = cfg_grupos$adapt_delta,
               max_treedepth = cfg_grupos$max_treedepth,
               divergencias = rstan::get_num_divergent(fit@stanfit),
               treedepth_max_atingido = rstan::get_num_max_treedepth(fit@stanfit),
               BFMI_min = min(rstan::get_bfmi(fit@stanfit)))
}

## Convergencia de um ajuste do brms (modelo de tamanho)
conv_brms <- function(m) {
    sm <- summary(m)
    tab <- bind_rows(
        mutate(as.data.frame(sm$fixed), grupo = "tamanho: efeitos fixos"),
        mutate(as.data.frame(sm$random$cnuc), grupo = "tamanho: desvios dos efeitos aleatorios"),
        mutate(as.data.frame(sm$spec_pars), grupo = "tamanho: forma (shape)"))
    sim <- m$fit@sim
    tab |>
        summarise(n_parametros = n(), Rhat_max = max(Rhat), n_eff_min = min(Bulk_ESS),
                  .by = grupo) |>
        mutate(cadeias = sim$chains, iter = sim$iter, aquecimento = sim$warmup,
               adapt_delta = cfg_tamanho$adapt_delta,
               max_treedepth = cfg_tamanho$max_treedepth,
               divergencias = rstan::get_num_divergent(m$fit),
               treedepth_max_atingido = rstan::get_num_max_treedepth(m$fit),
               BFMI_min = min(rstan::get_bfmi(m$fit)))
}

## Parametros de um submodelo do ubms em tabela
param_ubms <- function(fit, submodelo) {
    s <- summary(fit, submodelo)
    tibble(submodelo = submodelo, parametro = rownames(s)) |>
        bind_cols(as_tibble(s, .name_repair = "minimal"))
}

## Probabilidade de detectar um grupo presente na faixa [0, w] (transecto
## linear, funcao half-normal): p = sigma * sqrt(2*pi) * (Phi(w/sigma) - 1/2) / w
p_halfnorm <- function(sigma, w) sigma * sqrt(2 * pi) * (pnorm(w / sigma) - 0.5) / w

################################################################################
## Calculos por especie (com cache)
################################################################################
calcula_especie <- function() {
    ## Objetos criados pelo 04 (no ambiente global): dados, dados_todos,
    ## ucs_especie, dist_trunc, dist_breaks, site_covs, y_grp, ua, y_ua, dados_tam
    esc <- filter(escolha, especie == .env$especie)
    estrutura <- if (esc$inclinacao_uc) "incl" else "noincl"
    principal <- esc$modelo_principal
    fits <- map(c(g_sem = "g_sem", g_com = "g_com"), function(m)
        readRDS(arq_modelo(especie, m, estrutura, iter_final(especie, m))))
    fit <- fits[[principal]]
    pref <- file.path(dir_ajustes, gsub(" ", "_", especie))
    m_tam <- list(individuos    = readRDS(paste0(pref, "_tam_dist.rds")),
                  individuos_sc = readRDS(paste0(pref, "_tam_dist_sc.rds")))
    R_anos <- ano_max - ano_min
    w_km <- max(dist_breaks)

    ## S1. Dados de amostragem
    s1 <- tibble(
        especie = especie,
        ucs_ocorrencia = length(ucs_especie),
        ucs_com_esforco = nlevels(ua$cnuc),
        ucs_com_deteccao = n_distinct(dados$cnuc),
        visitas = nrow(site_covs),
        uc_anos = nrow(ua),
        anos = n_distinct(ua$ano),
        periodo = paste(min(ua$ano), max(ua$ano), sep = "-"),
        registros_antes_truncagem = sum(dados_todos$cnuc %in% ucs_especie),
        distancia_truncagem_m = dist_trunc,
        grupos_detectados = sum(y_grp),
        individuos_detectados = sum(dados$n_animais, na.rm = TRUE),
        grupos_sem_tamanho = sum(is.na(dados$n_animais) | dados$n_animais == 0),
        grupos_tamanho_anomalo = sum(dados$anomalo),
        ## Tamanho de grupo ponderado pelos individuos (sum n^2 / sum n): e a
        ## dispersao esperada se a contagem fosse de individuos e os grupos
        ## fossem detectados de forma independente (registros nao anomalos)
        tam_medio_obs = with(filter(dados, !anomalo, n_animais > 0), mean(n_animais)),
        tam_ponderado = with(filter(dados, !anomalo, n_animais > 0),
                             sum(n_animais^2) / sum(n_animais)),
        avistamentos_modelo_tamanho = nrow(dados_tam),
        contagens_parciais = sum(dados_tam$parcial),
        limites_faixas_m = paste(round(dist_breaks * 1000, 1), collapse = " | "),
        esforco_km = sum(ua$tlength))

    ## S2. Convergencia (modelo de grupos principal e modelos de tamanho)
    s2 <- bind_rows(
        mutate(conv_ubms(fit), modelo = principal),
        mutate(conv_brms(m_tam$individuos), modelo = "tamanho com censura"),
        mutate(conv_brms(m_tam$individuos_sc), modelo = "tamanho sem censura")) |>
        mutate(especie = especie, .before = 1) |>
        relocate(modelo, .after = especie)

    ## S3. Parametros dos modelos de grupos (os dois; 'principal' marca o escolhido)
    s3 <- imap_dfr(fits, function(f, m) bind_rows(param_ubms(f, "state"), param_ubms(f, "det")) |>
                       mutate(modelo = m, principal = m == principal)) |>
        mutate(especie = especie, .before = 1) |>
        relocate(modelo, principal, .after = especie)

    ## S4 e T2. Deteccao (modelo principal): log(sigma) ~ ano_c; sigma em km
    bd <- ubms::extract(fit, "beta_det")[[1]]
    bn <- ubms:::beta_names(fit["det"])
    ls_ini <- bd[, match("(Intercept)", bn)]
    ls_fim <- ls_ini + bd[, match("ano_c", bn)]
    sig <- list(ini = exp(ls_ini), fim = exp(ls_fim))          ## km
    det <- imap_dfr(sig, function(s, k) {
        p <- p_halfnorm(s, w_km)
        tibble(ano = if (k == "ini") ano_min else ano_max,
               medida = c("sigma_m", "ESW_m", "p_deteccao"),
               bind_rows(q3(s * 1000), q3(p * w_km * 1000), q3(p)))
    }) |>
        mutate(especie = especie, modelo = principal, truncagem_m = w_km * 1000,
               .before = 1)
    det_var <- q3(100 * (exp(ls_fim - ls_ini) - 1))           ## mudanca de sigma no periodo (%)
    det <- bind_rows(det, tibble(especie = especie, modelo = principal,
                                 truncagem_m = w_km * 1000, ano = NA,
                                 medida = "mudanca_sigma_periodo_pct",
                                 media = det_var[["media"]], ic_inf = det_var[["ic_inf"]],
                                 ic_sup = det_var[["ic_sup"]]))

    ## S5. Qualidade de ajuste (qui-quadrado do ubms, posterior preditiva)
    s5 <- imap_dfr(fits, function(f, m) {
        g <- tryCatch(ubms::gof(f, draws = 500, quiet = TRUE),
                      error = function(e) { message("gof falhou: ", conditionMessage(e)); NULL })
        tibble(especie = especie, modelo = m, principal = m == principal,
               qui_quadrado = if (is.null(g)) NA_real_ else g@estimate,
               valor_p_bayesiano = if (is.null(g)) NA_real_ else g@post_pred_p)
    })

    ## S6. Parametros dos modelos de tamanho
    s6 <- imap_dfr(m_tam, function(m, k) {
        sm <- summary(m)
        bind_rows(as.data.frame(sm$fixed) |> rownames_to_column("parametro"),
                  as.data.frame(sm$random$cnuc) |> rownames_to_column("parametro") |>
                      mutate(parametro = paste0(parametro, " | cnuc")),
                  as.data.frame(sm$spec_pars) |> rownames_to_column("parametro")) |>
            mutate(modelo = if (k == "individuos") "com censura" else "sem censura (so totais)",
                   n_avistamentos = nobs(m), .before = 1)
    }) |>
        mutate(especie = especie, .before = 1)

    ## T3 e S9. Densidade numa UC tipica (efeitos aleatorios = 0), ano a ano.
    ## Tamanho medio do grupo a distancia zero (correcao do vies de tamanho).
    anos <- ano_min:ano_max
    tc <- (anos - ano_min) / R_anos
    bs <- ubms::extract(fit, "beta_state")[[1]]
    bns <- ubms:::beta_names(fit["state"])
    g0 <- bs[, match("(Intercept)", bns)]; g1 <- bs[, match("ano_c", bns)]
    dist_c0 <- (0 - mean(dados_tam$dist.m)) / sd(dados_tam$dist.m)
    eta <- map(m_tam, ~ brms::posterior_linpred(.x, newdata = tibble(ano_c = tc, dist_c = dist_c0),
                                                re_formula = NA))
    set.seed(1)
    nd_s <- min(length(g0), map_int(eta, nrow))
    ig <- sample(length(g0), nd_s)
    it <- map(eta, ~ sample(nrow(.x), nd_s))
    Dg  <- exp(outer(g0[ig], rep(1, length(tc))) + outer(g1[ig], tc))   ## grupos/ha
    tam <- imap(eta, ~ 1 + exp(.x[it[[.y]], , drop = FALSE]))           ## tamanho medio
    por_ano <- list(grupos = Dg, individuos = Dg * tam$individuos,
                    individuos_sc = Dg * tam$individuos_sc,
                    tamanho = tam$individuos, tamanho_sc = tam$individuos_sc)
    ## Densidades em por km2 (o modelo estima por ha: x 100)
    fator <- c(grupos = 100, individuos = 100, individuos_sc = 100, tamanho = 1, tamanho_sc = 1)
    t3 <- imap_dfr(por_ano, ~ c(componente = .y, q3(rowMeans(.x) * fator[[.y]]))) |>
        mutate(across(-componente, as.numeric),
               unidade = if_else(startsWith(componente, "tamanho"), "individuos/grupo", "por km2"),
               especie = especie, modelo = principal, .before = 1)
    s9 <- imap_dfr(por_ano, function(m, k) map_dfr(seq_along(anos), function(j)
        c(componente = k, ano = anos[j], q3(m[, j] * fator[[k]])))) |>
        mutate(across(-componente, as.numeric),
               especie = especie, modelo = principal, .before = 1)

    ## S8. Fracao de contagens parciais por ano (dados do modelo de tamanho)
    s8 <- dados_tam |>
        summarise(n_avistamentos = n(), n_parciais = sum(parcial),
                  frac_parcial = mean(parcial), .by = ano) |>
        arrange(ano) |>
        mutate(especie = especie, .before = 1)

    ## T1 (complemento): coeficiente de distancia no modelo de tamanho
    ## (> 0: grupos maiores detectados mais longe)
    b_dist <- as_draws_df(m_tam$individuos)$b_dist_c
    t1 <- tibble(especie = especie,
                 avistamentos_tamanho = nrow(dados_tam),
                 pct_parciais = 100 * mean(dados_tam$parcial),
                 coef_distancia = mean(b_dist),
                 coef_distancia_ic_inf = unname(quantile(b_dist, 0.025)),
                 coef_distancia_ic_sup = unname(quantile(b_dist, 0.975)))

    list(s1 = s1, s2 = s2, s3 = s3, det = det, s5 = s5, s6 = s6, t3 = t3,
         s9 = s9, s8 = s8, t1 = t1)
}

for (especie in especies) {
    arq_cache <- file.path(dir_cache, paste0(gsub(" ", "_", especie), ".rds"))
    if (file.exists(arq_cache)) next
    message(sprintf("09: %s (%s)", especie, format(Sys.time(), "%H:%M")))
    source("04_prepara_dados_uc_ano.R")
    saveRDS(calcula_especie(), arq_cache)
}

################################################################################
## Junta as especies e grava as tabelas
################################################################################
res <- map(set_names(especies), ~ readRDS(file.path(dir_cache, paste0(gsub(" ", "_", .x), ".rds"))))
junta <- function(k) map_dfr(res, k)
grava <- function(tab, nome) write_csv(tab, file.path(dir_tab, nome))

## T1. Modelos selecionados
t1 <- escolha |>
    left_join(filter(comparacao, modelo == "g_sem") |>
                  select(especie, phi_g_sem = phi, p_eff_g_sem = p_eff,
                         teste_incl_Rhat, teste_incl_BFMI), by = "especie") |>
    mutate(formula_densidade = case_when(
               modelo_principal == "g_sem" &  inclinacao_uc ~ "ano_c + (ano_c || cnuc)",
               modelo_principal == "g_com" &  inclinacao_uc ~ "ano_c + (ano_c || cnuc) + (1 | cnuc_ano)",
               modelo_principal == "g_sem" & !inclinacao_uc ~ "ano_c + (1 | cnuc)",
               modelo_principal == "g_com" & !inclinacao_uc ~ "ano_c + (1 | cnuc) + (1 | cnuc_ano)"),
           formula_deteccao = "log(sigma) ~ ano_c (half-normal)",
           formula_tamanho = "n - 1 | cens ~ ano_c + dist_c + (ano_c || cnuc), binomial negativa",
           iter_final = map2_dbl(especie, modelo_principal, iter_final)) |>
    left_join(junta("t1"), by = "especie") |>
    select(especie, modelo_principal, inclinacao_uc, formula_densidade, formula_deteccao,
           iter_final, teste_incl_Rhat, teste_incl_BFMI, phi_g_sem, p_eff_g_sem,
           phi_aj_g_sem, limiar_phi, formula_tamanho, avistamentos_tamanho,
           pct_parciais, coef_distancia, coef_distancia_ic_inf, coef_distancia_ic_sup)
grava(t1, "T1_modelos_selecionados.csv")

## T2 e S4. Deteccao
det <- junta("det")
t2 <- det |>
    filter(medida != "mudanca_sigma_periodo_pct") |>
    mutate(coluna = paste(medida, ano, sep = "_")) |>
    select(especie, modelo, truncagem_m, coluna, media, ic_inf, ic_sup) |>
    pivot_wider(names_from = coluna, values_from = c(media, ic_inf, ic_sup),
                names_glue = "{coluna}_{.value}") |>
    left_join(filter(det, medida == "mudanca_sigma_periodo_pct") |>
                  select(especie, mudanca_sigma_pct = media,
                         mudanca_sigma_pct_ic_inf = ic_inf, mudanca_sigma_pct_ic_sup = ic_sup),
              by = "especie")
grava(t2, "T2_deteccao.csv")
s4 <- junta("s3") |>
    filter(principal, submodelo == "det") |>
    select(-principal) |>
    bind_rows(det |> transmute(especie, modelo, submodelo = "det (derivado)",
                               parametro = if_else(is.na(ano), medida, paste(medida, ano)),
                               mean = media, `2.5%` = ic_inf, `97.5%` = ic_sup)) |>
    arrange(especie)
grava(s4, "S4_deteccao_parametros.csv")

## T3. Densidade media no periodo
grava(junta("t3"), "T3_densidade_media.csv")
grava(junta("s9"), "S9_densidade_por_ano.csv")

## T4. Tendencia media (resultados do 07/08)
tend <- read_csv(file.path(dir_resultados, "todas_tend_media.csv"), show_col_types = FALSE)
t4 <- tend |>
    mutate(evidencia = case_when(ic_sup < 0 ~ "declinio", ic_inf > 0 ~ "aumento",
                                 TRUE ~ "inconclusivo")) |>
    arrange(especie, desc(principal), componente)
grava(t4, "T4_tendencia_media.csv")

## T5. Tendencia por UC (so especies com inclinacao por UC)
tend_uc <- read_csv(file.path(dir_resultados, "todas_tend_uc.csv"), show_col_types = FALSE)
t5 <- tend_uc |>
    filter(principal, especie %in% escolha$especie[escolha$inclinacao_uc]) |>
    arrange(especie, componente, media)
grava(t5, "T5_tendencia_uc.csv")

## Tabelas suplementares
grava(junta("s1"), "S1_amostragem.csv")
grava(junta("s2"), "S2_convergencia.csv")
grava(filter(junta("s3"), principal) |> select(-principal), "S3_parametros_grupos.csv")
grava(junta("s5") |>
          left_join(select(comparacao, especie, modelo, phi, phi_aj, p_eff),
                    by = c("especie", "modelo")),
      "S5_qualidade_ajuste.csv")
grava(junta("s6"), "S6_parametros_tamanho.csv")
grava(comparacao |>
          left_join(select(junta("s5"), especie, modelo, qui_quadrado, valor_p_bayesiano),
                    by = c("especie", "modelo")) |>
          left_join(select(escolha, especie, modelo_principal), by = "especie") |>
          mutate(principal = modelo == modelo_principal) |>
          select(-modelo_principal),
      "S7_comparacao_g_sem_g_com.csv")
grava(junta("s8"), "S8_contagens_parciais.csv")

message("09: tabelas gravadas em ", dir_tab)
