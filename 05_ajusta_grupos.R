## Modelo de distancia para a densidade de GRUPOS, dados agregados por UC-ano.
## Ajusta dois modelos e compara:
##   g_sem : abundancia ~ ano_c + (ano_c || cnuc)
##   g_com : abundancia ~ ano_c + (ano_c || cnuc) + (1 | cnuc_ano)
## Criterio de convergencia (00_config.R), decidido por especie no g_sem:
##   1. ajusta g_sem COM inclinacao aleatoria por UC;
##   2. se max Rhat > conv_max_rhat ou BFMI minimo < conv_min_bfmi, a especie
##      passa a usar os dois modelos SEM inclinacao por UC (so tendencia media):
##        g_sem : abundancia ~ ano_c + (1 | cnuc)
##        g_com : abundancia ~ ano_c + (1 | cnuc) + (1 | cnuc_ano)
## Deteccao (ambos): half-normal, log(sigma) ~ ano_c
##   Sem efeito de UC na deteccao: com poucas deteccoes por UC, deteccao e
##   abundancia de cada UC nao sao separaveis (divergencias concentradas em
##   sigma_det alto no teste com Tayassu pecari).
##
## Escolha do modelo principal (por especie):
##   g_sem, se phi ajustado do g_sem <= limiar_phi (00_config.R); senao g_com.
##   phi ajustado = Pearson / (n UC-anos - p_eff), p_eff = p_loo do g_sem.
##   O outro modelo fica como analise de sensibilidade.
##   A escolha e gravada em <dir_ajustes>/escolha_modelos.csv (lida pelo 07).
##
## Uso:
##   especie <- "Tayassu pecari"
##   source("05_ajusta_grupos.R")
## Salva em <dir_ajustes>/<Especie>_<modelo>_<incl|noincl>.rds

library(ubms)
source("04_prepara_dados_uc_ano.R")

## Configuracoes do amostrador
cfg <- list(chains = 6, cores = min(6, parallel::detectCores()),
            adapt_delta = 0.99, max_treedepth = 12)
## Iteracoes por modelo (aquecimento = metade). O g_com precisa de mais
## sorteios para os desvios dos efeitos aleatorios atingirem n_eff > 400.
iter_modelo <- c(g_sem = 4000, g_com = 6000)

dir.create(dir_ajustes, showWarnings = FALSE)
arq <- function(m) file.path(dir_ajustes,
                             paste0(gsub(" ", "_", especie), "_", m, ".rds"))

ajusta <- function(formula, iter) {
    stan_distsamp(formula, data = umf_ua,
                  keyfun = "halfnorm", output = "density", unitsOut = "ha",
                  chains = cfg$chains, iter = iter, warmup = iter / 2,
                  cores = cfg$cores, refresh = 0,
                  control = list(adapt_delta = cfg$adapt_delta,
                                 max_treedepth = cfg$max_treedepth))
}

## Formula dupla do ubms: ~ deteccao ~ abundancia
formulas <- list(
    incl   = list(g_sem = ~ ano_c ~ ano_c + (ano_c || cnuc),
                  g_com = ~ ano_c ~ ano_c + (ano_c || cnuc) + (1 | cnuc_ano)),
    noincl = list(g_sem = ~ ano_c ~ ano_c + (1 | cnuc),
                  g_com = ~ ano_c ~ ano_c + (1 | cnuc) + (1 | cnuc_ano)))

## 1. Convergencia. Criterios: Rhat < 1.01 e n_eff > 400 nos parametros
##    principais (beta_ e sigma_), nenhuma divergencia, BFMI > 0.3.
##    O n_eff minimo dos desvios (sigma_) e mostrado a parte: e o ponto
##    fraco esperado com dados esparsos (funil dos efeitos aleatorios).
convergencia <- function(fit) {
    s <- summary(fit@stanfit)$summary
    principais <- grepl("^beta_|^sigma_", rownames(s))
    sigmas     <- grepl("^sigma_", rownames(s))
    tibble(max_Rhat      = max(s[principais, "Rhat"], na.rm = TRUE),
           min_n_eff     = min(s[principais, "n_eff"], na.rm = TRUE),
           min_n_eff_sigma = if (any(sigmas)) min(s[sigmas, "n_eff"]) else NA_real_,
           divergencias  = rstan::get_num_divergent(fit@stanfit),
           treedepth_max = rstan::get_num_max_treedepth(fit@stanfit),
           min_BFMI      = min(rstan::get_bfmi(fit@stanfit)))
}

convergiu <- function(fit) {
    cv <- convergencia(fit)
    cv$max_Rhat <= conv_max_rhat && cv$min_BFMI >= conv_min_bfmi
}

## Ajusta (ou carrega, se ja existir um ajuste com o mesmo numero de iteracoes)
carrega_ou_ajusta <- function(m, estrutura) {
    nome  <- paste0(m, "_", estrutura)
    salvo <- if (file.exists(arq(nome))) readRDS(arq(nome)) else NULL
    if (!is.null(salvo) && salvo@stanfit@sim$iter == iter_modelo[[m]]) return(salvo)
    t0 <- Sys.time()
    fit <- ajusta(formulas[[estrutura]][[m]], iter_modelo[[m]])
    message(sprintf("%s %s: %.1f min", especie, nome,
                    as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    saveRDS(fit, arq(nome))
    fit
}

## Decisao da estrutura pela convergencia do g_sem com inclinacao por UC
g_sem <- carrega_ou_ajusta("g_sem", "incl")
conv_incl <- convergencia(g_sem)
inclinacao_uc <- convergiu(g_sem)
if (!inclinacao_uc) g_sem <- carrega_ou_ajusta("g_sem", "noincl")
estrutura <- if (inclinacao_uc) "incl" else "noincl"
g_com <- carrega_ou_ajusta("g_com", estrutura)
message(sprintf("%s: g_sem com inclinacao -> Rhat %.3f, BFMI %.2f => %s inclinacao por UC",
                especie, conv_incl$max_Rhat, conv_incl$min_BFMI,
                if (inclinacao_uc) "COM" else "SEM"))

################################################################################
## Comparacao
################################################################################
## 1. Tendencia media de grupos (%/ano), com IC 95%
R_anos <- ano_max - ano_min
tendencia <- function(fit) {
    ## coeficientes fixos de abundancia; a coluna de ano_c e localizada pelo nome
    b <- ubms::extract(fit, "beta_state")[[1]]
    b <- b[, match("ano_c", ubms:::beta_names(fit["state"]))]
    p <- 100 * (exp(b / R_anos) - 1)
    tibble(pct_ano = mean(p), ic_inf = quantile(p, 0.025),
           ic_sup = quantile(p, 0.975), P_declinio = mean(b < 0))
}

## 2. Dispersao de Pearson nos totais por UC-ano, com o valor esperado
##    NAO condicionado aos dados: E[y_i] = lambda_i * area_i * p_i, onde p_i e a
##    probabilidade de detectar um grupo dentro da faixa truncada, calculado em
##    cada sorteio e depois promediado.
##    (posterior_predict/gof do ubms sorteiam a abundancia latente condicionada
##    as contagens observadas, o que deixa o phi artificialmente baixo.)
##    phi ~ 1: variancia compativel com Poisson. No g_com o esperado inclui o
##    efeito de cada UC-ano, entao phi tende a ficar <= 1 por construcao; o phi
##    informativo para sobredispersao e o do g_sem.
dispersao <- function(fit, draws = 1000, p_eff = 1) {
    samp <- ubms:::get_samples(fit, draws)
    lam  <- t(ubms:::sim_lp(fit, "state", transform = TRUE, newdata = NULL,
                            re.form = NULL, samples = samp))           ## UC-ano x sorteio
    lam  <- lam * ubms:::get_area_adjust(fit@response)                 ## densidade -> numero na faixa
    p    <- ubms:::get_p_for_multinom(fit, samp)                       ## UC-ano x (J+1) x sorteio
    J    <- dim(p)[2] - 1
    pdet <- apply(p[, 1:J, , drop = FALSE], c(1, 3), sum)              ## UC-ano x sorteio
    esp  <- rowMeans(lam * pdet)
    obs  <- rowSums(y_ua)
    pearson <- sum((obs - esp)^2 / esp)
    c(phi = pearson / (length(obs) - 1),
      phi_aj = pearson / (length(obs) - p_eff))
}

## Numero efetivo de parametros (p_loo) para ajustar o phi
## ubms::loo explicito: se o brms estiver carregado, loo() seria o do brms
loo_sem <- ubms::loo(g_sem); loo_com <- ubms::loo(g_com)
p_eff <- c(g_sem = loo_sem$estimates["p_loo", "Estimate"],
           g_com = loo_com$estimates["p_loo", "Estimate"])

comparacao <-
    bind_rows(
        bind_cols(modelo = "g_sem", convergencia(g_sem), tendencia(g_sem),
                  p_eff = p_eff[["g_sem"]],
                  as_tibble_row(dispersao(g_sem, p_eff = p_eff[["g_sem"]]))),
        bind_cols(modelo = "g_com", convergencia(g_com), tendencia(g_com),
                  p_eff = p_eff[["g_com"]],
                  as_tibble_row(dispersao(g_com, p_eff = p_eff[["g_com"]])))) |>
    mutate(especie = especie, .before = 1)
print(as.data.frame(comparacao), digits = 3)

################################################################################
## Escolha do modelo principal
################################################################################
phi_aj_sem <- comparacao$phi_aj[comparacao$modelo == "g_sem"]
modelo_principal <- if (phi_aj_sem > limiar_phi) "g_com" else "g_sem"
message(sprintf("%s: phi ajustado do g_sem = %.2f -> modelo principal: %s",
                especie, phi_aj_sem, modelo_principal))

## Grava comparacao e escolha (substitui a linha da especie, se ja existir)
grava <- function(tab, arquivo) {
    if (file.exists(arquivo)) {
        ant <- read_csv(arquivo, show_col_types = FALSE)
        tab <- bind_rows(filter(ant, especie != .env$especie), tab)
    }
    write_csv(arrange(tab, especie), arquivo)
}
grava(mutate(comparacao, inclinacao_uc = inclinacao_uc,
             teste_incl_Rhat = conv_incl$max_Rhat, teste_incl_BFMI = conv_incl$min_BFMI),
      file.path(dir_ajustes, "comparacao_modelos.csv"))
grava(tibble(especie = especie, phi_aj_g_sem = phi_aj_sem,
             limiar_phi = limiar_phi, modelo_principal = modelo_principal,
             inclinacao_uc = inclinacao_uc),
      file.path(dir_ajustes, "escolha_modelos.csv"))

## 3. Ajuste (qui-quadrado posterior preditivo do ubms). Usa a mesma simulacao
##    condicional aos dados; e pouco sensivel a sobredispersao.
##    Em ajustes mal convergidos o sorteio pode falhar (probabilidades NA);
##    nesse caso o gof e pulado e o fluxo continua.
gof_seguro <- function(fit) tryCatch(gof(fit, draws = 500), error = function(e) {
    message("gof falhou: ", conditionMessage(e)); NULL })
gof_sem <- gof_seguro(g_sem); gof_com <- gof_seguro(g_com)
print(gof_sem); print(gof_com)

## 4. LOO. Atencao: com efeito por observacao (cnuc_ano) o PSIS-LOO costuma
##    ter muitos k de Pareto altos; nesse caso a comparacao nao e confiavel e
##    vale mais olhar convergencia, phi e o IC da tendencia.
ms <- loo::loo_compare(list(g_sem = loo_sem, g_com = loo_com))
print(ms)

## Desvio do efeito UC-ano (quanto maior, mais variacao extra-Poisson entre anos)
print(summary(g_com, "state"))
