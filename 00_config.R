## Configuracoes da analise (lidas por todos os scripts)
## Mudar a versao cria pastas novas de ajustes e resultados; as anteriores
## ficam intactas para comparacao.

versao <- "v2"
dir_ajustes   <- paste0("ajustes_", versao)     ## modelos ajustados (.rds) e escolhas
dir_resultados <- paste0("resultados_", versao) ## tabelas de resultados

## Faixas de distancia: a primeira faixa vai de 0 a 'primeira_faixa_m' metros,
## para absorver o acumulo de distancias registradas como 0 (~20% dos registros);
## as demais tres faixas sao tercis das distancias acima desse limite.
primeira_faixa_m <- 2

## Inclinacao aleatoria por UC: cada especie e ajustada primeiro COM
## inclinacao por UC (g_sem); se nao convergir (Rhat maximo acima de
## conv_max_rhat ou BFMI minimo abaixo de conv_min_bfmi), os modelos sao
## reajustados SEM inclinacao por UC: abundancia ~ ano_c + (1 | cnuc).
conv_max_rhat <- 1.01
conv_min_bfmi <- 0.2

## Escolha do modelo de grupos: g_com se o phi ajustado do g_sem passar disto
limiar_phi <- 1.5

## UC com tendencia propria classificada: minimo de grupos e de anos com deteccao
crit_min_grupos <- 10
crit_min_anos   <- 3
