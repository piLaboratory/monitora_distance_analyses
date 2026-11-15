# Notas para o relatório

Registro das decisões de análise e dos achados que precisam entrar no texto.
Atualizar ao longo das análises.

## Preparação dos dados

- **Esforço**: calculado a partir da tabela completa (todas as classes), uma
  linha por transecto e data. Na versão anterior o esforço era tirado só de
  linhas de mamíferos identificados até espécie, o que perdia ~40% do esforço
  e as visitas sem detecção de mamíferos.
- **Filtro taxonômico**: pelo nível validado (`class_taxon_val`), não pelo de campo.
- **Correções de ano/data**: 1 data com ano digitado errado (2014 em vez de
  2024, Parna da Amazônia); 2 casos de mesma UC e data com `ano` diferente
  (Ouro Preto, Tapajós-Arapiuns), corrigidos para o ano mais frequente no dia.
  `ano` é o ano da campanha (pode diferir do ano da data em campanhas que
  atravessam a virada do ano, ex. Viruá).
- **Período**: 2014–2025 (R = 11 anos).
- **Tamanho de grupo**: valores acima do máximo plausível por espécie
  (`data/tamanho_max_grupo.csv`, com referências) marcados como anômalos e
  excluídos da estimativa de tamanho, mas mantidos como 1 avistamento na
  contagem de grupos. Anômalos concentrados em poucas UCs (ex. Hadrosciurus em
  Cazumbá-Iracema; Mazama americana em Tapajós-Arapiuns) → possível diferença de
  registro entre equipes; vale consultar a equipe do Monitora. O padrão se
  repete nas espécies adicionais: anômalos de Dasyprocta iacki concentrados
  em Tapajós-Arapiuns (369 de 605), de D. fuliginosa e Microsciurus em
  Cazumbá-Iracema, de D. leporina em Uatumã. Tapajós-Arapiuns e
  Cazumbá-Iracema aparecem repetidamente.

## Modelo

- **Unidade**: grupo (avistamento), não indivíduo. Com indivíduos havia forte
  sobredispersão (φ até ~32 em Tayassu), explicada pelo agrupamento. Com grupos,
  φ ajustado ≈ 1–2.
- **Agregação**: contagens e esforço por UC-ano.
- **Abundância de grupos** (`ubms::stan_distsamp`, half-normal):
  `~ ano_c ~ ano_c + (ano_c || cnuc)` (g_sem) e o mesmo com `+ (1 | cnuc_ano)`
  (g_com). Principal: g_sem se φ ajustado ≤ 1,5, senão g_com; o outro é
  sensibilidade.
- **Detecção**: `~ ano_c`, sem efeito de UC (efeito de UC na detecção não é
  identificável com poucas detecções por UC: divergências concentradas em
  σ_det alto).
- **Inclinação por UC**: mantida se o g_sem com inclinação convergir
  (Rhat ≤ 1,01 e BFMI ≥ 0,2); senão, modelos só com tendência média.
- **Faixas de distância**: 0–2 m, depois tercis acima de 2 m. ~20% das
  distâncias são 0 (arredondamento); com faixas por quantis a 1ª faixa ficava
  estreita (0–0,4 a 0–1 m) e o modelo half-normal não ajustava.
- **Tamanho de grupo** (`brms`, binomial negativa em n − 1): com distância
  (tamanho avaliado à distância zero, corrige o viés de detecção de grupos
  grandes); contagens parciais censuradas à direita; sensibilidade sem censura
  (só contagens totais).
- **UC classificável**: ≥ 10 grupos em ≥ 3 anos com detecção.

## Achado importante: contagens parciais e tendência de tamanho

- **A tendência de GRUPOS é o resultado principal.** Não depende de
  `n_animais` nem do campo de contagem parcial.
- A fração de contagens marcadas como **parciais aumentou em todas as 15
  espécies** entre 2014–17 e 2022–25, inclusive nas solitárias (Tapirus
  0 → 13%, Eira 5 → 17%, Mazama americana 7 → 18%), e em 15 de 17 combinações
  espécie-UC com dados suficientes nos primatas.
- Nas espécies gregárias (Ateles, Cebus, Sapajus apella, S. macrocephalus,
  Dicotyles, Nasua, Tayassu) a fração subiu 16–32 pontos, enquanto a média das
  contagens totais **caiu** e o mínimo contado nas parciais **não subiu** (caiu
  na maioria).
- Interpretação: mudança na forma de registrar (protocolo/treinamento), não
  aumento real dos grupos. Com censura e distribuição de tamanho de cauda
  pesada, cada parcial "puxa" o tamanho estimado para cima; mais parciais ao
  longo do tempo → tendência de tamanho positiva espúria (v1: Ateles +6 %/ano,
  Cebus +6, S. apella +8; sem censura, Cebus −3).
- Alternativa não descartada: grupos maiores serem mais difíceis de contar
  (aí o modelo sem censura subestimaria). O padrão do mínimo contado pesa
  contra, mas não elimina.
- **No relatório**: indivíduos nas espécies gregárias reportados com as duas
  estimativas (com e sem censura) como faixa de sensibilidade; nas solitárias
  as duas devem coincidir.
- **A fazer**: perguntar à equipe do Monitora se houve mudança de orientação
  sobre marcar contagens como parciais (por volta de 2017–2018).

## Outros pontos para a discussão

- Tendências por UC: com poucas detecções a estimativa é encolhida para a
  média geral; só UCs classificáveis devem ser interpretadas.
- φ do `ubms` (`posterior_predict`, `gof`) é condicional aos dados (abundância
  latente sorteada dado y) → φ e `gof` pouco sensíveis; φ calculado com o
  valor esperado não condicional.
- LOO com efeito por observação (cnuc_ano) tem muitos k de Pareto altos; não
  usado para escolher modelos.

## Resultados da versão 2 (rodada de 14/11/2026)

Tendência de **grupos** (%/ano, IC 95%), modelo principal de cada espécie.
"Incl. UC" = modelo manteve inclinação aleatória por UC (critério de convergência).

| Espécie | Modelo | Incl. UC | Grupos | P(declínio) |
|---|---|---|---|---|
| Tayassu pecari | g_sem | sim | −14,0 (−20,8 a −7,7) | 1,00 |
| Guerlinguetus aestuans | g_sem | sim | −11,3 (−24,6 a −0,3) | 0,98 |
| Eira barbara | g_sem | não | −6,2 (−11,7 a −0,6) | 0,98 |
| Tapirus terrestris | g_sem | não | −4,9 (−11,4 a 2,0) | 0,92 |
| Mazama americana | g_sem | sim | −3,6 (−7,8 a 0,5) | 0,96 |
| Nasua nasua | g_sem | não | −3,5 (−8,2 a 1,2) | 0,93 |
| Hadrosciurus spadiceus | g_com | sim | −3,5 (−9,3 a 2,5) | 0,88 |
| Sapajus macrocephalus | g_com | sim | −3,4 (−8,1 a 1,7) | 0,91 |
| Tamandua tetradactyla | g_sem | não | −2,1 (−7,3 a 3,1) | 0,79 |
| Dicotyles tajacu | g_sem | sim | −2,0 (−6,0 a 2,1) | 0,84 |
| Sapajus apella | g_com | não | −0,9 (−4,2 a 2,4) | 0,72 |
| Sciurillus pusillus | g_com | não | 0,8 (−11,8 a 14,3) | 0,46 |
| Ateles chamek | g_com | sim | 1,8 (−4,7 a 9,1) | 0,31 |
| Mazama nemorivaga | g_sem | não | 2,1 (−2,0 a 6,3) | 0,17 |
| Cebus unicolor | g_com | sim | 2,7 (−3,7 a 9,7) | 0,22 |

- **Declínio com IC fora de zero**: Tayassu, Guerlinguetus, Eira (Eira com
  convergência fraca, ver abaixo). **Declínio provável** (P > 0,9): Tapirus,
  M. americana, Nasua, S. macrocephalus. Nenhuma espécie com aumento claro.
- **Robustez à mudança das faixas de distância**: a tendência de grupos quase
  não mudou entre v1 (faixas por quantis) e v2 (1ª faixa 0–2 m), ex. Tayassu
  −14,2 → −14,0; Guerlinguetus −11,2 → −11,3.
- **Inclinação por UC**: mantida em 8 espécies; 7 ficaram só com tendência
  média (Eira, M. nemorivaga, Nasua, S. apella, Sciurillus, Tamandua, Tapirus).
  S. apella, apesar de muitos dados, não convergiu com inclinação (BFMI 0,08).
- **Indivíduos, com vs. sem censura**: coincidem nas solitárias; divergem nas
  gregárias — Ateles +7,9 vs −0,9; Cebus +9,3 vs −0,2; S. apella +6,9 vs −4,9;
  Dicotyles −0,3 vs −6,1; Nasua −0,3 vs −5,6; S. macrocephalus −3,5 vs −7,5;
  Tayassu −13,0 vs −16,5. Confirma o efeito das contagens parciais (seção
  anterior). Atenção: a queda das contagens totais (que puxa o "sem censura"
  para baixo) também pode refletir protocolo.
- **Convergência dos modelos principais**: todos com Rhat ≤ 1,01 e n_eff > 400
  após o reajuste de 3 espécies com mais iterações (`iter_especie` em
  00_config.R; a decisão sobre a inclinação por UC continuou com as iterações
  padrão):
  - Eira (g_sem, 10.000 iter): Rhat 1,016 → 1,0005; n_eff 409 → 3.959;
    BFMI 0,09 → 0,40. Tendência inalterada: −6,2 (−11,8 a −0,5).
  - Hadrosciurus (g_com, 12.000): Rhat 1,022 → 1,003; tendência −3,5 → −3,4.
  - Cebus (g_com, 12.000): Rhat 1,011 → 1,005; tendência 2,7 → 2,6. Atenção:
    17% das transições atingiram a profundidade máxima da árvore (12), o que
    indica geometria difícil; não invalida os sorteios, mas registrar.
  - As tendências praticamente não mudaram: os problemas de convergência não
    estavam enviesando as estimativas.
  - BFMI entre 0,26 e 0,29 nos g_com de Ateles, Cebus e Hadrosciurus
    (ideal > 0,3): ressalva menor.
- **Técnico**: o R às vezes sofre segfault ao encerrar (liberação de objetos do
  Stan), depois de tudo gravado; não afeta os resultados.
