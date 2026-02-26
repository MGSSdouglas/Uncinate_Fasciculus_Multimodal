############################################################
## Mouse ↔ Human ortholog classification + DE overlap tests with Ximerakis et al. 2019, Nat Neuro, OL cluster only
############################################################

library(readr)
library(dplyr)
library(stringr)

#obtained from MGI website
#https://www.informatics.jax.org/homology.shtml

# MGI homology report (mouse ↔ human orthologs)
mgifile_hom   <- "~/Downloads/HOM_MouseHumanSequence.rpt.txt"

# MGI nomenclature batch report (Input → current Symbol, Input Type)
mgifile_nomen <- "~/Downloads/MGIBatchReport_20251116_160607.txt"

# Ximerakis OL gene list (must have column named "Gene")
ximerakis_file <- "~/Downloads/OL_ximerakis.csv"

# Mouse-only candidates (from manual checks)
mouseonly_file <- "~/Downloads/MouseONLYcanadidates.csv"

#  human DE tables
res_file       <- "~/AgeDE_celltype_subjectBatch_ge5cells_minSubj27_QWTRUE_GlobalAge_3pass_20260214_ALL_CELLTYPES.csv"
res_broad_file <- "~/AgeDE_BROAD_subjectBatch_ge10cells_minSubj27_QWTRUE_GlobalAge_3pass_20260214_ALL_CLUSTERS.csv"


############################################################
## 1. Read MGI homology file and build mouse ↔ human mapping
############################################################

ortho_raw <- read_tsv(mgifile_hom, show_col_types = FALSE)

# Mouse rows
mouse_tbl <- ortho_raw %>%
  filter(`NCBI Taxon ID` == 10090 | `Common Organism Name` == "mouse, laboratory") %>%
  transmute(
    DB_Class_Key   = `DB Class Key`,
    mouse_symbol   = Symbol,
    mouse_Entrez   = `EntrezGene ID`,
    mouse_MGI_ID   = `Mouse MGI ID`,
    mouse_location = `Genetic Location`,
    mouse_coords   = `Genome Coordinates (mouse: GRCm39 human: GRCh38)`
  )

# Human rows
human_tbl <- ortho_raw %>%
  filter(`NCBI Taxon ID` == 9606 | `Common Organism Name` == "human") %>%
  transmute(
    DB_Class_Key   = `DB Class Key`,
    human_symbol   = Symbol,
    human_Entrez   = `EntrezGene ID`,
    human_HGNC_ID  = `HGNC ID`,
    human_OMIM_ID  = `OMIM Gene ID`,
    human_location = `Genetic Location`,
    human_coords   = `Genome Coordinates (mouse: GRCm39 human: GRCh38)`
  )

# Mouse ↔ human pairs (many-to-many is OK; we collapse later)
mouse_human_pairs <- mouse_tbl %>%
  inner_join(human_tbl, by = "DB_Class_Key")

# Collapse so we have one row per mouse symbol, concatenating all human partners
mouse_human <- mouse_human_pairs %>%
  group_by(mouse_symbol) %>%
  summarise(
    mouse_Entrez   = paste(unique(mouse_Entrez), collapse = ";"),
    mouse_MGI_ID   = paste(unique(mouse_MGI_ID), collapse = ";"),
    mouse_location = paste(unique(mouse_location), collapse = ";"),
    mouse_coords   = paste(unique(mouse_coords), collapse = ";"),
    human_symbol   = paste(unique(human_symbol), collapse = ";"),
    human_Entrez   = paste(unique(human_Entrez), collapse = ";"),
    human_HGNC_ID  = paste(unique(human_HGNC_ID), collapse = ";"),
    human_OMIM_ID  = paste(unique(human_OMIM_ID), collapse = ";"),
    human_location = paste(unique(human_location), collapse = ";"),
    human_coords   = paste(unique(human_coords), collapse = ";"),
    .groups = "drop"
  )

############################################################
## 2. MGI nomenclature (Input → current Symbol)
############################################################

nomen <- read_tsv(mgifile_nomen, show_col_types = FALSE)

nomen2 <- nomen %>%
  transmute(
    Gene             = trimws(as.character(Input)),          
    input_type       = trimws(as.character(`Input Type`)),   # e.g. "current symbol", "old symbol"
    mouse_symbol_cur = trimws(as.character(Symbol))          # current official symbol
  ) %>%
  mutate(
    is_current_symbol = str_detect(tolower(input_type), "current"),
    is_old_symbol     = str_detect(tolower(input_type), "old")
  )

############################################################
## 3. Ximerakis OL list + mouse → human orthology classification
############################################################

Ximerakis_OL <- read.csv(ximerakis_file, stringsAsFactors = FALSE)

x_named <- Ximerakis_OL %>%
  mutate(Gene = trimws(as.character(Gene))) %>%
  left_join(nomen2, by = "Gene") %>%
  mutate(
    mouse_symbol_for_ortho = if_else(
      !is.na(mouse_symbol_cur) & mouse_symbol_cur != "",
      mouse_symbol_cur,
      Gene
    ),
    is_regular_symbol = is_current_symbol | is.na(input_type)
  )

x_with_ortho <- x_named %>%
  left_join(mouse_human,
            by = c("mouse_symbol_for_ortho" = "mouse_symbol")) %>%
  mutate(
    has_human_ortholog_curated = !is.na(human_symbol) & human_symbol != "",
    mouse_specific_curated     = !has_human_ortholog_curated
  )

x_regular <- x_with_ortho %>%
  filter(is_regular_symbol)

mouse_spec_regular <- x_regular %>%
  filter(mouse_specific_curated) %>%
  arrange(Gene) %>%
  dplyr::select(Gene, mouse_symbol_for_ortho, input_type,
                has_human_ortholog_curated, mouse_specific_curated)

############################################################
## 4. “Mouse-specific” regular genes: rescue by same human symbol
############################################################

human_symbol_table <- human_tbl %>%
  transmute(
    human_symbol_exact = human_symbol,
    human_symbol_upper = toupper(human_symbol),
    human_Entrez       = human_Entrez,
    human_HGNC_ID      = human_HGNC_ID
  ) %>%
  distinct()

mouse_spec_regular_checked <- mouse_spec_regular %>%
  mutate(mouse_symbol_upper = toupper(mouse_symbol_for_ortho)) %>%
  left_join(human_symbol_table,
            by = c("mouse_symbol_upper" = "human_symbol_upper")) %>%
  mutate(
    has_human_gene_same_symbol = !is.na(human_symbol_exact)
  )

reg_final <- x_regular %>%
  left_join(
    mouse_spec_regular_checked %>%
      dplyr::select(Gene,
                    has_human_gene_same_symbol,
                    human_symbol_exact,
                    human_HGNC_ID),
    by = "Gene"
  ) %>%
  mutate(
    has_human_gene_same_symbol = if_else(
      is.na(has_human_gene_same_symbol), FALSE, has_human_gene_same_symbol
    ),
    classification = case_when(
      has_human_ortholog_curated ~
        "MGI orthology: has human ortholog",
      !has_human_ortholog_curated & has_human_gene_same_symbol ~
        "No MGI ortholog, but human gene with same symbol exists",
      TRUE ~
        "No MGI ortholog and no human gene with same symbol (candidate mouse-specific)"
    )
  )

final_all <- x_with_ortho %>%
  left_join(
    reg_final %>%
      dplyr::select(Gene, classification),
    by = "Gene"
  )

write.csv(final_all,
          "~/Downloads/ximerakis_annot_homolog_classified.csv",
          row.names = FALSE)


############################################################
## 5. Remove candidate mouse-only genes from Ximerakis OL list
############################################################

mouseonly <- read.csv(mouseonly_file, stringsAsFactors = FALSE)
mouseonly <- mouseonly$Gene

Ximerakis_OL_no_mouseonly <- Ximerakis_OL %>%
  filter(!Gene %in% mouseonly)

Ximerakis_OL_no_mouseonly_DEG <- Ximerakis_OL_no_mouseonly %>%
  filter(padj < 0.1)


############################################################
## 6. Human → mouse orthology for DE tables
##    (find candidate human-specific genes)
############################################################

# Human DE 
res <- read.csv(res_file, stringsAsFactors = FALSE) %>%
  filter(cluster_id %in% c("OL1", "OL2", "OL3"))

res_broad <- read.csv(res_broad_file, stringsAsFactors = FALSE) %>%
  filter(cluster_id %in% "OL")

# Human → mouse mapping (collapse by human symbol)
human_mouse <- mouse_human_pairs %>%
  group_by(human_symbol) %>%
  summarise(
    human_Entrez   = paste(unique(human_Entrez), collapse = ";"),
    human_HGNC_ID  = paste(unique(human_HGNC_ID), collapse = ";"),
    human_location = paste(unique(human_location), collapse = ";"),
    human_coords   = paste(unique(human_coords), collapse = ";"),
    mouse_symbol   = paste(unique(mouse_symbol), collapse = ";"),
    mouse_Entrez   = paste(unique(mouse_Entrez), collapse = ";"),
    mouse_MGI_ID   = paste(unique(mouse_MGI_ID), collapse = ";"),
    mouse_location = paste(unique(mouse_location), collapse = ";"),
    mouse_coords   = paste(unique(mouse_coords), collapse = ";"),
    .groups = "drop"
  )

mouse_symbol_table <- mouse_tbl %>%
  transmute(
    mouse_symbol_exact = mouse_symbol,
    mouse_symbol_upper = toupper(mouse_symbol),
    mouse_Entrez       = mouse_Entrez,
    mouse_MGI_ID       = mouse_MGI_ID
  ) %>%
  distinct()

#  annotate a human gene df for human-specificity
annotate_human_specific <- function(df, gene_col = "gene") {
  
  genes_tbl <- df %>%
    mutate(
      gene_symbol = trimws(as.character(.data[[gene_col]]))
    )
  
  human_with_ortho <- genes_tbl %>%
    left_join(human_mouse,
              by = c("gene_symbol" = "human_symbol")) %>%
    mutate(
      has_mouse_ortholog_curated = !is.na(mouse_symbol) & mouse_symbol != ""
    )
  
  human_spec <- human_with_ortho %>%
    filter(!has_mouse_ortholog_curated) %>%
    distinct(gene_symbol)
  
  human_spec_checked <- human_spec %>%
    mutate(human_symbol_upper = toupper(gene_symbol)) %>%
    left_join(mouse_symbol_table,
              by = c("human_symbol_upper" = "mouse_symbol_upper")) %>%
    mutate(
      has_mouse_gene_same_symbol = !is.na(mouse_symbol_exact)
    ) %>%
    dplyr::select(
      gene_symbol,
      has_mouse_gene_same_symbol,
      mouse_symbol_exact,
      mouse_MGI_ID
    )
  
  human_final <- human_with_ortho %>%
    left_join(human_spec_checked, by = "gene_symbol") %>%
    mutate(
      has_mouse_gene_same_symbol = if_else(
        is.na(has_mouse_gene_same_symbol), FALSE, has_mouse_gene_same_symbol
      ),
      classification_human_side = case_when(
        has_mouse_ortholog_curated ~
          "MGI orthology: has mouse ortholog",
        !has_mouse_ortholog_curated & has_mouse_gene_same_symbol ~
          "No MGI ortholog, but mouse gene with same symbol exists",
        TRUE ~
          "No MGI ortholog and no mouse gene with same symbol (candidate human-specific)"
      )
    )
  
  human_final
}

# Annotate both DE tables
res_annot       <- annotate_human_specific(res,        gene_col = "gene")
res_broad_annot <- annotate_human_specific(res_broad,  gene_col = "gene")

# Candidate human-specific in each
human_only_res <- res_annot %>%
  filter(classification_human_side ==
           "No MGI ortholog and no mouse gene with same symbol (candidate human-specific)") %>%
  pull(gene_symbol)

human_only_res_broad <- res_broad_annot %>%
  filter(classification_human_side ==
           "No MGI ortholog and no mouse gene with same symbol (candidate human-specific)") %>%
  pull(gene_symbol)

# Filter DE tables to remove candidate human-only genes
res_filtered <- res %>%
  filter(!gene %in% human_only_res)

res_broad_filtered <- res_broad %>%
  filter(!gene %in% human_only_res_broad)


############################################################
## 7. Hypergeometric test for overlap:
##    Current dataset DE vs Ximerakis DE (mouse-only removed) for OL clusters
############################################################

# OL1,2,3

df_summary <- res_filtered %>%
  arrange(adj.P.Val) %>%
  filter(adj.P.Val < 0.1)

# Uppercase, trimmed
my_de   <- toupper(unique(trimws(df_summary$gene)))
xim_de  <- toupper(unique(trimws(Ximerakis_OL_no_mouseonly_DEG$Gene)))

universe <- union(res_filtered$gene,
                  Ximerakis_OL_no_mouseonly$Gene) %>%
  trimws() %>%
  toupper() %>%
  unique()

# Restrict DE lists to the universe
my_de   <- intersect(my_de, universe)
xim_de  <- intersect(xim_de, universe)

# Counts
a <- length(intersect(my_de, xim_de))
b <- length(setdiff(my_de, xim_de))
c <- length(setdiff(xim_de, my_de))
d <- length(setdiff(universe, union(my_de, xim_de)))

a + b + c + d  # should equal length(universe)

contingency <- matrix(
  c(a, b,
    c, d),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(
    "My_DE"     = c("DE", "Not_DE"),
    "Ximerakis" = c("DE", "Not_DE")
  )
)

fisher_res <- fisher.test(contingency, alternative = "greater")
fisher_res



## Broad OL

df_summary_broad <- res_broad_filtered %>%
  arrange(adj.P.Val) %>%
  filter(adj.P.Val < 0.1)

my_de_broad <- toupper(unique(trimws(df_summary_broad$gene)))
# same Ximerakis DE list:
xim_de_broad <- toupper(unique(trimws(Ximerakis_OL_no_mouseonly_DEG$Gene)))

universe_broad <- union(res_broad_filtered$gene,
                        Ximerakis_OL_no_mouseonly$Gene) %>%
  trimws() %>%
  toupper() %>%
  unique()

my_de_broad  <- intersect(my_de_broad, universe_broad)
xim_de_broad <- intersect(xim_de_broad, universe_broad)

a_broad <- length(intersect(my_de_broad, xim_de_broad))
b_broad <- length(setdiff(my_de_broad, xim_de_broad))
c_broad <- length(setdiff(xim_de_broad, my_de_broad))
d_broad <- length(setdiff(universe_broad, union(my_de_broad, xim_de_broad)))

a_broad + b_broad + c_broad + d_broad  # should equal length(universe_broad)

contingency_broad <- matrix(
  c(a_broad, b_broad,
    c_broad, d_broad),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(
    "My_DE"     = c("DE", "Not_DE"),
    "Ximerakis" = c("DE", "Not_DE")
  )
)

fisher_res_broad <- fisher.test(contingency_broad, alternative = "greater")
fisher_res_broad

####################################################################################################################################

overlapping_genes <- intersect(my_de, xim_de)
overlapping_genes_broad <- intersect(my_de_broad, xim_de_broad)