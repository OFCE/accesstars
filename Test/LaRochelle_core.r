library(sf)
library(stars)
library(glue)
library(tictoc)
library(devtools)

ville <- "La Rochelle"
localdata <- "~/files/la rochelle/r5"
DVFdata <- "~/files/DVFdata"

# choix d'un petit secteur de test. Core plutôt que commuting
FUA <- st_read("{DVFdata}/sources/FUA/FRA core/FRA_core.shp" |> glue(),
               stringsAsFactors=FALSE) |> st_transform(3035)
fua <- FUA |> dplyr::filter(fuaname == ville) |> dplyr::pull(geometry)


resol <- 200
jour_du_transit <- as.POSIXct("2022-03-09 08:00:00 UTC")

jMem <- "8G"
options(java.parameters = '-Xmx12G')

moteur_r5 <- r5r::setup_r5(localdata, verbose = FALSE, overwrite = FALSE)

# ---- CALCUL DE L'ACCESSIBILITE ----
# Construction des variables d'opportunités
# c200s <- qs::qread(file = glue("{DVFdata}/data_villes/c200_stars.rda"), nthreads = 4)
# larochelle_c200 <- c200s[fua]
# rm(c200s)
# qs::qsave(larochelle_c200, "~/files/la rochelle/la_rochelle_c200.rda", nthreads = 4)
larochelle_c200 <- qs::qread("~/files/la rochelle/la_rochelle_c200.rda", nthreads = 4)
larochelle_emp <- qs::qread("~/files/la rochelle/emp_LaRochelle.rda", nthreads = 4)

# Pour l'instant, test sur core et non commuting
larochelle_emp <- larochelle_emp[fua]
st_dimensions(larochelle_c200) <- st_dimensions(larochelle_emp)

opportunites_stars <- c(larochelle_emp, dplyr::select(larochelle_c200, Ind))

rm(larochelle_c200, larochelle_emp)

tic()
res1 <- st_accessibility(opportunites_stars, "emplois_total",
                         r5r_core = moteur_r5, mode = "TRANSIT",
                         cutoffs = seq(5, 30, 5), verbose = FALSE)
toc() # 14 sec

tic()
res1 <- st_accessibility2(opportunites_stars, scope = 5, "emplois_total",
                         r5r_core = moteur_r5, mode = "TRANSIT",
                         cutoffs = seq(5, 30, 5), verbose = FALSE)
toc() # 101 sec



idINS <- st_coordinates(opportunites_stars) |> coord2idINS(resolution = 200)
opportunites_stars[["idINS"]] <- matrix(idINS, nrow = dim(opportunites_stars)["x"], ncol = dim(opportunites_stars)["y"])

opp_data <- cbind(st_coordinates(opportunites_stars),
                 opportunites_stars[["idINS"]] |> as.vector(),
                 opportunites_stars[["emplois_total"]] |> as.vector())
names(opp_data) <- c("lon", "lat", "id", "emplois")

tic()
res2 <- r5r::accessibility(r5r_core = moteur_r5, origins = opp_data, destinations = opp_data,
                   opportunities_colname = "emplois", mode = "TRANSIT",
                   cutoffs = seq(5, 30, 5), verbose = FALSE)
toc() # 1 sec ?
