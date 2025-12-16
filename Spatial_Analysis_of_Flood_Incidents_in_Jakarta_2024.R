# =========================
# STEP R-1: Load Packages
# =========================
library(sf)
library(spatstat.geom)
library(spatstat.explore)
library(spdep)
library(tmap)
library(ggplot2)

# =========================
# STEP R-2: Set Working Directory
# =========================
setwd("D:/TUBES_SPASIAL/RStudio") 

# =========================
# STEP R-3: Load Data
# =========================
titik <- st_read("titik_banjir_utm48s.geojson")
kota  <- st_read("kota_banjir_join.geojson")

# =========================
# STEP R-4: Cek CRS & Geometry
# =========================
st_crs(titik)
st_crs(kota)

st_geometry_type(titik)
st_geometry_type(kota)


# =========================
# STEP R-5A: Mean Center
# =========================
coords <- st_coordinates(titik)

mean_center <- colMeans(coords)

mean_center

mean_center_sf <- st_as_sf(
  data.frame(x = mean_center[1], y = mean_center[2]),
  coords = c("x", "y"),
  crs = st_crs(titik)
)

# =========================
# STEP R-5C: Standard Distance
# =========================
distances <- sqrt(
  (coords[,1] - mean_center[1])^2 +
    (coords[,2] - mean_center[2])^2
)

std_distance <- sqrt(mean(distances^2))

std_distance

# =========================
# STEP R-6A: Create Quadrat Grid
# =========================
grid <- st_make_grid(
  kota,
  cellsize = 5000,   # 5 km x 5 km
  square = TRUE
)

grid_sf <- st_sf(geometry = grid)


# =========================
# STEP R-6B (FIX): Count Points in Grid
# =========================

# Hitung relasi grid - titik
intersections <- st_intersects(grid_sf, titik)

# Hitung jumlah titik di tiap grid
grid_sf$count <- lengths(intersections)

# Cek ringkasan
summary(grid_sf$count)


# =========================
# STEP R-6C: Variance Mean Ratio
# =========================
mean_q <- mean(grid_sf$count)
var_q  <- var(grid_sf$count)

VMR <- var_q / mean_q
VMR


# 7A Library
library(spatstat.geom)
library(spatstat.explore)

# 7B Konversi sf
# =========================
# STEP R-7B: Convert to ppp
# =========================

# Bounding box dari kota
bb <- st_bbox(kota)

# Window spatstat
window <- owin(
  xrange = c(bb["xmin"], bb["xmax"]),
  yrange = c(bb["ymin"], bb["ymax"])
)

pp <- ppp(
  x = st_coordinates(titik)[,1],
  y = st_coordinates(titik)[,2],
  window = window
)

summary(pp)

# 7C Ripley's K&L
# Ripley's K
K <- Kest(pp)
plot(K, main = "Ripley's K Function")

# Ripley's L
L <- Lest(pp)
plot(L, main = "Ripley's L Function")


#Step 8A - Global Moran's I
names(kota)

# 8B library(spdep)

# Neighbors (Queen contiguity)
nb <- poly2nb(kota, queen = TRUE)

# Spatial weights
lw <- nb2listw(nb, style = "W", zero.policy = TRUE)

library(spdep)

# Neighbors (Queen contiguity)
nb <- poly2nb(kota, queen = TRUE)

# Spatial weights
lw <- nb2listw(nb, style = "W", zero.policy = TRUE)


# 8C Hitung Moran's I
moran <- moran.test(kota$NUMPOINTS, lw, zero.policy = TRUE)
moran



# =========================
# STEP R-9: Spatial Regression
# =========================

library(spdep)

# -------------------------
# 9A. Ordinary Least Squares (OLS) — Baseline
# -------------------------
ols_model <- lm(NUMPOINTS ~ bps_Jumlah.Penduduk..Ribu. + bps_Kepadatan.Penduduk.per.km.persegi..Km2., 
                data = kota)

summary(ols_model)

# -------------------------
# 9B. Spatial Autoregressive Model (SAR) — Autokorelasi DV
# -------------------------

# Install 
install.packages("spatialreg")

# Load package
library(spatialreg)

sar_model <- lagsarlm(
  NUMPOINTS ~ bps_Jumlah.Penduduk..Ribu. + bps_Kepadatan.Penduduk.per.km.persegi..Km2.,
  data = kota,
  listw = lw,
  zero.policy = TRUE
)

summary(sar_model)

# -------------------------
# 9C. Spatial Error Model (SEM) — Autokorelasi Error
# -------------------------
sem_model <- errorsarlm(
  NUMPOINTS ~ bps_Jumlah.Penduduk..Ribu. + bps_Kepadatan.Penduduk.per.km.persegi..Km2.,
  data = kota,
  listw = lw,
  zero.policy = TRUE
)

summary(sem_model)

# -------------------------
# 9D. Bandingkan Model (AIC)
# -------------------------
aic_values <- c(
  OLS = AIC(ols_model),
  SAR = AIC(sar_model),
  SEM = AIC(sem_model)
)

print(aic_values)



# GRAFIK/PLOT
# Install 
install.packages("ggforce")

library(sf)
library(ggforce)
library(ggplot2)

# Mean Center & Standard Distance
ggplot() +
  geom_sf(data = kota, fill = NA, color = "black") +                    # Polygon kota
  geom_sf(data = titik, color = "red", size = 2) +                       # Titik banjir
  geom_sf(data = mean_center_sf, color = "green", size = 4) +            # Mean center
  geom_circle(aes(x0 = mean_center[1], y0 = mean_center[2], r = std_distance), 
              color = "blue", linetype = "dashed") +                     # Standard distance
  geom_sf_text(data = kota, aes(label = nm_dati2), size = 2, color = "black") + # Nama kota
  ggtitle("Mean Center & Standard Distance of Jakarta's Flood 2024") +
  theme_minimal()


# Plot Quadrad Analysis
ggplot(grid_sf, aes(x = count)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +
  ggtitle("Quadrat Count Distribution") +
  xlab("Jumlah Banjir per Grid") +
  ylab("Frekuensi")

ggplot() +
  geom_sf(data = grid_sf, aes(fill = count)) +
  geom_sf(data = titik, color = "red", size = 1) +
  scale_fill_viridis_c() +
  ggtitle("Quadrat Count Map") +
  theme_minimal()

# Plot Ripley's K&L
plot(K, main = "Ripley's K Function")
plot(L, main = "Ripley's L Function")

# Plot Moran's L
library(spdep)

moran.plot(
  kota$NUMPOINTS, 
  lw, 
  zero.policy = TRUE, 
  main = "Moran Scatterplot of Flood Incidents in Jakarta 2024",
  xlab = "Flood Count per City",   # label sumbu X
  ylab = "Spatial Lag of Flood Count" # label sumbu Y
)


library(sf)
library(spdep)
library(ggplot2)

# Hitung spatial lag dari NUMPOINTS
kota$lag_NUMPOINTS <- lag.listw(lw, kota$NUMPOINTS, zero.policy = TRUE)

# Buat dataframe untuk ggplot
df_moran <- data.frame(
  NUMPOINTS = kota$NUMPOINTS,
  lag_NUMPOINTS = kota$lag_NUMPOINTS,
  nm_dati2 = kota$nm_dati2
)

# Plot dengan ggplot
ggplot(df_moran, aes(x = NUMPOINTS, y = lag_NUMPOINTS, label = nm_dati2)) +
  geom_point(color = "blue", size = 3) +          # titik
  geom_text(vjust = -0.5, hjust = 0.5, size = 3) + # nama kota
  geom_smooth(method = "lm", se = FALSE, color = "red") + # garis regresi
  labs(
    title = "Moran Scatterplot of Flood Incidents in Jakarta",
    x = "Flood Count per City",
    y = "Spatial Lag of Flood Count"
  ) +
  theme_minimal()


library(sf)
library(spdep)
library(ggplot2)

# Hitung spatial lag dari NUMPOINTS
kota$lag_NUMPOINTS <- lag.listw(lw, kota$NUMPOINTS, zero.policy = TRUE)

# Buat dataframe untuk ggplot
df_moran <- data.frame(
  NUMPOINTS = kota$NUMPOINTS,
  lag_NUMPOINTS = kota$lag_NUMPOINTS,
  nm_dati2 = kota$nm_dati2
)

# Plot dengan warna berbeda per kota
ggplot(df_moran, aes(x = NUMPOINTS, y = lag_NUMPOINTS, color = nm_dati2)) +
  geom_point(size = 4) +                   # Titik kota
  geom_smooth(method = "lm", se = FALSE, color = "black") +  # Garis regresi Moran I
  labs(
    title = "Moran Scatterplot of Flood Incidents in Jakarta",
    x = "Flood Count per City",
    y = "Spatial Lag of Flood Count",
    color = "City"                          # Legend label
  ) +
  theme_minimal() +
  theme(legend.position = "right")



# Plot Spatial Regression
plot(ols_model$residuals, main = "OLS Residuals", ylab = "Residuals", xlab = "Observations")
abline(h = 0, col = "red")

# SAR
library(ggplot2)
library(viridis)

# SAR residuals map
ggplot() +
  geom_sf(data = kota, aes(fill = res_sar), color = "black") +  # Polygon kota
  geom_sf_text(data = kota, aes(label = nm_dati2), size = 2, color = "black") + # Nama kota
  scale_fill_viridis_c(option = "plasma") +
  ggtitle("SAR Model Residuals of Flood Incidents in Jakarta 2024") +
  theme_minimal()


# Pastikan package spatialreg sudah ter-load
library(spatialreg)

# Hitung residual SEM
kota$res_sem <- residuals(sem_model)


# SEM residuals map
ggplot() +
  geom_sf(data = kota, aes(fill = res_sem), color = "black") +
  geom_sf_text(data = kota, aes(label = nm_dati2), size = 2, color = "black") +
  scale_fill_viridis_c(option = "plasma") +
  ggtitle("SEM Model Residuals of Flood Incidents in Jakarta") +
  theme_minimal()




