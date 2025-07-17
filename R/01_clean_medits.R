library(tidyverse)
library(tidylog) # gives info on data wrangling
library(here)
library(janitor)
library(readxl)
library(rfishbase)
library(geosphere) # for mean trawl position
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)

# Read in data
TA <- read_csv(here('MEDBSsurvey/Demersal/TA.csv'))
TB <- read_csv(here('MEDBSsurvey/Demersal/TB.csv'))
# TC <- read_csv('MEDBSsurvey/Demersal/TC.csv') # Not used for now


# Trawl level data - TA ---------------------------------------------------

# Filter only MEDITS surveys
TA <- TA |> 
  filter(name_of_survey == 'MEDITS')

# Helper function to convert lat/lon format
format_position <- function(x) {
  sign_x <- sign(x)
  x <- abs(x)
  degrees <- floor(x / 100)
  minutes <- x - (degrees * 100)
  decimal <- degrees + (minutes / 60)
  return(sign_x * decimal)
}

# Clean and new columns in TA
TA <- TA |> 
  mutate(
    # Convert coordinates to decimal degrees
    across(c(hauling_longitude, hauling_latitude, shooting_longitude, shooting_latitude), format_position),
    
    # Calculate midpoint using geosphere::midPoint (lon-lat order!)
    midpoint = midPoint(
      p1 = cbind(shooting_longitude, shooting_latitude),
      p2 = cbind(hauling_longitude, hauling_latitude)
    ),
    
    longitude = midpoint[, 1],
    latitude  = midpoint[, 2],
    depth     = (shooting_depth + hauling_depth) / 2,
    
    bottom_salinity = rowMeans(cbind(
      na_if(bottom_salinity_beginning, -1),
      na_if(bottom_salinity_end, -1)
    ), na.rm = TRUE),
    
    bottom_temperature = rowMeans(cbind(
      na_if(bottom_temperature_beginning, -1),
      na_if(bottom_temperature_end, -1)
    ), na.rm = TRUE),
    
    swept_area = wing_opening / 10000000 * distance  # official
  ) |>
  mutate(
    bottom_temperature = ifelse(bottom_temperature == 0, NA, bottom_temperature)
  ) |>
  select(-codend_closing)

# Check that trawl positions make sense
# Download world coastline as 'sf' object
world <- ne_countries(scale = "medium", returnclass = "sf")

# Convert TA (trawl data) to sf object for plotting
TA_sf <- st_as_sf(TA, coords = c("longitude", "latitude"), crs = 4326)

# Plot
ggplot() +
  geom_sf(data = world, fill = "gray95", color = "gray70") +
  geom_sf(data = TA_sf, color = "red", size = 1, alpha = 0.7) +
  coord_sf(xlim = c(-10, 40), ylim = c(32, 50), expand = FALSE) +
  theme_bw() +
  labs(x = "Longitude", y = "Latitude")

### mm it seems Spain put the absolute values of longitude :( for now is safer to remove!


# Catch data - TB ---------------------------------------------------------


# Process TB (catch info)
TB <- TB |> mutate(
  kg = ptot / 100,
  n_individuals = nbtot
)

# Join TA and TB
combo <- left_join(TA, TB)

# Combine genus + species to match with external species names
combo <- combo |> mutate(
  sp      = paste0(genus, species),
  gsa     = area,
  haul_id = paste(year, country, gsa, haul_number, sep = "_")
)

# Load and clean species name mapping
spp_names <- read_excel(here('extra/sp_names.xlsx'), sheet = 1) |> 
  clean_names()

# Manually fix some names
species_data <- data.frame(
  species = c(
    "Trigloporus lastoviza", "Spicara flexuosa", "Torpedo nobiliana",
    "Hoplostethus mediterraneus mediterraneus", "Myctophidae", "Stomias boa boa",
    "Nettastoma melanurum", "Dasyatis centroura", "Blenniidae", "Facciolella oxyrhyncha",
    "Squalus uyato", "Diplecogaster bimaculata bimaculata", "Trachyscorpia cristulata echinata",
    "Diplodus cervinus cervinus", "Liza ramada", "Liza aurata", "Pteromylaeus bovinus",
    "Scomberesox saurus saurus", "Diplodus sargus sargus", "Oblada melanura",
    "Liza saliens", "Salmo trutta trutta", "Lestidiops jayakari jayakari",
    "Lagocephalus lagocephalus lagocephalus"
  ),
  scientific_name = c(
    "Chelidonichthys lastoviza", "Spicara flexuosum", "Tetronarce nobiliana",
    "Hoplostethus mediterraneus", "Myctophidae species", "Stomias boa",
    "Nettastoma melanura", "Bathytoshia centroura", "Blennius ocellaris", "Facciolella oxyrhyncha",
    "Centrophorus uyato", "Diplecogaster bimaculata", "Trachyscorpia cristulata",
    "Diplodus cervinus", "Chelon ramada", "Chelon aurata", "Aetomylaeus bovinus",
    "Scomberesox saurus", "Diplodus sargus", "Oblada melanura",
    "Chelon saliens", "Salmo trutta", "Lestidiops jayakari",
    "Lagocephalus lagocephalus"
  )
)

spp_names_updated <- spp_names %>%
  left_join(species_data, by = c("scientific_name_valid" = "species")) %>%
  mutate(scientific_name_valid = ifelse(!is.na(scientific_name), scientific_name, scientific_name_valid)) %>%
  select(-scientific_name)

# Add valid scientific name
combo <- combo |> mutate(
  sp_scientific = spp_names_updated$scientific_name_valid[match(combo$sp, spp_names_updated$medits_code)],
  sp_scientific = str_replace(sp_scientific, "^NO\\s+", "")
)

# Filter for valid hauls and fish only (catfau codes A, Aa, Ae, Ao)
combo <- combo |> 
  filter(catfau %in% c('A', 'Ao', 'Ae', 'Aa'), validity == 'V') |>
  mutate(
    kg_km2  = kg / swept_area,
    no_km2  = n_individuals / swept_area
  ) |>
  mutate(across(everything(), ~ifelse(is.nan(.), NA, .)))

# Final selected columns
combo <- combo |> select(
  country, gsa, vessel, year, month, day, haul_number, haul_id,
  latitude, longitude, depth,
  bottom_salinity, bottom_temperature, swept_area,
  sp_scientific, n_individuals, no_km2, kg, kg_km2
)

unique(combo$sp_scientific)

#  pelagic or demersal

# Lookup table to group by genus
lookup_table <- combo %>%
  mutate(genus = word(sp_scientific, 1)) %>%
  group_by(sp_scientific, genus) %>%
  summarise(sum_kg_km2 = sum(kg_km2, na.rm = TRUE)) %>%
  group_by(genus) %>%
  top_n(1, sum_kg_km2) %>%
  select(-sum_kg_km2) %>%
  mutate(genus = paste(genus, 'spp.'))

# Retrieve species habitat info
species <- unique(combo$sp_scientific)
habitat_info <- rfishbase::species(species)[, c("Species", "DemersPelag")]

# Join genus-level species for unknowns
join <- lookup_table |> left_join(habitat_info, by = c("sp_scientific" = "Species"))

habitat_info <- rbind(
  habitat_info,
  data.frame(Species = join$genus, DemersPelag = join$DemersPelag)
) |> transmute(sp_scientific = Species, habitat = DemersPelag)

# Clean habitat labels
habitat_info <- habitat_info %>%
  mutate(
    habitat = case_when(
      habitat %in% c("pelagic-oceanic", "pelagic-neritic", "bathypelagic") ~ "pelagic",
      habitat %in% c("demersal", "reef-associated", "bathydemersal", "benthopelagic") ~ "demersal",
      TRUE ~ NA_character_
    )
  )

# Manual overrides for some species
habitat_manual <- tibble::tibble(
  sp_scientific = c(
    "Blenniidae spp.",
    "Cyclothone spp.",
    "Diplecogaster spp.",
    "Facciolella spp.",
    "Hoplostethus spp.",
    "Liza spp.",
    "Mugilidae spp.",
    "Myctophidae spp.",
    "Nettastoma spp.",
    "Oblada spp.",
    "Pomatoschistus spp.",
    "Pteromylaeus spp.",
    "Salmo spp.",
    "Scomberesox spp.",
    "Stomias spp.",
    "Trachyscorpia spp.",
    "Triglidae spp.",
    "Trigloporus spp.",
    "Myctophidae",
    "Blenniidae",
    "Facciolella oxyrhyncha",
    "Myctophidae species",
    "Chelon aurata",
    "Oblada melanura",
    "Triglidae",
    "Mugilidae"
  ),
  habitat = c(
    "demersal",
    "pelagic",
    "demersal",
    "demersal",
    "demersal",
    "demersal",
    "demersal",
    "pelagic",
    "demersal",
    "demersal",
    "demersal",
    "pelagic",
    "demersal",
    "pelagic",
    "pelagic",
    "demersal",
    "demersal",
    "pelagic",
    "pelagic",
    "demersal",
    "demersal",
    "pelagic",
    "demersal",
    "demersal",
    "demersal",
    "demersal"
  )
)

# Final habitat assignment
habitat_info <- habitat_info %>%
  filter(!sp_scientific %in% habitat_manual$sp_scientific) |> 
  full_join(habitat_manual)

# Join habitat info back to combo
combo <- combo |> 
  left_join(habitat_info) |> 
  mutate(sp_category = habitat) |> 
  select(
    country, gsa, vessel, year, month, day, haul_number, haul_id,
    latitude, longitude, depth,
    bottom_salinity, bottom_temperature, swept_area, sp_category,
    sp_scientific, n_individuals, no_km2, kg, kg_km2
  )

# Remove Spain (due to longitude problems)
combo <- combo |> 
  filter(country != 'ESP')

# Save cleaned dataset
write_csv(combo, here('medits_clean.csv'))
