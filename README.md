# Accident-Prone Road Prediction — Central Kolkata: Model 1 + Model 2 merge

Study area: lat 22.550–22.580, lon 88.340–88.375 (Esplanade / Park Street / Maidan / BBD Bagh / Sealdah side).
All geometry is real OpenStreetMap data (© OpenStreetMap contributors), downloaded via Overpass Turbo.

## Files
- `data/Data_master_500m.csv` — **main table**: one row per 500 m grid cell (42 cells), lat/long of cell centre, road/building fractal features, thin/thick road features, climate summary columns.
- `data/Data_master_250m.csv` — same at 250 m (182 cells). RFD/BFD depend on cell size, so do not compare values across the two files.
- `data/grid_fractal_*.geojson`, `data/roads_classified.geojson` — same data as GIS layers (open in QGIS / geopandas; EPSG:4326).
- `data/climate_master_partial.csv` — ERA5-Land monthly 2016–2025 (temperature, wind speed, dew point, partial CSI).
- `maps/` — thin/thick road extraction map; RFD and BFD maps.
- `code/` — `fractal_pipeline.py` (grid + box counting), `merge_and_maps.py` (climate merge + maps).

## Method
1. Roads filtered to the drivable network (motorway…residential, service; alleys/driveways/private excluded; footways/paths/steps removed).
2. **Thick vs thin road (mentor's road-extraction step):** OSM road class is the proxy for visible road width. Thick = trunk, primary, secondary (+ links). Thin = tertiary, residential, living_street, service, unclassified.
3. Reprojected to UTM 45N (EPSG:32645); whole 500 m (or 250 m) cells centred over the study area, no edge slivers.
4. **RFD / BFD:** each cell rasterised to 128×128; box counting at box sizes 2, 4, 8, 16, 32, 64 px; D = slope of log N(ε) vs log(1/ε). Roads are drawn as lines, buildings as filled footprints. `RFD_R2` / `BFD_R2` = quality of the log-log fit. Function checked on shapes with known answers: filled square = 2.00, straight line = 1.00.
5. `RFD_thick` / `RFD_thin` = the same measure on thick-only / thin-only roads. `RFD_thick` is only meaningful where there are several thick roads (otherwise it is close to 1).
6. `cell_road_class` = Thick-dominant if thick roads are ≥ 40 % of the cell's road length, else Thin-dominant. The 40 % cut is an arbitrary, adjustable choice (`THICK_SHARE_CUT` in `fractal_pipeline.py`).

## Column guide (Data_master)
grid_id, row, col (row 0 = south), lat, long · coverage (share of cell inside study area) · road_length_m, thick_length_m, thin_length_m, thick_share · n_road_segments, n_buildings, built_fraction · RFD, RFD_thick, RFD_thin, BFD · RFD_R2, BFD_R2 · low_data_road (road length < 1000 m per 500 m cell), low_data_built (< 20 buildings) — flagged values are kept but should not be over-interpreted · cell_road_class · clim_* (climate summary).

## Important limitations
- **Climate has no spatial detail.** ERA5-Land is ~9 km resolution, so the whole study area sits in one climate cell. The `clim_*` columns are therefore identical in every row; they are attached so the table is complete, but they cannot explain differences between cells. Climate varies over time (see `climate_master_partial.csv`), not over space.
- **CSI is partial:** temperature, wind speed and dew point only. Rainfall and surface pressure are not yet included. The source file's "Humidity" column was actually dew-point temperature in Kelvin and was renamed.
- Thin/thick uses OSM road class, not measured road width; no accident records are used yet.
- RBFCI, CSRFI and PPS are not computed here because their formulas should be confirmed with the mentor first.
## How to run
python code/fractal_pipeline.py data/raw/roads.geojson data/raw/buildings.geojson out 500
python code/merge_and_maps.py out data/climate_master_partial.csv data/raw/buildings.geojson
