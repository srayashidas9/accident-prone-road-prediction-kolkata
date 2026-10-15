"""
Real-data pipeline: OSM roads + buildings (GeoJSON from Overpass Turbo) -> grid -> RFD / BFD
usage: python3 -I fractal_pipeline.py ROADS.geojson BUILDINGS.geojson OUTDIR CELL_M
"""
import sys, os
import numpy as np
import pandas as pd
import geopandas as gpd
import shapely
from shapely.geometry import box
from pyproj import Transformer

roads_fp, bld_fp, outdir, cell_m = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4])
os.makedirs(outdir, exist_ok=True)

# study area = the bbox used in the Overpass query (south, west, north, east)
S, W, N_, E = 22.550, 88.340, 22.580, 88.375
CRS_M = 32645  # UTM 45N, metres

roads = gpd.read_file(roads_fp)
bld = gpd.read_file(bld_fp)

# ---------- road filtering: mimic osmnx 'drive' network ----------
drive = {"motorway", "trunk", "primary", "secondary", "tertiary", "unclassified", "residential",
         "living_street", "service", "motorway_link", "trunk_link", "primary_link",
         "secondary_link", "tertiary_link"}
roads = roads[roads.geometry.type == "LineString"]
roads = roads[roads["highway"].isin(drive)]
if "service" in roads.columns:
    roads = roads[~roads["service"].isin(["alley", "driveway", "emergency_access", "parking",
                                          "parking_aisle", "private"])]
if "access" in roads.columns:
    roads = roads[~roads["access"].isin(["private", "no"])]

# thick vs thin: OSM functional class as the proxy for how wide / visible the road is
THICK = {"motorway", "trunk", "primary", "secondary",
         "motorway_link", "trunk_link", "primary_link", "secondary_link"}
roads = roads.assign(road_class=np.where(roads["highway"].isin(THICK), "thick", "thin"))

roads = roads[["highway", "road_class", "geometry"]].to_crs(CRS_M)
bld = bld[bld.geometry.type.isin(["Polygon", "MultiPolygon"])][["geometry"]].to_crs(CRS_M)
bld["geometry"] = shapely.make_valid(bld.geometry.values)

# study polygon in metres
tr = Transformer.from_crs(4326, CRS_M, always_xy=True)
corners = [tr.transform(x, y) for x, y in [(W, S), (E, S), (E, N_), (W, N_)]]
study = shapely.Polygon(corners)
roads["geometry"] = roads.geometry.intersection(study)
roads = roads[~roads.geometry.is_empty]
bld["geometry"] = bld.geometry.intersection(study)
bld = bld[~bld.geometry.is_empty]
print(f"after clip: {len(roads)} road segments, {len(bld)} buildings")

# ---------- box counting ----------
NPX = 128                      # raster resolution inside every cell
BOXES = [64, 32, 16, 8, 4, 2]  # box edge in pixels  -> cell/2, /4, /8, /16, /32, /64


def box_count_dim(mask):
    """fractal dimension of a boolean NPXxNPX raster. returns (D, R2, n_pixels_set)"""
    if mask.sum() < 4:
        return np.nan, np.nan, int(mask.sum())
    xs, ys = [], []
    for b in BOXES:
        n = mask.reshape(NPX // b, b, NPX // b, b).any(axis=(1, 3)).sum()
        xs.append(np.log(NPX / b))
        ys.append(np.log(n))
    xs, ys = np.array(xs), np.array(ys)
    slope, icpt = np.polyfit(xs, ys, 1)
    pred = slope * xs + icpt
    ss_res = ((ys - pred) ** 2).sum()
    ss_tot = ((ys - ys.mean()) ** 2).sum()
    return float(slope), float(1 - ss_res / ss_tot) if ss_tot > 0 else np.nan, int(mask.sum())


def line_mask(lines, x0, y0, px):
    m = np.zeros((NPX, NPX), bool)
    if len(lines) == 0:
        return m
    geoms = shapely.segmentize(np.asarray(lines, dtype=object), px / 3)
    xy = shapely.get_coordinates(geoms)
    if len(xy) == 0:
        return m
    ix = np.clip(((xy[:, 0] - x0) / px).astype(int), 0, NPX - 1)
    iy = np.clip(((xy[:, 1] - y0) / px).astype(int), 0, NPX - 1)
    m[iy, ix] = True
    return m


def poly_mask(polys, x0, y0, px):
    m = np.zeros((NPX, NPX), bool)
    if len(polys) == 0:
        return m
    u = shapely.union_all(np.asarray(polys, dtype=object))
    cx = x0 + (np.arange(NPX) + 0.5) * px
    cy = y0 + (np.arange(NPX) + 0.5) * px
    X, Y = np.meshgrid(cx, cy)
    m = shapely.contains_xy(u, X.ravel(), Y.ravel()).reshape(NPX, NPX)
    return m


minx, miny, maxx, maxy = study.bounds
# whole cells only, centred in the study area (no thin sliver cells at the edges)
ncol = int(np.floor((maxx - minx) / cell_m))
nrow = int(np.floor((maxy - miny) / cell_m))
minx += ((maxx - minx) - ncol * cell_m) / 2
miny += ((maxy - miny) - nrow * cell_m) / 2
print(f"grid {ncol} x {nrow} @ {cell_m:.0f} m")

r_tree = shapely.STRtree(roads.geometry.values)
b_tree = shapely.STRtree(bld.geometry.values)
r_geom = roads.geometry.values
r_cls = roads["road_class"].values
b_geom = bld.geometry.values

back = Transformer.from_crs(CRS_M, 4326, always_xy=True)
rows = []
geoms_out = []
for j in range(nrow):
    for i in range(ncol):
        x0, y0 = minx + i * cell_m, miny + j * cell_m
        cell = box(x0, y0, x0 + cell_m, y0 + cell_m)
        px = cell_m / NPX
        coverage = cell.intersection(study).area / cell.area

        ridx = r_tree.query(cell, predicate="intersects")
        rg = shapely.intersection(r_geom[ridx], cell)
        keep = ~shapely.is_empty(rg)
        rg, rc = rg[keep], r_cls[ridx][keep]
        # intersection of a line with a box can give multi/collections - fine for coordinate sampling
        m_all = line_mask(rg, x0, y0, px)
        m_thk = line_mask(rg[rc == "thick"], x0, y0, px)
        m_thn = line_mask(rg[rc == "thin"], x0, y0, px)

        bidx = b_tree.query(cell, predicate="intersects")
        bg = shapely.intersection(b_geom[bidx], cell)
        bg = bg[~shapely.is_empty(bg)]
        m_b = poly_mask(bg, x0, y0, px)

        rfd, rfd_r2, rpx = box_count_dim(m_all)
        rfd_thk, _, thk_px = box_count_dim(m_thk)
        rfd_thn, _, thn_px = box_count_dim(m_thn)
        bfd, bfd_r2, bpx = box_count_dim(m_b)

        len_all = float(shapely.length(rg).sum())
        len_thk = float(shapely.length(rg[rc == "thick"]).sum())
        len_thn = float(shapely.length(rg[rc == "thin"]).sum())
        built_frac = float(m_b.mean())

        c = cell.centroid
        geoms_out.append(cell)
        lon, lat = back.transform(c.x, c.y)
        rows.append(dict(
            grid_id=f"R{j:02d}C{i:02d}", row=j, col=i, lat=round(lat, 6), long=round(lon, 6),
            coverage=round(coverage, 3),
            road_length_m=round(len_all, 1), thick_length_m=round(len_thk, 1), thin_length_m=round(len_thn, 1),
            thick_share=round(len_thk / len_all, 3) if len_all > 0 else np.nan,
            n_road_segments=int(len(rg)), n_buildings=int(len(bg)),
            built_fraction=round(built_frac, 3),
            RFD=round(rfd, 4) if rfd == rfd else np.nan, RFD_R2=round(rfd_r2, 3) if rfd_r2 == rfd_r2 else np.nan,
            RFD_thick=round(rfd_thk, 4) if rfd_thk == rfd_thk else np.nan,
            RFD_thin=round(rfd_thn, 4) if rfd_thn == rfd_thn else np.nan,
            BFD=round(bfd, 4) if bfd == bfd else np.nan, BFD_R2=round(bfd_r2, 3) if bfd_r2 == bfd_r2 else np.nan,
        ))

df = pd.DataFrame(rows)
# reliability flags (values are kept, but flagged so they are not over-interpreted)
SCALE = (cell_m / 500.0) ** 2  # thresholds scale with cell area
df["low_data_road"] = df["road_length_m"] < 1000 * SCALE
df["low_data_built"] = df["n_buildings"] < 20 * SCALE
THICK_SHARE_CUT = 0.40  # arbitrary 2-region split at cell level (adjustable)
df["cell_road_class"] = np.where(df["thick_share"] >= THICK_SHARE_CUT, "Thick-dominant", "Thin-dominant")
gdf = gpd.GeoDataFrame(df.copy(), geometry=geoms_out, crs=CRS_M).to_crs(4326)
gdf.to_file(os.path.join(outdir, f"grid_fractal_{int(cell_m)}m.geojson"), driver="GeoJSON")
df.to_csv(os.path.join(outdir, f"grid_fractal_{int(cell_m)}m.csv"), index=False)
roads.to_crs(4326).to_file(os.path.join(outdir, "roads_classified.geojson"), driver="GeoJSON")
print(df.describe().T[["count", "mean", "min", "max"]])
