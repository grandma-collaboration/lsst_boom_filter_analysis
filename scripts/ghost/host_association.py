"""
host_association.py
-------------------
Three independent host association methods per transient:

1. GHOST        - astro_ghost PS1-based DLR (full sky, deep).
2. NED full DLR - astroquery remote cone search on ALL of NED (not just
                  NED-LVS local volume) + per-object diameter fetch,
                  then the shared DLR math below.
3. Legacy DLR   - NOIRLab DataLab TAP query on ls_dr10.tractor,
                  Tractor morphology, same DLR math. No sweep files.

NED DLR is primary for redshift (spectroscopic catalog, full depth).
Legacy DLR adds Tractor morphology depth.
GHOST is the PS1 second opinion.

Dependencies (pip install):
    astroquery>=0.4.7
    astro-datalab        ->  pip install --ignore-installed astro-datalab
"""

from __future__ import annotations

import math
import os
import time
import warnings
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.cosmology import FlatLambdaCDM
from astropy.table import Table
import astropy.units as u
from astropy.units import UnitsWarning

warnings.simplefilter("ignore", UnitsWarning)
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=pd.errors.SettingWithCopyWarning)
pd.options.mode.chained_assignment = None

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

COSMO = FlatLambdaCDM(H0=70, Om0=0.3)

MAX_DLR                    = 5.0
NEDLVS_SEARCH_RADIUS_ARCSEC = 120.0
MIN_AXIS_ARCSEC            = 0.5
AGREEMENT_THRESHOLD_ARCSEC = 3.0

NON_HOST_OBJTYPES = frozenset({"QSO", "AbLS", "EmLS", "EmObj", "Q_Lens", "G_Lens"})


def _non_host_mask(objtypes: pd.Series, non_host_types=NON_HOST_OBJTYPES) -> pd.Series:
    """True for NED rows that cannot host a transient: the NON_HOST_OBJTYPES list,
    any star (every NED stellar type contains '*': *, **, V*, WD*, *Cl, ...), and
    Galactic objects (NED prefixes those with '!'). Stars were not excluded before,
    so a star's velocity could become the host redshift."""
    t = objtypes.astype(str).str.strip()
    return t.isin(non_host_types) | t.str.contains("*", regex=False) | t.str.startswith("!")

_NED_RETRY_WAIT  = 2
_NED_MAX_RETRIES = 2

# Optional speed-up, OFF by default. NED sizes are fetched one
# request per galaxy, so deep fields (COSMOS: ~1,000 galaxies in 120") take hours.
# A caller can set these to fetch sizes only for the nearest galaxies.

NED_DIAM_MAX_SEP_ARCSEC = None   
NED_DIAM_MAX_CANDIDATES = None 

# Optional BOOM-matched Legacy DLR, OFF by default.
# True -> Legacy galaxy sizes follow BOOM (boom-astro/boom main, src/utils/host/
# sersic.rs, ellipse.rs, catalog.rs, config.rs, dlr.rs) instead of 1.5 x r_e:
#   size  = semi-major axis where the Tractor Sersic model reaches 25 mag/arcsec^2
#           in r (keeps r_e when the galaxy is fainter than that at its centre)
#   angle = 0.5*atan2(e2, e1) east of north (checked against HyperLEDA PAs)
#   floor = 0.05" on both axes; PSF/DUP and marginal REX rows dropped; b/a < 0.05 dropped
LEGACY_DLR_BOOM = False
# Rows the Legacy box query may return. 500 truncates deep fields (COSMOS boxes hold
# ~750-800 galaxies, so ~1/3 are dropped at random); a caller can raise it.
LEGACY_QUERY_LIMIT = 500


# ---------------------------------------------------------------------------
# Scalar helpers
# ---------------------------------------------------------------------------

def _clean_scalar(value: Any) -> Any:
    if value is None:
        return None
    try:
        if pd.isna(value):
            return None
    except Exception:
        pass
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        value = float(value)
        return value if math.isfinite(value) else None
    if isinstance(value, (bytes, bytearray)):
        return value.decode("utf-8", errors="ignore").strip()
    return value


def _as_float(value: Any) -> float | None:
    value = _clean_scalar(value)
    if value is None:
        return None
    try:
        x = float(value)
    except Exception:
        return None
    return x if math.isfinite(x) else None


# ---------------------------------------------------------------------------
# Shared DLR math  (same as dlr_ranker.py + nedlvs_morphology.py
#                   + tractor_morphology.py)
# ---------------------------------------------------------------------------

def _diam_to_ellipse(
    diam_arcsec: Any,
    diam_ba: Any,
    diam_pa_ned: Any,
) -> tuple[float, float, float, float]:
    """
    NED D25 diameter columns -> (a_arcsec, b_arcsec, axis_ratio, pa_deg).
    diam_arcsec is the full major-axis diameter in arcsec; a = diam/2.
    NED PA is North-through-East; DLR wants CCW-from-East: pa = (90 - pa_ned) % 180.
    """
    diam = _as_float(diam_arcsec)
    if diam is None or diam <= 0:
        return np.nan, np.nan, np.nan, np.nan
    a = diam / 2.0
    q = _as_float(diam_ba)
    q = float(np.clip(q, 1e-4, 1.0)) if q is not None else 1.0
    b = max(a * q, MIN_AXIS_ARCSEC)
    pa_ned = _as_float(diam_pa_ned)
    pa = 0.0 if pa_ned is None else float((90.0 - pa_ned) % 180.0)
    return a, b, q, pa


def _tractor_to_ellipse(
    shape_r: Any,
    shape_e1: Any,
    shape_e2: Any,
    size_scale: float = 1.5,
) -> tuple[float, float, float, float]:
    """
    Tractor shape columns -> (a_arcsec, b_arcsec, axis_ratio, pa_deg).
    Mirrors tractor_morphology.py exactly.

    shape_r is the Tractor effective (half-light) radius in arcsec for
    de Vaucouleurs / exponential profiles, or the Gaussian sigma for REX.
    size_scale = 1.5 converts effective radius to a DLR semi-major axis
    proxy: for an exponential disk, the D25 isophote (used in NED) falls
    at roughly 1.5 * r_eff, making the scale factors comparable between
    the NED-DLR and Legacy-DLR methods (Sérsic n=1 geometry; ref:
    Graham & Driver 2005, PASA 22, 118).  Value inherited from
    dlr_ranker.py / tractor_morphology.py.
    """
    r = _as_float(shape_r)
    if r is None:
        return np.nan, np.nan, np.nan, np.nan
    e1 = _as_float(shape_e1) or 0.0
    e2 = _as_float(shape_e2) or 0.0
    e  = min(max(math.hypot(e1, e2), 0.0), 0.999999)
    q  = float(np.clip((1 - e) / (1 + e), 1e-6, 1.0))
    pa = float(math.degrees(0.5 * math.atan2(e2, e1)) % 180.0)
    a  = r * size_scale
    b  = max(a * q, MIN_AXIS_ARCSEC)
    return a, b, q, pa


# ---------------------------------------------------------------------------
# BOOM-matched Legacy size (used only when LEGACY_DLR_BOOM is True).
# Numbers are BOOM's HostGalaxyConfig defaults (config.rs).
# ---------------------------------------------------------------------------

BOOM_MIN_AXIS_ARCSEC  = 0.05   # min_axis_arcsec: floor on both semi-axes
BOOM_MIN_AXIS_RATIO   = 0.05   # min_axis_ratio: below this b/a the fit failed, row dropped
BOOM_ISOPHOTE_MAG     = 25.0   # isophote_mag: r-band mag/arcsec^2 (Legacy fluxes are AB)
BOOM_REX_MIN_SHAPE_R  = 0.3    # rex_min_shape_r_arcsec
BOOM_REX_MIN_SNR      = 5.0    # rex_min_snr (r band)
BOOM_REX_MAX_FRACFLUX = 0.5    # rex_max_fracflux (r band)
BOOM_STAR_TYPES       = frozenset({"PSF", "DUP"})   # star_type_values


def _sersic_index_for_type(objtype: Any, fitted: Any) -> float | None:
    """sersic.rs sersic_index_for_type: only SER carries a fitted index."""
    t = str(objtype).strip().upper() if objtype is not None else ""
    if t in ("REX", "EXP"):
        return 1.0
    if t == "DEV":
        return 4.0
    if t == "SER":
        n = _as_float(fitted)
        return n if (n is not None and math.isfinite(n) and n > 0) else None
    return None


def _isophotal_semi_major(
    r_e: float, q: float, n: float, m_tot: float, mu_target: float = BOOM_ISOPHOTE_MAG,
) -> float | None:
    """
    sersic.rs isophotal_semi_major: the semi-major axis (arcsec) at which a Sersic
    model of half-light semi-major axis r_e, axis ratio q, index n and total
    magnitude m_tot reaches surface brightness mu_target. None when the model is
    already fainter than mu_target at its centre (no size is invented).

      <mu>_e = m_tot + 2.5 log10(2 pi r_e^2 q)               mean SB inside r_e
      mu_e   = <mu>_e + 2.5 log10(n e^b Gamma(2n) / b^(2n))   SB at r_e
      mu(r)  = mu_e + (2.5 b / ln10) [(r/r_e)^(1/n) - 1]      Sersic profile
    solved for mu(a25) = mu_target. b = b_n is the median of Gamma(2n), the same
    number BOOM finds by bisection (sersic_b).
    """
    from scipy.special import gammaincinv, gammaln

    if not (
        all(math.isfinite(v) for v in (r_e, q, n, m_tot))
        and r_e > 0 and q > 0 and n > 0
    ):
        return None
    b = float(gammaincinv(2.0 * n, 0.5))
    mu_bar_e = m_tot + 2.5 * math.log10(2.0 * math.pi * r_e * r_e * q)
    # log10(n e^b b^(-2n) Gamma(2n)), in logs so large n cannot overflow
    log10_f = (math.log(n) + b - 2.0 * n * math.log(b) + float(gammaln(2.0 * n))) / math.log(10.0)
    mu_e = mu_bar_e + 2.5 * log10_f
    bracket = 1.0 + (math.log(10.0) / (2.5 * b)) * (mu_target - mu_e)
    if bracket <= 0:
        return None
    a25 = r_e * bracket ** n
    return a25 if (math.isfinite(a25) and a25 > 0) else None


def _boom_marginal_rex(row: Any, shape_r: float) -> bool:
    """catalog.rs rejected_as_marginal_rex: too small, too faint, or a blend fragment."""
    if shape_r < BOOM_REX_MIN_SHAPE_R:
        return True
    # Absent columns mean no S/N was measured, which must not reject the row.
    flux = _as_float(row.get("flux_r"))
    ivar = _as_float(row.get("flux_ivar_r"))
    if flux is not None and ivar is not None:
        if flux <= 0:
            return True
        if ivar > 0 and flux * math.sqrt(ivar) < BOOM_REX_MIN_SNR:
            return True
    frac = _as_float(row.get("fracflux_r"))
    return frac is not None and frac > BOOM_REX_MAX_FRACFLUX


def _tractor_row_to_boom_ellipse(row: Any) -> tuple[float, float, float, float, bool]:
    """
    One ls_dr10.tractor row -> (a_arcsec, b_arcsec, axis_ratio, pa_deg, size_is_isophotal),
    following BOOM's galaxy_from_ls_dr10 + Ellipse::from_tractor / scaled_to_semi_major /
    from_candidate. All-NaN when BOOM would drop the row.

    pa_deg is returned in THIS file's DLR frame (counter-clockwise from East, as
    _compute_dlr_values expects); BOOM's angle is 0.5*atan2(e2, e1) east of north,
    so pa_deg = 90 - that.
    """
    dropped = (np.nan, np.nan, np.nan, np.nan, False)
    objtype = str(row.get("type") or "").strip().upper()
    if objtype in BOOM_STAR_TYPES:
        return dropped
    shape_r = _as_float(row.get("shape_r"))
    if shape_r is None or not (shape_r > 0):
        return dropped
    if objtype == "REX" and _boom_marginal_rex(row, shape_r):
        return dropped
    e1 = _as_float(row.get("shape_e1"))
    e2 = _as_float(row.get("shape_e2"))
    e1 = 0.0 if e1 is None else e1
    e2 = 0.0 if e2 is None else e2
    if not (math.isfinite(e1) and math.isfinite(e2)):
        return dropped

    # Ellipse::from_tractor
    e = min(math.hypot(e1, e2), 0.999)
    q = (1.0 - e) / (1.0 + e)
    a = shape_r
    b = max(a * q, BOOM_MIN_AXIS_ARCSEC)
    axis_ratio = b / a
    pa_east_of_north = math.degrees(0.5 * math.atan2(e2, e1))
    if not math.isfinite(axis_ratio) or axis_ratio < BOOM_MIN_AXIS_RATIO:
        return dropped   # bounded_axis_ratio: a shape fit failure

    # isophotal_semi_major_for + scaled_to_semi_major
    n = _sersic_index_for_type(objtype, row.get("sersic"))
    flux_r = _as_float(row.get("flux_r"))
    m_tot = (22.5 - 2.5 * math.log10(flux_r)) if (flux_r is not None and math.isfinite(flux_r) and flux_r > 0) else None
    a25 = _isophotal_semi_major(a, axis_ratio, n, m_tot) if (n is not None and m_tot is not None) else None
    size_is_isophotal = a25 is not None
    if size_is_isophotal:
        a = max(a25, BOOM_MIN_AXIS_ARCSEC)
        b = max(a * axis_ratio, BOOM_MIN_AXIS_ARCSEC)

    # Ellipse::from_candidate -> Ellipse::new: floor b, larger axis becomes a
    b = max(b, BOOM_MIN_AXIS_ARCSEC)
    if b > a:
        a, b = b, a
    return a, b, b / a, (90.0 - pa_east_of_north) % 180.0, size_is_isophotal


def _compute_dlr_values(
    transient_ra_deg: float,
    transient_dec_deg: float,
    cand: pd.DataFrame,
    min_axis_arcsec: float | None = None,
) -> pd.DataFrame:
    required = {"ra", "dec", "a_arcsec", "b_arcsec", "pa_deg"}
    if required - set(cand.columns):
        return cand.copy()
    min_axis = MIN_AXIS_ARCSEC if min_axis_arcsec is None else min_axis_arcsec
    out      = cand.copy()
    host_ra  = out["ra"].to_numpy(dtype=float)
    host_dec = out["dec"].to_numpy(dtype=float)
    a        = out["a_arcsec"].to_numpy(dtype=float)
    b        = out["b_arcsec"].to_numpy(dtype=float)
    pa       = out["pa_deg"].to_numpy(dtype=float)
    a = np.where(np.isfinite(a), np.maximum(a, min_axis), np.nan)
    b = np.where(np.isfinite(b), np.maximum(b, min_axis), np.nan)
    cos_dec  = np.cos(np.deg2rad(host_dec))
    dx       = (transient_ra_deg  - host_ra)  * cos_dec * 3600.0
    dy       = (transient_dec_deg - host_dec) * 3600.0
    theta    = np.deg2rad(pa)
    x_major  = dx * np.cos(theta) + dy * np.sin(theta)
    y_minor  = -dx * np.sin(theta) + dy * np.cos(theta)
    out["dlr"]                = np.sqrt((x_major / a) ** 2 + (y_minor / b) ** 2)
    out["transient_sep_arcsec"] = np.sqrt(dx**2 + dy**2)
    return out


def _rank_candidates(
    ra_deg: float,
    dec_deg: float,
    cand: pd.DataFrame,
    non_host_types: frozenset | None = None,
    min_axis_arcsec: float | None = None,
) -> pd.DataFrame:
    work = cand.copy()
    if non_host_types and "objtype" in work.columns:
        work = work.loc[
            ~_non_host_mask(work["objtype"], non_host_types)
        ].copy()
    # Drop Tractor PSF sources (stars)
    if "type" in work.columns:
        work = work.loc[
            work["type"].astype(str).str.upper().str.strip() != "PSF"
        ].copy()
    if len(work) == 0:
        return work
    work = _compute_dlr_values(ra_deg, dec_deg, work, min_axis_arcsec=min_axis_arcsec)
    num  = work.select_dtypes(include=[np.number]).columns
    work[num] = work[num].replace([np.inf, -np.inf], np.nan)
    work = work.loc[np.isfinite(work["dlr"]) & (work["dlr"] <= MAX_DLR)].copy()
    if len(work) == 0:
        return work
    work = work.sort_values(
        ["dlr", "transient_sep_arcsec"], kind="stable"
    ).reset_index(drop=True)
    work["dlr_rank"] = np.arange(1, len(work) + 1, dtype=int)
    return work


# ---------------------------------------------------------------------------
# NED full-sky remote query  (replaces local NED-LVS FITS)
# ---------------------------------------------------------------------------

def _ned_cone_search(
    ra_deg: float,
    dec_deg: float,
    radius_arcsec: float,
) -> pd.DataFrame:
    """
    Cone search ALL of NED via astroquery.ipac.ned.
    Returns DataFrame with: ra, dec, objname, objtype, z, transient_sep_arcsec.
    """
    try:
        from astroquery.ipac.ned import Ned
        coord = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg)

        for attempt in range(_NED_MAX_RETRIES + 1):
            try:
                result = Ned.query_region(coord, radius=radius_arcsec * u.arcsec)
                break
            except Exception as e:
                if attempt < _NED_MAX_RETRIES:
                    time.sleep(_NED_RETRY_WAIT)
                else:
                    print(f"    [NED] cone query failed after retries: {e}")
                    return pd.DataFrame()

        if result is None or len(result) == 0:
            return pd.DataFrame()

        df = result.to_pandas()
        # Normalise column names — astroquery uses spaces and mixed case
        df.columns = [str(c).strip().lower().replace(" ", "_") for c in df.columns]

        # Map to standard names used downstream.
        # Use exact matches first to avoid substring collisions (e.g. "decra"
        # or "search_offset_arcsec" falsely matching "ra"/"dec").
        _RA_EXACT  = {"ra", "ra(deg)", "ra_deg"}
        _DEC_EXACT = {"dec", "dec(deg)", "dec_deg"}
        col_map: dict[str, str] = {}
        for col in df.columns:
            lc = col.lower()
            if lc in _RA_EXACT and "ra" not in col_map.values():
                col_map[col] = "ra"
            elif lc in _DEC_EXACT and "dec" not in col_map.values():
                col_map[col] = "dec"
            elif "object_name" in lc and "objname" not in col_map.values():
                col_map[col] = "objname"
            elif lc in ("type", "object_type") and "objtype" not in col_map.values():
                col_map[col] = "objtype"
            elif lc == "redshift" and "z" not in col_map.values():
                col_map[col] = "z"
            elif lc == "velocity" and "z" not in col_map.values():
                # NED velocity column is cz in km/s — store raw, convert below
                col_map[col] = "z_raw_velocity"
        # Fallback: suffix match for astroquery versions that prefix with table name
        if "ra" not in col_map.values():
            for col in df.columns:
                lc = col.lower()
                if (lc.endswith("_ra") or lc.endswith(".ra")) and "ra" not in col_map.values():
                    col_map[col] = "ra"
        if "dec" not in col_map.values():
            for col in df.columns:
                lc = col.lower()
                if (lc.endswith("_dec") or lc.endswith(".dec")) and "dec" not in col_map.values():
                    col_map[col] = "dec"
        df = df.rename(columns=col_map)

        # Convert NED cz (km/s) to dimensionless redshift where needed.
        # NED returns velocity in km/s when the column is called "Velocity"
        # rather than "Redshift".  Values > 600 km/s are unambiguously cz.
        _C_KMS = 299792.458
        if "z_raw_velocity" in df.columns and "z" not in df.columns:
            df["z"] = df["z_raw_velocity"].apply(
                lambda v: _as_float(v) / _C_KMS if _as_float(v) is not None else None
            )
        elif "z_raw_velocity" in df.columns:
            # z column already exists (from redshift col); drop the raw velocity
            df = df.drop(columns=["z_raw_velocity"])

        for col in df.columns:
            if df[col].dtype == object:
                df[col] = df[col].apply(_clean_scalar)

        if "ra" not in df.columns or "dec" not in df.columns:
            return pd.DataFrame()

        center = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg)
        coords = SkyCoord(
            ra=df["ra"].to_numpy(dtype=float) * u.deg,
            dec=df["dec"].to_numpy(dtype=float) * u.deg,
        )
        df["transient_sep_arcsec"] = center.separation(coords).arcsec
        return df.sort_values("transient_sep_arcsec").reset_index(drop=True)

    except Exception as e:
        print(f"    [NED] cone search error: {e}")
        return pd.DataFrame()


def _ned_fetch_diameters(objname: str) -> tuple[float, float, float]:
    """
    Fetch D25 diameter, axis ratio, PA for one NED object via astroquery.
    Returns (diam_arcsec, diam_ba, diam_pa_ned) or (nan, nan, nan).
    NED diameter table gives values in arcmin -> we convert to arcsec here.
    """
    try:
        from astroquery.ipac.ned import Ned
        for attempt in range(_NED_MAX_RETRIES + 1):
            try:
                table = Ned.get_table(objname, table="diameters")
                break
            except Exception as e:
                if attempt < _NED_MAX_RETRIES:
                    time.sleep(_NED_RETRY_WAIT)
                else:
                    return np.nan, np.nan, np.nan

        if table is None or len(table) == 0:
            return np.nan, np.nan, np.nan

        df = table.to_pandas()
        df.columns = [str(c).strip().lower().replace(" ", "_") for c in df.columns]

        # Find columns by partial name match.
        # NED diameter tables can return either an explicit axis ratio column
        # ("axis_ratio") or separate major/minor axis columns.  We must
        # distinguish the two cases: if we have a raw minor-axis value we must
        # divide by the major axis ourselves; if we have an axis ratio we use
        # it directly.  Mixing the two would feed a raw arcmin value where a
        # dimensionless ratio is expected and produce wildly wrong ellipses.
        diam_col      = next((c for c in df.columns if "major_axis" in c), None)
        ratio_col     = next((c for c in df.columns if "axis_ratio" in c), None)
        minor_col     = next((c for c in df.columns if "minor_axis" in c and "major" not in c), None)
        pa_col        = next((c for c in df.columns if "position_angle" in c or c == "pa"), None)

        if diam_col is None:
            return np.nan, np.nan, np.nan

        diam_arcmin = _as_float(df[diam_col].median())
        if diam_arcmin is None:
            return np.nan, np.nan, np.nan

        # Derive axis ratio (b/a, dimensionless, in [0,1]).
        ba: float | None = None
        if ratio_col is not None:
            ba = _as_float(df[ratio_col].median())
        elif minor_col is not None:
            # minor_axis is in the same unit as major_axis (arcmin); ratio = minor/major
            minor_arcmin = _as_float(df[minor_col].median())
            if minor_arcmin is not None and diam_arcmin > 0:
                ba = minor_arcmin / diam_arcmin
                ba = float(np.clip(ba, 1e-4, 1.0))

        pa = _as_float(df[pa_col].median()) if pa_col else None

        # NED diameter table is in arcmin -> convert to arcsec
        diam_arcsec = diam_arcmin * 60.0

        return diam_arcsec, ba if ba is not None else np.nan, pa if pa is not None else np.nan

    except Exception:
        return np.nan, np.nan, np.nan


def _add_ned_ellipse_columns(cand: pd.DataFrame) -> pd.DataFrame:
    """
    For each candidate, fetch NED diameters and add ellipse columns.
    One HTTP call per candidate — sequential with retries.

    Fallback ellipse for missing diameters
    --------------------------------------
    NED's per-object diameter table is sparsely populated: the majority of
    faint/uncatalogued galaxies in a cone have no D25 entry at all.  When the
    fetch returns nan, _rank_candidates drops those rows (non-finite DLR),
    which is why the cone finds 60+ candidates but ranks 0.

    For objects with no measured diameter we assign a point-source fallback:
    a circular aperture of POINT_SOURCE_FALLBACK_ARCSEC radius.  This is
    consistent with the NED-LVS (local FITS catalog) version, which gets ellipse data
    inline from the FITS catalog and only falls back for objects with missing
    diam columns.  With the fallback, objects without D25 get a finite DLR
    (= transient_sep_arcsec / FALLBACK), rank behind properly-sized objects,
    but are not silently discarded.
    """
    # Nominal circular aperture for objects with no NED diameter entry.
    # 3 arcsec is a conservative lower bound for any galaxy NED catalogs.
    POINT_SOURCE_FALLBACK_ARCSEC = 3.0

    out = cand.copy()
    a_list, b_list, q_list, pa_list, has_diam_list = [], [], [], [], []

    for _, row in out.iterrows():
        objname = _clean_scalar(row.get("objname", ""))
        if objname:
            diam, ba, pa_ned = _ned_fetch_diameters(str(objname))
        else:
            diam = ba = pa_ned = np.nan

        a, b, q, pa = _diam_to_ellipse(diam, ba, pa_ned)

        if not math.isfinite(a):
            # Fallback: treat as circular point source
            a  = POINT_SOURCE_FALLBACK_ARCSEC
            b  = POINT_SOURCE_FALLBACK_ARCSEC
            q  = 1.0
            pa = 0.0
            has_diam_list.append(False)
        else:
            has_diam_list.append(True)

        a_list.append(a)
        b_list.append(b)
        q_list.append(q)
        pa_list.append(pa)

    out["a_arcsec"]     = a_list
    out["b_arcsec"]     = b_list
    out["axis_ratio"]   = q_list
    out["pa_deg"]       = pa_list
    out["has_ned_diam"] = has_diam_list  # True = real D25, False = fallback
    return out


def _run_ned_dlr(transient_name: str, ra_deg: float, dec_deg: float) -> dict[str, Any]:
    """Full NED remote cone -> fetch diameters -> DLR rank."""
    print(f"    [NED] querying ra={ra_deg:.5f} dec={dec_deg:.5f} r={NEDLVS_SEARCH_RADIUS_ARCSEC}\"")
    cone = _ned_cone_search(ra_deg, dec_deg, NEDLVS_SEARCH_RADIUS_ARCSEC)

    if len(cone) == 0:
        return {"nedlvs_status": "no_candidates", "nedlvs_n_candidates": 0, "nedlvs_n_ranked": 0}

    # Filter non-host types before fetching diameters (saves HTTP calls)
    work = cone.copy()
    if "objtype" in work.columns:
        work = work.loc[
            ~_non_host_mask(work["objtype"])
        ].copy()

    if len(work) == 0:
        return {
            "nedlvs_status": "no_host_within_dlr_threshold",
            "nedlvs_n_candidates": int(len(cone)),
            "nedlvs_n_ranked": 0,
        }

    n_all = len(work)
    if "transient_sep_arcsec" in work.columns:
        if NED_DIAM_MAX_SEP_ARCSEC is not None:
            work = work.loc[work["transient_sep_arcsec"] <= NED_DIAM_MAX_SEP_ARCSEC].copy()
        if NED_DIAM_MAX_CANDIDATES is not None and len(work) > NED_DIAM_MAX_CANDIDATES:
            work = work.sort_values("transient_sep_arcsec", kind="stable").head(NED_DIAM_MAX_CANDIDATES).copy()
    if len(work) == 0:
        return {
            "nedlvs_status": "no_host_within_dlr_threshold",
            "nedlvs_n_candidates": int(len(cone)),
            "nedlvs_n_ranked": 0,
        }
    if len(work) < n_all:
        print(f"    [NED] {n_all} galaxy candidates, fetching diameters for the nearest {len(work)} …")
    else:
        print(f"    [NED] {len(work)} galaxy candidates, fetching diameters …")
    work = _add_ned_ellipse_columns(work)
    ranked = _rank_candidates(ra_deg, dec_deg, work, non_host_types=NON_HOST_OBJTYPES)

    if len(ranked) == 0:
        return {
            "nedlvs_status": "no_host_within_dlr_threshold",
            "nedlvs_n_candidates": int(len(cone)),
            "nedlvs_n_ranked": 0,
        }

    best = ranked.iloc[0]
    # nedlvs_diam_arcsec: reconstruct from a_arcsec (semi-major axis = diam/2)
    # so the CSV column matches what was fetched from NED's diameter table.
    a_best = _as_float(best.get("a_arcsec"))
    diam_best = (a_best * 2.0) if a_best is not None else None
    return {
        "nedlvs_status":        "host_found",
        "nedlvs_ra":            _as_float(best.get("ra")),
        "nedlvs_dec":           _as_float(best.get("dec")),
        "nedlvs_dlr":           _as_float(best.get("dlr")),
        "nedlvs_sep_arcsec":    _as_float(best.get("transient_sep_arcsec")),
        "nedlvs_objname":       _clean_scalar(best.get("objname")),
        "nedlvs_objtype":       _clean_scalar(best.get("objtype")),
        "nedlvs_z":             _as_float(best.get("z")),
        "nedlvs_distmpc":       None,   # NED cone doesn't return distmpc directly
        "nedlvs_diam_arcsec":   diam_best,
        "nedlvs_a_arcsec":      a_best,
        "nedlvs_b_arcsec":      _as_float(best.get("b_arcsec")),
        "nedlvs_axis_ratio":    _as_float(best.get("axis_ratio")),
        "nedlvs_pa_deg":        _as_float(best.get("pa_deg")),
        "nedlvs_n_candidates":  int(len(cone)),
        "nedlvs_n_ranked":      int(len(ranked)),
    }


# ---------------------------------------------------------------------------
# Legacy Survey DR10 remote query via NOIRLab DataLab TAP
# ---------------------------------------------------------------------------

def _run_legacy_dlr(transient_name: str, ra_deg: float, dec_deg: float) -> dict[str, Any]:
    """
    Query ls_dr10.tractor via NOIRLab DataLab TAP using raw HTTP (no dl client).
    Uses a box pre-filter then exact cone cut. No sweep files needed.
    """
    import requests as _requests
    import io

    radius_deg = NEDLVS_SEARCH_RADIUS_ARCSEC / 3600.0
    cos_dec    = max(abs(math.cos(math.radians(dec_deg))), 1e-6)
    ra_margin  = radius_deg / cos_dec

    ra_lo  = ra_deg  - ra_margin
    ra_hi  = ra_deg  + ra_margin
    dec_lo = dec_deg - radius_deg
    dec_hi = dec_deg + radius_deg

    if LEGACY_DLR_BOOM:
        # BOOM's size needs the Sersic index and r-band flux; REX cuts need S/N and fracflux.
        query = (
            f"SELECT ra, dec, type, shape_r, shape_e1, shape_e2, "
            f"sersic, flux_r, flux_ivar_r, fracflux_r "
            f"FROM ls_dr10.tractor "
            f"WHERE ra BETWEEN {ra_lo} AND {ra_hi} "
            f"AND dec BETWEEN {dec_lo} AND {dec_hi} "
            f"AND type != 'PSF' AND type != 'DUP' "
            f"LIMIT {int(LEGACY_QUERY_LIMIT)}"
        )
    else:
        query = (
            f"SELECT ra, dec, type, shape_r, shape_e1, shape_e2 "
            f"FROM ls_dr10.tractor "
            f"WHERE ra BETWEEN {ra_lo} AND {ra_hi} "
            f"AND dec BETWEEN {dec_lo} AND {dec_hi} "
            f"AND type != 'PSF' "
            f"LIMIT {int(LEGACY_QUERY_LIMIT)}"
        )

    try:
        resp = _requests.post(
            "https://datalab.noirlab.edu/tap/sync",
            data={"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": query},
            timeout=30,
        )
        resp.raise_for_status()
    except Exception as e:
        print(f"    [Legacy TAP] request failed: {e}")
        return {"legacy_status": f"query_error: {e}"}

    # Check for VOTable error response
    if b"QUERY_STATUS" in resp.content and b"ERROR" in resp.content:
        print(f"    [Legacy TAP] server error: {resp.text[:200]}")
        return {"legacy_status": "server_error"}

    try:
        result = pd.read_csv(io.StringIO(resp.text))
    except Exception as e:
        print(f"    [Legacy TAP] CSV parse failed: {e}")
        return {"legacy_status": f"parse_error: {e}"}

    if result is None or len(result) == 0:
        return {"legacy_status": "no_candidates", "legacy_n_candidates": 0}
    legacy_truncated = len(result) >= int(LEGACY_QUERY_LIMIT)
    if legacy_truncated:
        print(f"    [Legacy TAP] WARNING: hit the {int(LEGACY_QUERY_LIMIT)}-row limit; "
              f"some galaxies in the box were not returned")

    # Exact cone cut — compute separations once, reset index so assignment is safe
    center = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg)
    coords = SkyCoord(
        ra=result["ra"].to_numpy(dtype=float) * u.deg,
        dec=result["dec"].to_numpy(dtype=float) * u.deg,
    )
    seps = center.separation(coords).arcsec
    mask = seps <= NEDLVS_SEARCH_RADIUS_ARCSEC
    result = result.loc[mask].copy().reset_index(drop=True)
    result["transient_sep_arcsec"] = seps[mask]  # boolean-mask preserves original index
    # Re-index the separation array to match the reset DataFrame index
    result["transient_sep_arcsec"] = result["transient_sep_arcsec"].to_numpy()

    if len(result) == 0:
        return {"legacy_status": "no_candidates", "legacy_n_candidates": 0}

    # Add Tractor ellipse columns
    if LEGACY_DLR_BOOM:
        ellipse = result.apply(
            _tractor_row_to_boom_ellipse, axis=1, result_type="expand",
        )
        ellipse.columns = ["a_arcsec", "b_arcsec", "axis_ratio", "pa_deg", "size_is_isophotal"]
        result[["a_arcsec", "b_arcsec", "axis_ratio", "pa_deg", "size_is_isophotal"]] = ellipse
        min_axis = BOOM_MIN_AXIS_ARCSEC
        dlr_method = "boom_isophote_25"
    else:
        ellipse = result.apply(
            lambda row: _tractor_to_ellipse(
                row.get("shape_r"), row.get("shape_e1"), row.get("shape_e2")
            ),
            axis=1,
            result_type="expand",
        )
        ellipse.columns = ["a_arcsec", "b_arcsec", "axis_ratio", "pa_deg"]
        result[["a_arcsec", "b_arcsec", "axis_ratio", "pa_deg"]] = ellipse
        min_axis = None
        dlr_method = "fixed_1.5re"

    ranked = _rank_candidates(ra_deg, dec_deg, result, min_axis_arcsec=min_axis)
    if len(ranked) == 0:
        out = {
            "legacy_status":       "no_host_within_dlr_threshold",
            "legacy_n_candidates": int(len(result)),
            "legacy_n_ranked":     0,
        }
        if LEGACY_DLR_BOOM:
            out.update({"legacy_dlr_method": dlr_method, "legacy_truncated": legacy_truncated})
        return out

    best = ranked.iloc[0]
    out = {
        "legacy_status":       "host_found",
        "legacy_ra":           _as_float(best.get("ra")),
        "legacy_dec":          _as_float(best.get("dec")),
        "legacy_dlr":          _as_float(best.get("dlr")),
        "legacy_sep_arcsec":   _as_float(best.get("transient_sep_arcsec")),
        "legacy_type":         _clean_scalar(best.get("type")),
        "legacy_a_arcsec":     _as_float(best.get("a_arcsec")),
        "legacy_b_arcsec":     _as_float(best.get("b_arcsec")),
        "legacy_n_candidates": int(len(result)),
        "legacy_n_ranked":     int(len(ranked)),
    }
    if LEGACY_DLR_BOOM:
        out.update({
            "legacy_dlr_method":        dlr_method,
            "legacy_size_is_isophotal": bool(best.get("size_is_isophotal")),
            "legacy_shape_r":           _as_float(best.get("shape_r")),
            "legacy_pa_deg":            _as_float(best.get("pa_deg")),   # this file's DLR frame (CCW from East)
            "legacy_truncated":         legacy_truncated,
        })
    return out


# ---------------------------------------------------------------------------
# GHOST wrapper
# ---------------------------------------------------------------------------


def _run_ghost(
    transient_name: str,
    ra_deg: float,
    dec_deg: float,
    ghost_savepath: str = "/tmp/ghost_out/",
) -> dict[str, Any]:
    """Run astro_ghost for one transient and return a flat dict."""
    try:
        from astro_ghost.ghostHelperFunctions import getTransientHosts
        import shutil

        ghost_savepath = os.path.join(ghost_savepath, "")

        # Wipe any stale state from a previous failed run before starting.
        # When SkyMapper times out mid-call, GHOST writes the transient CSV
        # but never writes the hosts CSV.  The next object then crashes with
        # FileNotFoundError trying to read a file that doesn't exist in the
        # stale directory.  Clearing it first breaks the cascade.
        if os.path.exists(ghost_savepath):
            shutil.rmtree(ghost_savepath)
        os.makedirs(ghost_savepath, exist_ok=True)

        coord = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg)
        host_df = getTransientHosts(
            transientName=[str(transient_name)],
            transientCoord=[coord],
            verbose=False,
            starcut="normal",
            savepath=ghost_savepath,
            redo_search=True,
        )

        if host_df is None or len(host_df) == 0:
            return {"ghost_status": "no_host_found"}

        row = host_df.iloc[0]

        def get_first(names: list[str]) -> Any:
            for name in names:
                if name in row.index:
                    value = _clean_scalar(row[name])
                    if value is not None:
                        return value
            return None

        ghost_ra  = _as_float(get_first(["raMean", "ra", "RAJ2000", "RA"]))
        ghost_dec = _as_float(get_first(["decMean", "dec", "DEJ2000", "DEC"]))
        ghost_z   = _as_float(get_first(["z_best", "photo_z", "z_phot", "photo_z_best", "z"]))

        ghost_sep = None
        if ghost_ra is not None and ghost_dec is not None:
            c1 = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg)
            c2 = SkyCoord(ra=ghost_ra * u.deg, dec=ghost_dec * u.deg)
            ghost_sep = float(c1.separation(c2).arcsec)

        if ghost_ra is None or ghost_dec is None:
            return {"ghost_status": "host_found_missing_coordinates"}

        return {
            "ghost_status":     "host_found",
            "ghost_ra":         ghost_ra,
            "ghost_dec":        ghost_dec,
            "ghost_z":          ghost_z,
            "ghost_sep_arcsec": ghost_sep,
            "ghost_objid":      _clean_scalar(get_first(["objID", "objid", "id"])),
            "ghost_type":       _clean_scalar(get_first(["type", "sourceType", "hostType"])),
        }

    except Exception as exc:
        import traceback
        msg = str(exc)
        if "index 0 is out of bounds" in msg or "No potential hosts" in msg:
            return {"ghost_status": "no_host_found", "ghost_error": msg}
        print(f"    [GHOST] failed for {transient_name}: {msg}")
        traceback.print_exc()
        return {"ghost_status": f"error: {msg}", "ghost_error": msg}


# ---------------------------------------------------------------------------
# Agreement and derived quantities
# ---------------------------------------------------------------------------


def _angular_sep_arcsec(ra1: Any, dec1: Any, ra2: Any, dec2: Any) -> float | None:
    ra1  = _as_float(ra1)
    dec1 = _as_float(dec1)
    ra2  = _as_float(ra2)
    dec2 = _as_float(dec2)
    if None in (ra1, dec1, ra2, dec2):
        return None
    c1 = SkyCoord(ra=ra1 * u.deg, dec=dec1 * u.deg)
    c2 = SkyCoord(ra=ra2 * u.deg, dec=dec2 * u.deg)
    return float(c1.separation(c2).arcsec)


def _check_agreement(
    ghost: dict[str, Any],
    nedlvs: dict[str, Any],
    legacy: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """
    Compute pairwise angular separations between all three host methods.

    GHOST vs NED, NED vs Legacy, and GHOST vs Legacy are all reported so that
    the agreement field is meaningful regardless of which method is ultimately
    chosen as host_source.  Previously only GHOST vs NED was compared, which
    gave no information when host_source == LEGACY_DLR.
    """
    legacy = legacy or {}
    g_ok = ghost.get("ghost_status") == "host_found"
    n_ok = nedlvs.get("nedlvs_status") == "host_found"
    l_ok = legacy.get("legacy_status") == "host_found"

    def _sep_or_none(ra1_key_pair, ra2_key_pair):
        ra1, dec1 = ra1_key_pair
        ra2, dec2 = ra2_key_pair
        return _angular_sep_arcsec(ra1, dec1, ra2, dec2)

    ghost_nedlvs_sep  = None
    ghost_legacy_sep  = None
    nedlvs_legacy_sep = None

    if g_ok and n_ok:
        ghost_nedlvs_sep = _angular_sep_arcsec(
            ghost.get("ghost_ra"), ghost.get("ghost_dec"),
            nedlvs.get("nedlvs_ra"), nedlvs.get("nedlvs_dec"),
        )
    if g_ok and l_ok:
        ghost_legacy_sep = _angular_sep_arcsec(
            ghost.get("ghost_ra"), ghost.get("ghost_dec"),
            legacy.get("legacy_ra"), legacy.get("legacy_dec"),
        )
    if n_ok and l_ok:
        nedlvs_legacy_sep = _angular_sep_arcsec(
            nedlvs.get("nedlvs_ra"), nedlvs.get("nedlvs_dec"),
            legacy.get("legacy_ra"), legacy.get("legacy_dec"),
        )

    # Summarise overall agreement across all available method pairs.
    available_seps = [s for s in (ghost_nedlvs_sep, ghost_legacy_sep, nedlvs_legacy_sep) if s is not None]
    if not available_seps:
        n_ok_methods = sum([g_ok, n_ok, l_ok])
        if n_ok_methods == 0:
            host_agreement = "all_failed"
        elif n_ok_methods == 1:
            host_agreement = "only_one_method_succeeded"
        else:
            host_agreement = "unknown"
    elif all(s <= AGREEMENT_THRESHOLD_ARCSEC for s in available_seps):
        host_agreement = "agree"
    elif any(s <= AGREEMENT_THRESHOLD_ARCSEC for s in available_seps):
        host_agreement = "partial_agree"
    else:
        host_agreement = "disagree"

    return {
        "host_agreement":           host_agreement,
        "ghost_nedlvs_sep_arcsec":  round(ghost_nedlvs_sep,  3) if ghost_nedlvs_sep  is not None else None,
        "ghost_legacy_sep_arcsec":  round(ghost_legacy_sep,  3) if ghost_legacy_sep  is not None else None,
        "nedlvs_legacy_sep_arcsec": round(nedlvs_legacy_sep, 3) if nedlvs_legacy_sep is not None else None,
    }


def _distance_modulus_from_distmpc(dist_mpc: float) -> float:
    return 5.0 * (math.log10(dist_mpc * 1e6) - 1.0)


def _distance_modulus_from_z(z: float) -> float:
    dl_pc = COSMO.luminosity_distance(z).to(u.pc).value
    return 5.0 * (math.log10(dl_pc) - 1.0)


def _sep_kpc(sep_arcsec: float, z: float) -> float:
    d_a_kpc = COSMO.angular_diameter_distance(z).to(u.kpc).value
    return sep_arcsec * d_a_kpc / 206265.0


def _ned_redshift_near_position(
    ra_deg: float,
    dec_deg: float,
    radius_arcsec: float = 15.0,
) -> dict[str, Any]:
    """
    Look up the closest NED object with z > 0 near an already-selected host
    position. This is used when Legacy/Tractor gives the best morphology/host
    position, but we still want a NED redshift for DM.

    Default radius is 15 arcsec (not AGREEMENT_THRESHOLD_ARCSEC = 3 arcsec):
    Tractor photometric centroids can disagree with NED catalog centers by
    several arcsec for extended or irregular sources, so 3 arcsec is too tight
    and would systematically miss valid redshifts.
    """
    try:
        from astroquery.ipac.ned import Ned

        coord = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg)
        result = Ned.query_region(coord, radius=radius_arcsec * u.arcsec)
        if result is None or len(result) == 0:
            return {}

        df = result.to_pandas()
        df.columns = [str(c).strip().lower().replace(" ", "_") for c in df.columns]

        col_map: dict[str, str] = {}
        for col in df.columns:
            lc = col.lower()
            if (lc == "ra" or lc == "ra(deg)" or lc.endswith("_ra")) and "ra" not in col_map.values():
                col_map[col] = "ra"
            elif (lc == "dec" or lc == "dec(deg)" or lc.endswith("_dec")) and "dec" not in col_map.values():
                col_map[col] = "dec"
            elif "object_name" in lc and "objname" not in col_map.values():
                col_map[col] = "objname"
            elif lc in ("type", "object_type") and "objtype" not in col_map.values():
                col_map[col] = "objtype"
            elif lc == "redshift" and "z" not in col_map.values():
                col_map[col] = "z"
            elif lc == "velocity" and "z_raw_velocity" not in col_map.values():
                col_map[col] = "z_raw_velocity"
        df = df.rename(columns=col_map)

        # NED sometimes returns cz (km/s) in a "velocity" column instead of
        # a dimensionless "redshift" column.  Convert to z if needed.
        _C_KMS = 299792.458
        if "z_raw_velocity" in df.columns and "z" not in df.columns:
            df["z"] = df["z_raw_velocity"].apply(
                lambda v: _as_float(v) / _C_KMS if _as_float(v) is not None else None
            )

        if "ra" not in df.columns or "dec" not in df.columns or "z" not in df.columns:
            return {}

        for col in df.columns:
            if df[col].dtype == object:
                df[col] = df[col].apply(_clean_scalar)

        df["z"] = df["z"].apply(_as_float)
        df = df.loc[df["z"].notna() & (df["z"] > 0)].copy()

        if "objtype" in df.columns:
            df = df.loc[
                ~_non_host_mask(df["objtype"])
            ].copy()

        if len(df) == 0:
            return {}

        center = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg)
        coords = SkyCoord(
            ra=df["ra"].to_numpy(dtype=float) * u.deg,
            dec=df["dec"].to_numpy(dtype=float) * u.deg,
        )
        df["sep_arcsec"] = center.separation(coords).arcsec
        best = df.sort_values("sep_arcsec", kind="stable").iloc[0]

        return {
            "z": _as_float(best.get("z")),
            "sep_arcsec": _as_float(best.get("sep_arcsec")),
            "objname": _clean_scalar(best.get("objname")),
            "objtype": _clean_scalar(best.get("objtype")),
        }

    except Exception as exc:
        print(f"    [NED z near selected host] failed: {exc}")
        return {}


def _ned_agrees_with_legacy(
    nedlvs: dict[str, Any],
    legacy: dict[str, Any],
    threshold_arcsec: float = AGREEMENT_THRESHOLD_ARCSEC,
) -> bool:
    if nedlvs.get("nedlvs_status") != "host_found":
        return False
    if legacy.get("legacy_status") != "host_found":
        return False

    sep = _angular_sep_arcsec(
        nedlvs.get("nedlvs_ra"), nedlvs.get("nedlvs_dec"),
        legacy.get("legacy_ra"), legacy.get("legacy_dec"),
    )
    return bool(sep is not None and sep <= threshold_arcsec)


def _choose_host_for_boom(
    ghost: dict[str, Any],
    nedlvs: dict[str, Any],
    legacy: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """
    Choose one host for BOOM absolute-magnitude columns.

    New practical priority:
      1. Legacy/Tractor DLR for host position when available.
      2. NED DLR for redshift only if it agrees with Legacy, otherwise a small
         NED cone search around the Legacy host position is used for z/DM.
      3. NED DLR host if Legacy failed.
      4. GHOST only if both DLR methods failed.

    GHOST still always runs upstream and all GHOST columns are kept.
    """
    legacy = legacy or {}

    # ------------------------------------------------------------
    # 1) Prefer Legacy/Tractor for the chosen host position.
    # ------------------------------------------------------------
    if legacy.get("legacy_status") == "host_found":
        sep = _as_float(legacy.get("legacy_sep_arcsec"))
        z = dist_mpc = None
        dm = dm_source = None
        z_tech = None

        # Use NED DLR redshift only if it is the same host as Legacy.
        if _ned_agrees_with_legacy(nedlvs, legacy):
            z = _as_float(nedlvs.get("nedlvs_z"))
            dist_mpc = _as_float(nedlvs.get("nedlvs_distmpc"))
            if dist_mpc is not None and dist_mpc > 0:
                dm = _distance_modulus_from_distmpc(dist_mpc)
                dm_source = "ned_dlr_agrees_with_legacy_distmpc"
                z_tech = "NED_SPEC"
            elif z is not None and z > 0:
                dm = _distance_modulus_from_z(z)
                dm_source = "ned_dlr_agrees_with_legacy_z_cosmo"
                z_tech = "NED_SPEC"

        # If the NED-DLR ellipse chose a different huge/offset object, still
        # try to get a redshift at the Legacy host position.
        if dm is None:
            z_lookup = _ned_redshift_near_position(
                _as_float(legacy.get("legacy_ra")),
                _as_float(legacy.get("legacy_dec")),
            ) if legacy.get("legacy_ra") is not None and legacy.get("legacy_dec") is not None else {}
            z = _as_float(z_lookup.get("z"))
            if z is not None and z > 0:
                dm = _distance_modulus_from_z(z)
                dm_source = "ned_near_legacy_z_cosmo"
                z_tech = "NED_SPEC"

        return {
            "host_source":     "LEGACY_DLR",
            "host_z":          z,
            "host_sep_arcsec": sep,
            "host_sep_kpc":    _sep_kpc(sep, z) if sep is not None and z is not None and z > 0 else None,
            "dm":              dm,
            "dm_source":       dm_source,
            "z_tech":          z_tech,
            "z_unc":           None,
            "z_qual":          None,
            "has_nearby_host": bool(sep is not None),
        }

    # ------------------------------------------------------------
    # 2) If Legacy failed, use the NED DLR host.
    # ------------------------------------------------------------
    if nedlvs.get("nedlvs_status") == "host_found":
        z        = _as_float(nedlvs.get("nedlvs_z"))
        dist_mpc = _as_float(nedlvs.get("nedlvs_distmpc"))
        sep      = _as_float(nedlvs.get("nedlvs_sep_arcsec"))
        dm = dm_source = None
        if dist_mpc is not None and dist_mpc > 0:
            dm        = _distance_modulus_from_distmpc(dist_mpc)
            dm_source = "ned_dlr_distmpc"
        elif z is not None and z > 0:
            dm        = _distance_modulus_from_z(z)
            dm_source = "ned_dlr_z_cosmo"
        return {
            "host_source":     "NED_DLR",
            "host_z":          z,
            "host_sep_arcsec": sep,
            "host_sep_kpc":    _sep_kpc(sep, z) if sep is not None and z is not None and z > 0 else None,
            "dm":              dm,
            "dm_source":       dm_source,
            "z_tech":          "NED_SPEC" if z is not None else None,
            "z_unc":           None,
            "z_qual":          None,
            "has_nearby_host": bool(sep is not None),
        }

    # ------------------------------------------------------------
    # 3) Last fallback: GHOST. It still runs and is recorded either way.
    # ------------------------------------------------------------
    if ghost.get("ghost_status") == "host_found":
        z   = _as_float(ghost.get("ghost_z"))
        sep = _as_float(ghost.get("ghost_sep_arcsec"))
        dm  = _distance_modulus_from_z(z) if z is not None and z > 0 else None
        return {
            "host_source":     "GHOST",
            "host_z":          z,
            "host_sep_arcsec": sep,
            "host_sep_kpc":    _sep_kpc(sep, z) if sep is not None and z is not None and z > 0 else None,
            "dm":              dm,
            "dm_source":       "ghost_z_cosmo" if dm is not None else None,
            "z_tech":          "GHOST_PHOTO" if z is not None else None,
            "z_unc":           None,
            "z_qual":          None,
            "has_nearby_host": bool(sep is not None),
        }

    return {
        "host_source": None, "host_z": None, "host_sep_arcsec": None,
        "host_sep_kpc": None, "dm": None, "dm_source": None,
        "z_tech": None, "z_unc": None, "z_qual": None, "has_nearby_host": False,
    }


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------


def associate_host(
    transient_name: str,
    ra_deg: float,
    dec_deg: float,
    ghost_savepath: str = "/tmp/ghost_out/",
    run_legacy: bool = True,
    obj_idx: int = 0,
    obj_total: int = 0,
) -> dict[str, Any]:
    """
    Run GHOST + NED full DLR (+ optionally Legacy DLR) for one transient.

    Parameters
    ----------
    run_legacy : bool
        If True, also run the Legacy TAP DLR query.
        Requires astro-datalab installed. Set False to skip.
    obj_idx : int
        1-based index of this object in the processing order.
    obj_total : int
        Total number of objects to process.
    """
    counter = f" [{obj_idx}/{obj_total}]" if obj_total else ""
    print(f"  [HOST]{counter} {transient_name}  ra={ra_deg:.5f}  dec={dec_deg:.5f}", flush=True)

    t0 = time.time()

    ghost  = _run_ghost(transient_name, ra_deg, dec_deg, ghost_savepath=ghost_savepath)
    nedlvs = _run_ned_dlr(transient_name, ra_deg, dec_deg)
    legacy = _run_legacy_dlr(transient_name, ra_deg, dec_deg) if run_legacy else {"legacy_status": "skipped"}

    agreement = _check_agreement(ghost, nedlvs, legacy)
    chosen    = _choose_host_for_boom(ghost, nedlvs, legacy)

    elapsed = time.time() - t0
    print(f"    GHOST:   {ghost.get('ghost_status')}  sep={ghost.get('ghost_sep_arcsec')}  z={ghost.get('ghost_z')}", flush=True)
    print(f"    NED DLR: {nedlvs.get('nedlvs_status')}  dlr={nedlvs.get('nedlvs_dlr')}  sep={nedlvs.get('nedlvs_sep_arcsec')}  z={nedlvs.get('nedlvs_z')}", flush=True)
    print(f"    Legacy:  {legacy.get('legacy_status')}  dlr={legacy.get('legacy_dlr')}  sep={legacy.get('legacy_sep_arcsec')}", flush=True)
    print(f"    Agreement: {agreement.get('host_agreement')}  GN={agreement.get('ghost_nedlvs_sep_arcsec')}\"  GL={agreement.get('ghost_legacy_sep_arcsec')}\"  NL={agreement.get('nedlvs_legacy_sep_arcsec')}\"  Chosen: {chosen.get('host_source')}  ({elapsed:.1f}s)", flush=True)

    return {**ghost, **nedlvs, **legacy, **agreement, **chosen}