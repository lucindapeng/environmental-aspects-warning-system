import ee

def get_environmental_risk(aoi, date_str):
    """
    Calculates the 90-day conflict risk based on 6 environmental factors.
    """
    # 1. Set Time Windows
    now = ee.Date(date_str)
    three_months_ago = now.advance(-3, 'month')
    baseline_start = '2005-01-01'
    baseline_end = '2025-01-01'

    # 2. Load Collections
    ndvi = ee.ImageCollection("MODIS/061/MOD13Q1").select('NDVI')
    chirps = ee.ImageCollection("UCSB-CHG/CHIRPS/PENTAD")
    s1 = ee.ImageCollection("COPERNICUS/S1_GRD")\
        .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))\
        .select('VV')
    fire = ee.ImageCollection("FIRMS").select('T21')
    burned = ee.ImageCollection("MODIS/061/MCD64A1").select('BurnDate')

    # 3. CALCULATE ROOT CAUSES
    veg_baseline_mean = ndvi.filterDate(baseline_start, baseline_end).mean()
    veg_baseline_std = ndvi.filterDate(baseline_start, baseline_end).reduce(ee.Reducer.stdDev())
    veg_current = ndvi.filterDate(three_months_ago, now).mean()
    z_score = veg_current.subtract(veg_baseline_mean).divide(veg_baseline_std)

    rain_baseline = chirps.filterDate(baseline_start, baseline_end).sum().divide(20)
    rain_current = chirps.filterDate(three_months_ago, now).sum()
    rain_ratio = rain_current.divide(rain_baseline.divide(4)) 

    # 4. CALCULATE INDICATORS
    bare_ground = s1.filterDate(three_months_ago, now).mean()
    recent_fire = fire.filterDate(now.advance(-1, 'month'), now).max()
    burn_scar = burned.filterDate(now.advance(-6, 'month'), now).max()

    # 5. THE INTEGRATED RED ZONE LOGIC
    stress_condition = z_score.lt(-1.5).or(rain_ratio.lt(0.7))
    evidence_condition = bare_ground.gt(-10).or(burn_scar.gt(0)).or(recent_fire.gt(310))

    final_red_zone = stress_condition.and(evidence_condition)

    return final_red_zone.clip(aoi).rename('environmental_risk')
