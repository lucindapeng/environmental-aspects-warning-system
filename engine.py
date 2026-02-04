import ee

def get_environmental_risk(aoi, date_str):
    # 1. Setup Dates
    now = ee.Date(date_str)
    three_months_ago = now.advance(-3, 'month')
    baseline_start = '2005-01-01'
    baseline_end = '2025-01-01'

    # 2. Load Multi-Sensor Collections
    ndvi = ee.ImageCollection("MODIS/061/MOD13Q1").select('NDVI')
    chirps = ee.ImageCollection("UCSB-CHG/CHIRPS/PENTAD")
    s1 = ee.ImageCollection("COPERNICUS/S1_GRD").filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV')).select('VV')
    fire = ee.ImageCollection("FIRMS").select('T21')
    burned = ee.ImageCollection("MODIS/061/MCD64A1").select('BurnDate')

    # 3. Calculate Primary Drivers (Causes)
    veg_mean = ndvi.filterDate(baseline_start, baseline_end).mean()
    veg_std = ndvi.filterDate(baseline_start, baseline_end).reduce(ee.Reducer.stdDev())
    z_score = ndvi.filterDate(three_months_ago, now).mean().subtract(veg_mean).divide(veg_std)

    rain_baseline = chirps.filterDate(baseline_start, baseline_end).sum().divide(20)
    rain_current = chirps.filterDate(three_months_ago, now).sum()
    rain_ratio = rain_current.divide(rain_baseline.divide(4))

    # 4. Calculate Physical Indicators (Results)
    bare_ground = s1.filterDate(three_months_ago, now).mean()
    recent_fire = fire.filterDate(now.advance(-1, 'month'), now).max()
    burn_scar = burned.filterDate(now.advance(-6, 'month'), now).max()

    # 5. Combined Threshold Logic
    stress = z_score.lt(-1.5).or(rain_ratio.lt(0.7))
    evidence = bare_ground.gt(-10).or(burn_scar.gt(0)).or(recent_fire.gt(310))
    
    return stress.and(evidence).clip(aoi).rename('risk_score')
