import ee

def get_environmental_risk(aoi, date_str):
    # 1. Setup Dates
    now = ee.Date(date_str)
    three_months_ago = now.advance(-3, 'month')
    baseline_start = '2005-01-01'
    baseline_end = '2024-01-01'

    # 2. Load Collections
    ndvi = ee.ImageCollection("MODIS/061/MOD13Q1").select('NDVI')
    chirps = ee.ImageCollection("UCSB-CHG/CHIRPS/PENTAD")
    s1 = ee.ImageCollection("COPERNICUS/S1_GRD").filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV')).select('VV')
    fire = ee.ImageCollection("FIRMS").select('T21')

    # 3. Calculate Percentile Flags (Current vs. 20-Year History)
    
    # --- VEGETATION (Worst 20%) ---
    veg_hist = ndvi.filterDate(baseline_start, baseline_end)
    veg_current = ndvi.filterDate(three_months_ago, now).mean()
    veg_p20 = veg_hist.reduce(ee.Reducer.percentile([20]))
    veg_flag = veg_current.lt(veg_p20)

    # --- RAINFALL (Worst 25%) ---
    rain_hist = chirps.filterDate(baseline_start, baseline_end).sum().divide(19)
    rain_current = chirps.filterDate(three_months_ago, now).sum()
    rain_p25 = rain_hist.multiply(0.75) # Using 75% of mean as p25 proxy for speed
    rain_flag = rain_current.lt(rain_p25)

    # --- BARE GROUND (Top 30% VV brightness) ---
    vv_hist = s1.filterDate(baseline_start, baseline_end).reduce(ee.Reducer.percentile([70]))
    vv_current = s1.filterDate(three_months_ago, now).mean()
    bare_flag = vv_current.gt(vv_hist)

    # --- FIRE (Binary detection) ---
    fire_current = fire.filterDate(three_months_ago, now).max()
    fire_flag = fire_current.gt(310) # Simple thermal threshold for fire existence

    # 4. STEP 3: GRADED STRESS SCORE (Additive)
    # Sum the flags: 0 (None) to 4 (All indicators extreme)
    stress_score = veg_flag.add(rain_flag).add(bare_flag).add(fire_flag)

    return stress_score.clip(aoi).rename('stress_score')
