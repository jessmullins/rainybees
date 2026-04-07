// ============================================================================
// evi_extraction.js
//
// Purpose : Extract MODIS MOD13Q1 Enhanced Vegetation Index (EVI) for
//           natural wildland vegetation in western San Diego County during
//           the growing season (March–August), 2000–2024.
//           Applies a custom wildlands mask (uploaded GeoTIFF asset) to
//           restrict analysis to natural vegetation only.
//
// Inputs  : MODIS/061/MOD13Q1                          — 16-day EVI composites, 250 m
//           users/jessmullins/Wildlands_mask_sitesBBox  — Custom wildlands mask asset
//               Asset URL: https://code.earthengine.google.com/?asset=users/jessmullins/Wildlands_mask_sitesBBox
//
// Output  : WildlandsOnly_EVI_MarchAugust_2000_2024.csv
//           — One row per year; columns: year, mean_EVI
//           — Export to Google Drive, then move to data/raw/wild_evi.csv
//
// Usage   : Open in the Google Earth Engine Code Editor at
//           https://code.earthengine.google.com/
//           The wildlands mask asset is public. Click Run, then go to the
//           Tasks tab and click Run again next to the export task.
//
// Author  : Jessica Mullins
// Date    : 2026
// ============================================================================

// --- 1) Load custom wildlands mask (GeoTIFF uploaded as EE asset) -----------
var wildlands = ee.Image("users/jessmullins/Wildlands_mask_sitesBBox");

// --- 2) Build list of growing seasons (2000–2024) ---------------------------
var years = ee.List.sequence(2000, 2024);

var growingSeasonYears = years.map(function(y) {
  y = ee.Number(y);
  return {
    'year' : y,
    'start': ee.Date.fromYMD(y, 3, 1),  // March 1
    'end'  : ee.Date.fromYMD(y, 8, 31)  // August 31
  };
});

// --- 3) Load MODIS EVI collection -------------------------------------------
// MOD13Q1: 16-day composite, 250 m resolution, Version 6.1
var modisEVI = ee.ImageCollection("MODIS/061/MOD13Q1")
                 .select('EVI');

// --- 4) Compute wildlands-masked mean EVI per growing season ----------------
var eviResults = growingSeasonYears.map(function(d) {
  d = ee.Dictionary(d);

  var filtered = modisEVI
    .filterDate(d.get('start'), d.get('end'))
    .filterBounds(wildlands.geometry());

  var count = filtered.size();
  var valid = count.gt(0);

  return ee.Algorithms.If(
    valid,
    ee.Feature(null, {
      'year'    : d.get('year'),
      'mean_EVI': filtered.mean()
          .multiply(0.0001)          // apply MODIS scale factor (stored as ×10000)
          .updateMask(wildlands)     // restrict to wildland pixels only
          .reduceRegion({
            reducer  : ee.Reducer.mean(),
            geometry : wildlands.geometry(),
            scale    : 250,
            maxPixels: 1e13
          }).get('EVI')
    }),
    // If no images found for that year, return null
    ee.Feature(null, {
      'year'    : d.get('year'),
      'mean_EVI': null
    })
  );
});

// --- 5) Preview in console --------------------------------------------------
var eviTable = ee.FeatureCollection(eviResults);
print('Wildland-only March–August mean EVI 2000–2024:', eviTable);

// --- 6) Export to Google Drive ----------------------------------------------
// Go to the Tasks tab and click Run next to the export task.
// Download the CSV and save it as data/raw/wild_evi.csv in this repository.
Export.table.toDrive({
  collection : eviTable,
  description: 'WildlandsOnly_EVI_MarchAugust_2000_2024',
  fileFormat : 'CSV'
});
