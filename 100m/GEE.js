var mask = ee.Image('users/arthfen/Mask');
var bbx = mask.geometry().bounds();

var imgVs = ee.ImageCollection('COPERNICUS/S1_GRD').filterBounds(bbx)
        .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
        .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
        .filter(ee.Filter.eq('instrumentMode', 'IW'))
        .select(['VV', 'VH'])
        .map(function(image) {
          var edge = image.lt(-30.0);
          var maskedImage = image.mask().and(edge.not());
          return image.updateMask(maskedImage);
        });
var proj = imgVs.first().projection();
// First, 2015

var start_date = ee.Date('2016-08-10');
var end_date = 0;

var i = 0;
while((ee.Date('2017-08-10').difference(start_date, 'day')).getInfo() > 0){ // I use 2017-08-01 so the last and the first one are in the same day, approximately
  i = i + 1;
  end_date = start_date.advance(12, 'day');
//  print(start_date.format('yyyyMMdd'));

  var mo = ee.Filter.date(start_date, end_date);
  var out = ee.Image.cat(imgVs.filter(mo).median());
  var outVV = out.select('VV');
  var outVH = out.select('VH');
  var outCR = outVH.subtract(outVV).rename('CR');
  if(i == 1) Map.addLayer(outCR);
  
  Export.image.toDrive({  
   image: outCR.multiply(1000).int(),
   folder: 'JRC',
   description: 'IndexesSAll-' + i + '-' + start_date.format('yyyyMMdd').getInfo(),
   scale: 100,
   region: bbx,
   crs: 'EPSG:3035',
   fileFormat: 'GeoTIFF',
   skipEmptyTiles: true,
   maxPixels: 1e13
  });

  start_date = end_date;
}
